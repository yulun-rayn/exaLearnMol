import os
import shutil
import subprocess
import sys
import gym
import copy
import numpy as np
from collections import deque, OrderedDict

import torch
import torch.multiprocessing as mp
import torch.nn as nn
from torch.distributions import MultivariateNormal
from torch.utils.tensorboard import SummaryWriter

from torch_geometric.data import Data, Batch
from torch_geometric.utils import dense_to_sparse

from .gcpn_policy import GCPN_CReM

from utils.general_utils import get_current_datetime
from utils.graph_utils import mol_to_pyg_graph, get_batch_shift

from gnn_surrogate.model import GNN_MyGAT

from rdkit import Chem

class Memory:
    def __init__(self):
        self.actions = []
        self.states = []
        self.logprobs = []
        self.rewards = []
        self.is_terminals = []

    def extend(self, memory):
        self.actions.extend(memory.actions)
        self.states.extend(memory.states)
        self.logprobs.extend(memory.logprobs)
        self.rewards.extend(memory.rewards)
        self.is_terminals.extend(memory.is_terminals)

    def clear(self):
        del self.actions[:]
        del self.states[:]
        del self.logprobs[:]
        del self.rewards[:]
        del self.is_terminals[:]


#################################################
#                   GCPN PPO                    #
#################################################

class GCPN_Critic(nn.Module):
    def __init__(self, emb_dim, nb_layers, nb_hidden):
        super(GCPN_Critic, self).__init__()
        layers = [nn.Linear(emb_dim, nb_hidden)]
        for _ in range(nb_layers-1):
            layers.append(nn.Linear(nb_hidden, nb_hidden))

        self.layers = nn.ModuleList(layers)
        self.final_layer = nn.Linear(nb_hidden, 1)
        self.act = nn.ReLU()

    def forward(self, X):
        for i, l in enumerate(self.layers):
            X = self.act(l(X))
        return self.final_layer(X).squeeze(1)


class ActorCriticGCPN(nn.Module):
    def __init__(self,
                 input_dim,
                 emb_dim,
                 nb_edge_types,
                 gnn_nb_layers,
                 gnn_nb_hidden,
                 mlp_nb_layers,
                 mlp_nb_hidden):
        super(ActorCriticGCPN, self).__init__()

        # action mean range -1 to 1
        self.actor = GCPN_CReM(input_dim,
                               emb_dim,
                               nb_edge_types,
                               gnn_nb_layers,
                               gnn_nb_hidden,
                               mlp_nb_layers,
                               mlp_nb_hidden)
        # critic
        self.critic = GCPN_Critic(emb_dim, mlp_nb_layers, mlp_nb_hidden)

    def forward(self):
        raise NotImplementedError

    def select_action(self, states, candidates, surrogate_model, batch_idx):
        return self.actor.select_action(states, candidates, surrogate_model, batch_idx)

    def evaluate(self, states, candidates, actions):   
        probs = self.actor.evaluate(candidates, actions)

        action_logprobs = torch.log(probs)
        state_value = self.critic(states)

        entropy = probs * action_logprobs

        return action_logprobs, state_value, entropy


def wrap_state(ob):
    adj = ob['adj']
    nodes = ob['node'].squeeze()

    adj = torch.Tensor(adj)
    nodes = torch.Tensor(nodes)

    adj = [dense_to_sparse(a) for a in adj]
    data = Data(x=nodes, edge_index=adj[0][0], edge_attr=adj[0][1])
    return data


class PPO_GCPN(nn.Module):
    def __init__(self,
                 lr,
                 betas,
                 eps,
                 gamma,
                 eta,
                 upsilon,
                 K_epochs,
                 eps_clip,
                 input_dim,
                 emb_dim,
                 nb_edge_types,
                 gnn_nb_layers,
                 gnn_nb_hidden,
                 mlp_nb_layers,
                 mlp_nb_hidden):
        super(PPO_GCPN, self).__init__()
        self.lr = lr
        self.betas = betas
        self.eps = eps
        self.gamma = gamma
        self.eta = eta
        self.upsilon = upsilon
        self.eps_clip = eps_clip
        self.K_epochs = K_epochs
        
        self.policy = ActorCriticGCPN(input_dim,
                                      emb_dim,
                                      nb_edge_types,
                                      gnn_nb_layers,
                                      gnn_nb_hidden,
                                      mlp_nb_layers,
                                      mlp_nb_hidden)
        self.optimizer = torch.optim.Adam(self.policy.parameters(), lr=lr, betas=betas, eps=eps)
        
        self.policy_old = ActorCriticGCPN(input_dim,
                                          emb_dim,
                                          nb_edge_types,
                                          gnn_nb_layers,
                                          gnn_nb_hidden,
                                          mlp_nb_layers,
                                          mlp_nb_hidden)
        self.policy_old.load_state_dict(self.policy.state_dict())
        
        self.MseLoss = nn.MSELoss()

        self.device = torch.device("cpu")

    def to_device(self, device):
        self.policy.to(device)
        self.policy_old.to(device)
        self.device = device

    def forward(self):
        raise NotImplementedError

    def select_action(self, states, candidates, surrogate_model, batch_idx=None, return_shifted=False):
        if not isinstance(states, list):
            states = [states]
        if batch_idx is None:
            batch_idx = torch.empty(len(candidates), dtype=torch.long).fill_(0)
        batch_idx = batch_idx.to(self.device)

        g = Batch.from_data_list(states).to(self.device)
        g_candidates = Batch.from_data_list(candidates).to(self.device)
        with torch.autograd.no_grad():
            g_emb, X_states, action_logprobs, actions, shifted_actions = self.policy_old.select_action(g, g_candidates, surrogate_model, batch_idx)

        if return_shifted:
            return [g_emb, X_states], action_logprobs, actions, shifted_actions
        else:
            return [g_emb, X_states], action_logprobs, actions

    def update(self, memory, save_dir):
        # Monte Carlo estimate of rewards:
        rewards = []
        discounted_reward = 0
        for reward, is_terminal in zip(reversed(memory.rewards), reversed(memory.is_terminals)):
            if is_terminal:
                discounted_reward = 0
            discounted_reward = reward + (self.gamma * discounted_reward)
            rewards.insert(0, discounted_reward)

        # Normalizing the rewards:
        rewards = torch.tensor(rewards).to(self.device)
        rewards = (rewards - rewards.mean()) / (rewards.std() + 1e-5)

        # convert list to tensor
        old_states = torch.cat(([m[0] for m in memory.states]),dim=0).to(self.device)
        old_candidates = Batch().from_data_list([Data(x=m[1]) for m in memory.states]).to(self.device)
        old_actions = torch.tensor(memory.actions).to(self.device)
        old_logprobs = torch.tensor(memory.logprobs).to(self.device)

        # Optimize policy for K epochs:
        print("Optimizing...")

        for i in range(self.K_epochs):
            # Evaluating old actions and values :
            logprobs, state_values, entropies = self.policy.evaluate(old_states,
                                                                     old_candidates,
                                                                     old_actions)

            # Finding the ratio (pi_theta / pi_theta__old):
            ratios = torch.exp(logprobs - old_logprobs)
            advantages = rewards - state_values.detach()

            # loss
            ## policy
            surr1 = ratios * advantages
            surr2 = torch.clamp(ratios, 1-self.eps_clip, 1+self.eps_clip) * advantages
            loss = -torch.min(surr1, surr2)
            if torch.isnan(loss).any():
                print("found nan in loss")
                exit()
            ## entropy
            loss += self.eta * entropies
            ## baseline
            loss = loss.mean() + self.upsilon*self.MseLoss(state_values, rewards)

            ## take gradient step
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
            if (i%10)==0:
                print("  {:3d}: Loss: {:7.3f}".format(i, loss))

        # save running model
        torch.save(self.policy.state_dict(), os.path.join(save_dir, 'running_gcpn.pth'))
        # Copy new weights into old policy:
        self.policy_old.load_state_dict(self.policy.state_dict())

    def __repr__(self):
        return "{}\n{}".format(repr(self.policy), repr(self.optimizer))


#####################################################
#                   FINAL REWARDS                   #
#####################################################

def get_reward(states, surrogate_model, device, done_idx=None):
    if not isinstance(states, list):
        states = [states]
    if done_idx is None:
        done_idx = range(len(states))

    #print(states)
    smiles = [Chem.MolToSmiles(mol) for mol in states]
    #print(smiles)

    #Print smiles to a file as input to ADT-gpu
    adttmp="./src/adtgpu/autodockgpu"
    if not os.path.exists(adttmp):
        os.makedirs(adttmp)
    if not os.path.exists(adttmp+"/ligands"):
        os.makedirs(adttmp+"/ligands")
    with open(adttmp+"/smiles.txt","w") as smiles_file:
        smiles_file.writelines("%s\n" % sm for sm in smiles)

    #Run ADT-gpu
    #cmd="python /gpfs/alpine/syb105/proj-shared/Personal/mcashman/Projects/GCPN_MLDrugDiscovery/exaLearnMol_crem_parallel_GPU/scripts/run_adtgpu.py -r /gpfs/alpine/syb105/proj-shared/Projects/StructPred/Software/progs/adtgpu_sa/gcpn/scripts/test_run.Ver4/receptor/NSP15_6W01_A_1_F_receptor.pdb -s \"" + str(smiles) + "\" -d autodockgpu"    
    cmd="python ./src/adtgpu/run_adtgpu.py -r ./src/adtgpu/receptor/NSP15_6W01_A_1_F_receptor.pdb -s \"" + str(smiles) + "\" -d "+adttmp    
    subprocess.Popen(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).wait()

    #Fetch output scores
    pred_docking_score=[]
    with open(adttmp+"/adt_affinity_vals","r") as scores:
        for mol in scores.readlines():
            pred_docking_score.append(-float(mol.strip()))
    #print("\npred_docking_score:\n{}".format(pred_docking_score))

    shutil.rmtree(adttmp)
    print("Reward Scores (-dock): {}".format(pred_docking_score))
    return (pred_docking_score)


#####################################################
#                   TRAINING LOOP                   #
#####################################################

################## Process ##################
tasks = mp.JoinableQueue()
results = mp.Queue()

# logging variables
running_reward = mp.Value("f", 0)
avg_length = mp.Value("i", 0)

class Worker(mp.Process):
    def __init__(self, env, task_queue, result_queue, max_timesteps):
        super(Worker, self).__init__()
        self.task_queue = task_queue
        self.result_queue = result_queue

        self.env = env

        self.max_timesteps = max_timesteps
        self.timestep_counter = 0

        np.random.seed(self.pid)

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task == None:
                # Poison pill means shutdown
                print('%s: Exiting' % proc_name)
                self.task_queue.task_done()
                break

            index, mol, done = next_task
            if index is None:
                self.result_queue.put((None, None, None, True))
                self.task_queue.task_done()
                continue
            #print('%s: Working' % proc_name)
            if done:
                self.timestep_counter = 0
                mol, candidates, done = self.env.reset()
            else:
                self.timestep_counter += 1
                mol, candidates, done = self.env.reset(mol)
                if self.timestep_counter >= self.max_timesteps:
                    done = True

            self.result_queue.put((index, mol, candidates, done))
            self.task_queue.task_done()
        return
'''
class Task(object):
    def __init__(self, index, action, done):
        self.index = index
        self.action = action
        self.done = done
    def __call__(self):
        return (self.action, self.done)
    def __str__(self):
        return '%d' % self.index

class Result(object):
    def __init__(self, index, state, candidates, done):
        self.index = index
        self.state = state
        self.candidates = candidates
        self.done = done
    def __call__(self):
        return (self.state, self.candidates, self.done)
    def __str__(self):
        return '%d' % self.index
'''
#############################################

def train_ppo(args, surrogate_model, env):
    ############## Hyperparameters ##############
    render = True
    solved_reward = 100         # stop training if avg_reward > solved_reward
    log_interval = 20           # print avg reward in the interval
    save_interval = 100         # save model in the interval

    max_episodes = 50000        # max training episodes
    max_timesteps = 6           # max timesteps in one episode
    update_timesteps = 500      # update policy every n timesteps

    K_epochs = 80               # update policy for K epochs
    eps_clip = 0.2              # clip parameter for PPO
    gamma = 0.99                # discount factor

    lr = 0.0001                 # parameters for Adam optimizer
    betas = (0.9, 0.999)
    eps = 0.01
    print("lr:", lr, "beta:", betas, "eps:", eps)

    #############################################

    print('Creating %d processes' % args.nb_procs)
    workers = [Worker(env, tasks, results, max_timesteps) for i in range(args.nb_procs)]
    for w in workers:
        w.start()

    # logging variables
    dt = get_current_datetime()
    writer = SummaryWriter(log_dir=os.path.join(args.artifact_path, 'runs/' + args.name + '_' + dt))
    save_dir = os.path.join(args.artifact_path, 'saves/' + args.name + '_' + dt)
    os.makedirs(save_dir, exist_ok=True)

    device = torch.device("cpu") if args.use_cpu else torch.device(
        'cuda:' + str(args.gpu) if torch.cuda.is_available() else "cpu")

    ppo = PPO_GCPN(lr,
                   betas,
                   eps,
                   gamma,
                   args.eta,
                   args.upsilon,
                   K_epochs,
                   eps_clip,
                   args.input_size,
                   args.emb_size,
                   args.nb_edge_types,
                   args.layer_num_g,
                   args.num_hidden_g,
                   args.mlp_num_layer,
                   args.mlp_num_hidden)
    ppo.to_device(device)
    print(ppo)

    surrogate_model.to(device)
    surrogate_model.eval()
    print(surrogate_model)

    sample_count = 0
    episode_count = 0
    update_count = 0 # for adversarial
    save_counter = 0
    log_counter = 0

    avg_length = 0
    running_reward = 0

    memory = Memory()
    memories = [Memory() for _ in range(args.nb_procs)]
    rewbuffer_env = deque(maxlen=100)
    # training loop
    i_episode = 0
    while i_episode < max_episodes:
        print("collecting rollouts")
        for i in range(args.nb_procs):
            tasks.put((i, None, True))
        tasks.join()
        # unpack results
        mols = [None]*args.nb_procs
        done_idx = []
        notdone_idx, candidates, batch_idx = [], [], []
        for i in range(args.nb_procs):
            index, mol, cands, done = results.get()

            mols[index] = mol

            notdone_idx.append(index)
            candidates.extend(cands)
            batch_idx.extend([index]*len(cands))
        batch_idx = torch.LongTensor(batch_idx)
        while True:
            # action selections (for not done)
            if len(notdone_idx) > 0:
                state_embs, action_logprobs, actions, shifted_actions = ppo.select_action(
                    [mol_to_pyg_graph(mols[idx])[0] for idx in notdone_idx], 
                    [mol_to_pyg_graph(cand)[0] for cand in candidates], 
                    surrogate_model, batch_idx, return_shifted=True)
            else:
                if sample_count >= update_timesteps:
                    break

            for i, idx in enumerate(notdone_idx):
                tasks.put((idx, candidates[shifted_actions[i]], False))

                memories[idx].states.append([state_embs[0][[i], :], state_embs[1][batch_idx == idx, :]])
                memories[idx].actions.append(actions[i])
                memories[idx].logprobs.append(action_logprobs[i])
            for idx in done_idx:
                if sample_count >= update_timesteps:
                    tasks.put((None, None, True))
                else:
                    tasks.put((idx, None, True))
            tasks.join()
            # unpack results
            mols = [None]*args.nb_procs
            new_done_idx = []
            new_notdone_idx, candidates, batch_idx = [], [], []
            for i in range(args.nb_procs):
                index, mol, cands, done = results.get()

                if index is not None:
                    mols[index] = mol
                if done:
                    new_done_idx.append(index)
                else:
                    new_notdone_idx.append(index)
                    candidates.extend(cands)
                    batch_idx.extend([index]*len(cands))
            batch_idx = torch.LongTensor(batch_idx)
            # get rewards (for done)
            nowdone_idx = [idx for idx in notdone_idx if idx in new_done_idx]
            stillnotdone_idx = [idx for idx in notdone_idx if idx in new_notdone_idx]
            if len(nowdone_idx) > 0:
                surr_rewards = get_reward(
                    [mols[idx] for idx in nowdone_idx],
                    surrogate_model, device)

            for i, idx in enumerate(nowdone_idx):
                i_episode += 1
                episode_count += 1
                avg_length += 1
                running_reward += surr_rewards[i]
                writer.add_scalar("EpSurrogate", -1*surr_rewards[i], i_episode-1)
                rewbuffer_env.append(surr_rewards[i])
                writer.add_scalar("EpRewEnvMean", np.mean(rewbuffer_env), i_episode-1)

                memories[idx].rewards.append(surr_rewards[i])
                memories[idx].is_terminals.append(True)
            for idx in stillnotdone_idx:
                avg_length += 1
                running_reward += 0

                memories[idx].rewards.append(0)
                memories[idx].is_terminals.append(False)

            sample_count += len(notdone_idx)
            done_idx = new_done_idx
            notdone_idx = new_notdone_idx

        for m in memories:
            memory.extend(m)
            m.clear()
        # update model
        print("\n\nupdating ppo @ episode %d..." % i_episode)
        ppo.update(memory, save_dir)
        memory.clear()

        update_count += 1
        save_counter += episode_count
        log_counter += episode_count

        # stop training if avg_reward > solved_reward
        if np.mean(rewbuffer_env) > solved_reward:
            print("########## Solved! ##########")
            torch.save(ppo.policy.state_dict(), os.path.join(save_dir, 'PPO_continuous_solved_{}.pth'.format('test')))
            break

        # save every 500 episodes
        if save_counter > save_interval:
            torch.save(ppo.policy.state_dict(), os.path.join(save_dir, '{:05d}_gcpn.pth'.format(i_episode)))
            save_counter -= save_interval

        if log_counter > log_interval:
            avg_length = int(avg_length/log_counter)
            running_reward = running_reward/log_counter
            
            print('Episode {} \t Avg length: {} \t Avg reward: {:5.3f}'.format(i_episode, avg_length, running_reward))
            running_reward = 0
            avg_length = 0
            log_counter -= log_interval

        episode_count = 0
        sample_count = 0

    # Add a poison pill for each process
    for i in range(args.nb_procs):
        tasks.put(None)
    tasks.join()

