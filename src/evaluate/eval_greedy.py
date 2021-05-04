import os
import numpy as np

import time
from datetime import datetime

from rdkit import Chem
from rdkit.Chem.Descriptors import MolLogP

import torch
from torch_geometric.data import Batch

from utils.graph_utils import mol_to_pyg_graph

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

from logp.get_reward import get_logp_scores


def get_rewards(g_batch, surrogate_model):
    with torch.autograd.no_grad():
        scores = surrogate_model(g_batch.to(DEVICE))
    return scores.cpu().numpy()*-1

def greedy_rollout(save_path, env, reward_type, surrogate_guide, surrogate_eval, K, max_rollout=6):
    mol, mol_candidates, done = env.reset()
    mol_start = mol
    mol_best = mol

    if reward_type == 'surr':
        g = Batch().from_data_list([mol_to_pyg_graph(mol)[0]]).to(DEVICE)
        new_rew = get_rewards(g, surrogate_guide)
    elif reward_type == 'logp':
        new_rew = get_logp_rewards(mol)

    start_rew = new_rew
    best_rew = new_rew
    steps_remaining = K

    for i in range(max_rollout):
        print(f"  {i+1} {steps_remaining} {new_rew}")
        steps_remaining -= 1
        
        if reward_type == 'surr':
            g_candidates = Batch().from_data_list([mol_to_pyg_graph(cand)[0] for cand in mol_candidates]).to(DEVICE)
            next_rewards = get_rewards(g_candidates, surrogate_guide)
        elif reward_type == 'logp':
            next_rewards = get_logp_rewards(mol_candidates)
        
        action = np.argmax(next_rewards)

        try:
            new_rew = next_rewards[action]
        except Exception as e:
            print(e)
            break
        
        mol, mol_candidates, done = env.step(action)
        
        if reward_type == 'surr':
            g = Batch().from_data_list([mol_to_pyg_graph(mol)[0]]).to(DEVICE)

        if new_rew > best_rew:
            mol_best = mol
            best_rew = new_rew
            steps_remaining = K

        if (steps_remaining == 0) or done:
            break

    with open(save_path, 'a') as f:
        print("Writing SMILE molecules!")

        smile = Chem.MolToSmiles(mol_best, isomericSmiles=False)
        print(smile, new_rew)
        row = ''.join(['{},'] * 2)[:-1] + '\n'
        f.write(row.format(smile, new_rew))

    return start_rew, best_rew

def eval_greedy(artifact_path, reward_type, surrogate_guide, surrogate_eval, env, N=30, K=1):
    # logging variables
    dt = datetime.now().strftime("%Y.%m.%d_%H:%M:%S")
    save_path = os.path.join(artifact_path, dt + '_greedy.csv')
    
    if reward_type == 'surr':
        surrogate_guide = surrogate_guide.to(DEVICE)
        surrogate_guide.eval()
        surrogate_eval  = surrogate_eval.to(DEVICE)
        surrogate_eval.eval()

    print("\nStarting greedy...\n")
    avg_improvement = []
    avg_best = []
    for i in range(N):
        start_rew, best_rew = greedy_rollout(save_path, env, reward_type, surrogate_guide, surrogate_eval, K)
        improvement = best_rew - start_rew
        print(f"{i+1}: {start_rew} {best_rew} {improvement}\n")
        avg_improvement.append(improvement)
        avg_best.append(best_rew)
    avg_improvement = sum(avg_improvement) / len(avg_improvement)
    avg_best = sum(avg_best) / len(avg_best)
    print(f"Avg improvement over {N} samples: {avg_improvement}")
    print(f"Avg best        over {N} samples: {avg_best}")
