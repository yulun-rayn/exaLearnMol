from collections import deque
import numpy as np

import torch
import torch.nn as nn

import torch_geometric as pyg
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import remove_self_loops, add_self_loops, softmax, degree

from gnn_embed.model import MyGNN

def init_network(model, method='uniform'):
    if model is not None:
        if method == 'uniform':
            model.weight.data.uniform_()
            model.bias.data.uniform_()
        elif method == 'normal':
            model.weight.data.normal_()
            model.bias.data.normal_()
        else:
            pass


class RNDistillation(nn.Module):
    def __init__(self,
                 lr,
                 betas,
                 eps,
                 input_dim,
                 nb_edge_types,
                 use_3d,
                 gnn_nb_layers,
                 gnn_nb_hidden,
                 rnd_nb_layers,
                 rnd_nb_hidden,
                 rnd_nb_output):
        super(RNDistillation, self).__init__()

        self.f = RandomNetwork(input_dim,
                               nb_edge_types,
                               use_3d,
                               gnn_nb_layers,
                               gnn_nb_hidden,
                               rnd_nb_layers,
                               rnd_nb_hidden,
                               rnd_nb_output,
                               init_method='uniform')
        self.f_hat = RandomNetwork(input_dim,
                                   nb_edge_types,
                                   use_3d,
                                   gnn_nb_layers,
                                   gnn_nb_hidden,
                                   rnd_nb_layers,
                                   rnd_nb_hidden,
                                   rnd_nb_output,
                                   init_method='normal')

        self.optimizer = torch.optim.Adam(self.f_hat.parameters(), lr=lr, betas=betas, eps=eps)

        self.running_error = deque(maxlen=5000)

    def forward(self, g_next_emb):
        errors = torch.norm(self.f(g_next_emb).detach() - self.f_hat(g_next_emb), dim=1)
        return errors

    def get_score(self, g_next_emb, out_min=-5., out_max=5., min_running=100, eps=0.01):
        with torch.autograd.no_grad():
            errors = self(g_next_emb).detach().cpu().numpy()

        if len(self.running_error) < min_running:
            return np.zeros_like(errors)
        scores = (errors - np.mean(self.running_error)) / (np.std(self.running_error) + eps)
        return np.clip(scores, out_min, out_max)

    def update(self, g_next_emb):
        errors = self(g_next_emb)
        loss = errors.mean()

        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()

        self.running_error.extend(errors.detach().cpu().numpy())
        return loss.item()


class RandomNetwork(nn.Module):
    def __init__(self,
                 input_dim,
                 nb_edge_types,
                 use_3d,
                 gnn_nb_layers,
                 gnn_nb_hidden,
                 rnd_nb_layers,
                 rnd_nb_hidden,
                 rnd_nb_output,
                 init_method=None):
        super(RandomNetwork, self).__init__()
        self.gnn = MyGNN(input_dim, gnn_nb_hidden, gnn_nb_layers, nb_edge_types, use_3d=use_3d, init_method=init_method)
        if gnn_nb_layers == 0:
            in_dim = input_dim
        else:
            in_dim = gnn_nb_hidden

        layers = []
        for _ in range(rnd_nb_layers):
            curr_layer = nn.Linear(in_dim, rnd_nb_hidden)
            if init_method is not None:
                init_network(curr_layer, init_method)
            layers.append(curr_layer)
            in_dim = rnd_nb_hidden

        self.layers = nn.ModuleList(layers)
        self.final_layer = nn.Linear(in_dim, rnd_nb_output)
        self.act = nn.ReLU()

    def forward(self, g_emb):
        X = self.gnn.get_embedding(g_emb, detach=False)
        for i, l in enumerate(self.layers):
            X = self.act(l(X))
        return self.final_layer(X)

