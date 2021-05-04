from rdkit.Chem.Descriptors import MolLogP

def get_logp_scores(states):
    if not isinstance(states, list):
        states = [states]
    scores = [MolLogP(state) for state in states]
    return scores
    