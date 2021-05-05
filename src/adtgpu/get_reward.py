import os
import shutil
import subprocess

from rdkit import Chem

def get_reward(states):
    if not isinstance(states, list):
        states = [states]

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
    cmd="python ./src/adtgpu/run_adtgpu.py -r ./src/adtgpu/receptor/NSP15_6W01_A_1_F_receptor.pdbqt -s \"" + str(smiles) + "\" -d "+adttmp    
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

