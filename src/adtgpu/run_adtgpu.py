import os
import re
import sys
import argparse
import subprocess

if __name__ == '__main__':
    #Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--receptor',dest='receptor_file',
        help='Path to receptor file')
    parser.add_argument('-s','--smiles', dest='smiles',
        help='List of smiles strings.\nIn single quotes and seperated by commas')
    parser.add_argument('-d','--rundir', dest='run_dir',
        help='Run directory')
    args = parser.parse_args()

    #Debug mode
    DEBUG=False

    #Check that input file path exist
    if not os.path.exists(args.receptor_file):
        exit("Receptor file does not exist: {}".format(args.receptor_file))   
    #Create output dirs
    ligands_dir="/ligands"
    if not os.path.exists(args.run_dir): 
        os.makedirs(args.run_dir)
    if not os.path.exists(args.run_dir+ligands_dir):
        os.makedirs(args.run_dir+ligands_dir)
    #TODO ^ to leave here or keep in gcpn code???? 

    if(DEBUG): print("Inside adt script")

    #Executable paths
    obabel="/gpfs/alpine/syb105/proj-shared/Personal/manesh/BIN/openbabel/summit/build/bin/obabel"
    adt="/gpfs/alpine/syb105/proj-shared/Personal/gabrielgaz/Apps/summit/autoDockGPU2/bin/autodock_gpu_64wi"

    #Parse smiles input into array
    smiles=re.sub('\ |\'', '', args.smiles[1:-1]).split(",")
    if(DEBUG): print("List of smiles:\n{}".format('\n'.join(smiles)))

    #Loop over smile strings to convert to pdbqt
    ligs_list=[]
    sm_counter=1
    for smile in smiles:
        #Create temp directory needed for obabel
        tmp_file=args.run_dir+ligands_dir+"/ligand"+str(sm_counter)+".smile"
        with open(tmp_file,'w') as f:
            f.write(smile+'\n')
        
        #Create name for output pdbqt file
        ligand_out=args.run_dir+ligands_dir+"/ligand"+str(sm_counter)+".pdbqt"

        #Create run command and execute
        cmd=obabel+" --gen3d --partialcharge gasteiger --addfilename -ismi "
        cmd+=tmp_file+" -opdbqt > "+ligand_out
        if(DEBUG): print("\nCmd to run:\n{}".format(cmd))
        subprocess.Popen(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).wait()

        #Clean up and increment smile counter
        os.remove(tmp_file)
        ligand_store_file=ligand_out.split('/')[-1][:-6]
        ligs_list.append(ligand_store_file)
        sm_counter+=1
    
    #Get stub name of receptor and field file
    receptor_dir='/'.join(args.receptor_file.split('/')[:-1])
    receptor_stub=args.receptor_file.split('/')[-1][:-6] #rm .pdbqt=6
    if(DEBUG): print("\nReceptor dir:  {}".format(receptor_dir))
    if(DEBUG): print("Receptor stub: {}".format(receptor_stub))
    receptor_field=receptor_stub+".maps.fld"

    #Create run file for Autodock-gpu
    run_file=args.run_dir+"/ligs_list.runfile"
    run_file_lbl="ligs_list.runfile"
    with open(run_file,'w') as f:
        f.write(receptor_field+'\n')
        for lig in ligs_list:
            f.write("ligands/"+lig+".pdbqt\n")
            f.write("ligands/"+lig+'\n')

    #Copy map files to run dir
    cmd="cp "+receptor_dir+"/"+receptor_stub+"* "+args.run_dir
    if(DEBUG): print("\nCopy cmd to run: {}".format(cmd))
    subprocess.Popen(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).wait()

    #Set up autodock-gpu run command
    cmd="/gpfs/alpine/syb105/proj-shared/Personal/gabrielgaz/Apps/summit/autoDockGPU2/bin/autodock_gpu_64wi -filelist "+run_file_lbl+" -nrun 10"
    if(DEBUG): print("\nAutodock cmd to run: {}".format(cmd))

    #Run autodock-gpu (in run_dir and move back)
    cur_dir=os.getcwd()
    os.chdir(args.run_dir)
    subprocess.Popen(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).wait()
    os.chdir(cur_dir)

    #Create output file for final scores
    output_file=args.run_dir+"/adt_affinity_vals"
    with open(output_file,'w') as f:
        for lig in ligs_list:
            #Parse for final score
            lig_path=args.run_dir+ligands_dir+"/"+lig+".dlg"
            if not os.path.exists(lig_path):
                print("ERROR: No such file {}".format(lig_path))
            else: 
                grep_cmd = "grep -2 \"^Rank \" "+lig_path+" | head -5 | tail -1 | cut -d \'|\' -f2 | sed \'s/ //g\'"
                grep_out=os.popen(grep_cmd).read()
                f.write(grep_out)

