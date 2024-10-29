# Automatic exploration of a general PES from CREST metadynamics or MSREACT workflow output
# 29 Ott 24
# Andrea Della Libera & Sarah N. Elliott

# IN PROGRESS: #
# 1. Use EStokTP grid search for TS finding as fallback to GSM
# 2. Setup AMech working dirs

# TO DO: #
# 1. Interface MOLGEN output

# Modules import
import os
import sys
import subprocess
import re
#import csv
from copy import deepcopy
from scipy.signal import find_peaks
from auxiliary_funcs import *
from old_bondcheck import *
from find_from_traj import finder
from mechanalyzer.builder._names import functional_group_name

#############################

# Functions
#0. Setup CREST calculation
def crest_calc(filename,crest_out_name,calcs,charge,spin,path_to_log):
    # Setup crest subfolder
    crest_dir_prefix = "crest_calc"
    dirs_lst = [dir for dir in os.listdir() if crest_dir_prefix in dir]
    folder_nums = []
    if not dirs_lst: crest_dir = crest_dir_prefix+"_1"
    else:
        for direc in dirs_lst:
            folder_nums.append(int(direc.split("_")[2])) 
        crest_dir = f"{crest_dir_prefix}_{max(folder_nums) + 1}"
    os.system(f"mkdir -p {crest_dir}") 
    write_log(f"\n####\nWorking in {crest_dir}\n####\n")

    os.system(f'cp {filename} {crest_dir}')
    os.system(f'echo {charge} > {crest_dir}/.CHRG')
    os.system(f'echo {spin} > {crest_dir}/.UHF')

    for calc_type in calcs:
        if calc_type == 'opt':
            command = 'crest INPUT --opt --gnf2'
            command = command.replace('INPUT',f'{filename}')
            outfile = 'crestopt.xyz'
        elif calc_type == 'msreact':
            command = 'crest INPUT --msreact --msnshifts 300 --msnshifts2 15 --msmolbar --T 30 –ewin 500. &> msreact.out'
            command = command.replace('INPUT',f'{filename}')
            outfile = crest_out_name
        elif calc_type == 'ensemble':
            write_log(f"Pre processing the CREST ensemble file {outfile}",path_to_log)
            process_ensemble_file(crest_out_name)
            os.system(f"cat crestopt.xyz >> {crest_out_name}")
            command = 'crest crestopt.xyz --cregen INPUT --ewin 50. --notopo --T 30 &> ensemble.out'
            command = command.replace('INPUT',f'{filename}')
            outfile = 'crest_ensemble.xyz'
        elif calc_type == 'sp':
            command = 'crest --for INPUT --prop singlepoint --ewin 50. --notopo &> energy.out'
            command = command.replace('INPUT',f'{filename}')
            outfile = 'crest_ensemble.xyz'
        else:
            print(f'CREST calculation type {calc_type} is not currently defined')

        crest_run = f'''cd {crest_dir}
                    {command}
                    '''
        filename = outfile
        write_log(f"Command is: {command}\nOutput is: {outfile}\n")
        with subprocess.Popen(crest_run, stdout=subprocess.PIPE, shell=True) as p:
            p.communicate()
            p.wait()

    os.system(f"cp {crest_dir}/{outfile} .") 
    # if 'sp' not in calcs:
    #     os.system(f"cat {crest_dir}/crestopt.xyz >> {outfile}") 

    return

# 1.b Get structures from MOLGEN output and also consider bimolecular products
# TODO, test needed, I believe there is still something wrong
def unite_molgen_xyz(xyz_folder,suffix="opt.xyz"):   
    isomers = []
    n_isom = 0
    
    xyz_list = [fil for fil in os.listdir(xyz_folder) if suffix in fil]

    with open(xyz_folder+xyz_list[0],'r') as f:
        lines = [line.strip() for line in f.readlines()]
    n_atoms = int(lines[0])
    # Write logfile here
    write_log(xyz_folder+xyz_list[0])
    write_log(f'n atoms: {n_atoms}')
    with open('natoms_unite.txt','w') as f: f.write(str(n_atoms))

    for xyzfil in xyz_list:
        write_log('Working on: '+xyzfil)
        with open(xyz_folder+xyzfil,'r') as f:
            lines = [line.strip() for line in f.readlines()]
        for i,line in enumerate(lines):
            try:
                n_at = int(line)
                if n_at == n_atoms:
                    isomers.append([li for li in lines[i:i+n_atoms+2]])
                    n_isom += 1
                    write_log('New isomer! n isom '+str(n_isom))
                else:
                    pass
            except:
                pass
    write_xyzlist(isomers,'allcrestprods_unite.xyz')


# 1. Unite metadynamics results in one single text file and create list of isomers
#    Names of md folders should all have the same root name
def unite_xyz(md_path,md_name): # COmputes high enes for bimoleculars as they are considered together
    # isomers = []
    # n_isom = 0
    # filename = 'crest_products.xyz' # Output file CREST MD

    # if float(crest_vers) >=3.0:
    os.system(f"cp crest_ensemble.xyz allcrestprods_unite.xyz")
    with open("allcrestprods_unite.xyz",'r') as f:
        lines = [line.strip() for line in f.readlines()]
    n_atoms = int(lines[0])
    write_log(f"Crest output was found in {md_path}{md_name}\n")
    write_log(f'n atoms: {n_atoms}')
    with open('natoms_unite.txt','w') as f: f.write(str(n_atoms))

    # else:
    #     if 'allcrestprods_unite.xyz' not in os.listdir():
    #     #if already run, don't recreate allprods (it takes a while...)
    #         md_list = [fold+'/' for fold in os.listdir(md_path) if md_name in fold]

    #         with open(md_path+md_list[0]+filename,'r') as f:
    #             lines = [line.strip() for line in f.readlines()]
    #         n_atoms = int(lines[0])
    #         # Write logfile here
    #         write_log(md_path+md_list[0]+filename)
    #         write_log(f'n atoms: {n_atoms}')
    #         with open('natoms_unite.txt','w') as f: f.write(str(n_atoms))

    #         for metadyn in md_list:
    #             write_log('Working on: '+metadyn)
    #             if filename not in os.listdir(md_path+metadyn): 
    #                 write_log('Skipping, not found crest_products.xyz')
    #                 continue
    #             else:
    #                 with open(md_path+metadyn+filename,'r') as f:
    #                     lines = [line.strip() for line in f.readlines()]
    #                 for i,line in enumerate(lines):
    #                     try:
    #                         n_at = int(line)
    #                         if n_at == n_atoms:
    #                             isomers.append([li for li in lines[i:i+n_atoms+2]])
    #                             n_isom += 1
    #                             write_log('New isomer! n isom '+str(n_isom))
    #                         else:
    #                             pass
    #                     except:
    #                         pass
    #         write_xyzlist(isomers,'allcrestprods_unite.xyz')

# 2. Reorder isomers based on energy
def sort_xyz(filinp,crest_out_name,charge,spin):
    isomers = []
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(filinp,n_atoms)
    # Sort by energy
    try:
        isomers.sort(key=lambda x:float(x[1].split()[0]))
    except:
        crest_calc(filinp,crest_out_name,['sp'],charge,spin)
        os.system(f'cp crest_ensemble.xyz {filinp}_sp')
        isomers = read_xyzlist(f'{filinp}_sp',n_atoms)
        isomers.sort(key=lambda x:float(x[1].split()[0]))
    write_xyzlist(isomers,'allcrestprods_sort.xyz')


# 3. Use find.py to get possible reactions
def find_reactions():
    geos,well_dct,reac_dct = finder(os.getcwd())
    write_pickle(geos,"geos")
    write_pickle(well_dct,"wells")
    write_pickle(reac_dct,"rxns")
    write_log(f"Number of species considered: {len(geos)}")
    write_log(f"Number of reactions found: {len(reac_dct.keys())}")

# 4. Convert data structures and setup GSM calculations
def setup_gsm(filinp,model_gsm,spin,charge):
    import automol

    # geos,well_dct,reac_dct = amech_data_structure(filinp)
    geos = read_pickle("geos")
    well_dct = read_pickle("wells")
    reac_dct = read_pickle("rxns")

    os.system('mkdir -p GSM_FOLDS')
    gsm_counters = {spc:0 for spc,_ in well_dct.items()}

    open("reactions.csv","w").close()
    for reaction,stuff in reac_dct.items():
        reac_prod = [spc.strip() for spc in reaction.split("=")]
        reac_prod_idxs = [spc.split("_")[1] for spc in reac_prod]
        reac_idx = int(reac_prod_idxs[0])
        prod_idx = int(reac_prod_idxs[1])
        atom_order = stuff[0]
        b_formed = stuff[1]
        b_broken = stuff[2]
        geo_r = deepcopy(geos[reac_idx])
        geo_p = deepcopy(geos[prod_idx])

        geo_p = swap_atoms(atom_order,geo_p)

        write_log(f"Working on reaction {reac_prod}")
        fold_name = f'GSM_FOLDS/{reac_prod[0]}_gsm_fold'
        write_log(f"Making folder {fold_name}")
        if not os.path.exists(fold_name): 
            os.mkdir(fold_name)
        if not os.path.exists(f'{fold_name}/scratch'): 
            os.mkdir(f'{fold_name}/scratch')
        gsm_path = os.getcwd()+f'/{fold_name}'
        # Copy model data folder
        os.system(f'cp -r {model_gsm} {gsm_path}')

        # Write initial file for each gsm calculation
        gsm_counters[reac_prod[0]] += 1
        gsm_num = str(gsm_counters[reac_prod[0]]).zfill(4)
        write_log(f"Current counter for reactions of {reac_prod[0]}: {gsm_num}")
        initial_str = automol.geom.xyz_trajectory_string([geo_r,geo_p])
        with open(f'{gsm_path}/scratch/initial{gsm_num}.xyz', 'w') as f:
            f.write(initial_str)
        with open(f'{gsm_path}/initial{gsm_num}.xyz', 'w') as f:
            f.write(initial_str)
        with open(f'{gsm_path}/ISOMERS{gsm_num}', 'w') as f:
            for a,b in iter(b_formed):
                f.write(f"ADD  {a} {b}\n")
            f.write("\n")
            for a,b in iter(b_broken):
                f.write(f"BREAK  {a} {b}\n")
        # Also write charge and spin where useful
        with open(f'{gsm_path}/scratch/.UHF','w') as f:
            f.write(spin)
        with open(f'{gsm_path}/scratch/.CHRG','w') as f:
            f.write(charge)
        os.system(f"sed -i 's/SPIN/{int(spin)+1}/g' {gsm_path}/gstart")
        os.system(f"sed -i 's/CHARGE/{charge}/g' {gsm_path}/gstart")
        
        # Write csv file that collect info on reactions found
        with open("reactions.csv","a") as f:
            f.write(f"{reac_prod[0]}, {reac_prod[1]}, {fold_name}, {gsm_num}\n")

# 5. run_gsm()
def run_gsm(is_ssm,gsm_theory,path_to_log):
    os.chdir('GSM_FOLDS')
    gsm_folds = [el for el in os.listdir() if 'gsm_fold' in el]
    gsm_folds.sort()

    if is_ssm:
        what_am_i = 'SSM'
    else:
        what_am_i = 'GSM'

    # Alternative for cycles for running sets of subfolders
    #for fold in gsm_folds[-2:]:
    #for fold in gsm_folds[-4:-2]:
    for fold in gsm_folds:
        if os.listdir(f'{fold}/scratch/') is not []:
            os.chdir(fold)
            print(f"In folder {fold}, is ssm on? {is_ssm}")
            write_log(f"In folder {fold}, is ssm on? {is_ssm}",path_to_log)


            os.system("cp inpfileq.gsm inpfileq")
            if is_ssm:
                os.system("cp inpfileq.ssm inpfileq")

            inputss = [el for el in os.listdir('scratch') if el.startswith('initial')]
            for i in range(len(inputss)):
                gsm_num = str(i+1).zfill(4)
                write_log(f'Running gsm.{gsm_theory} {gsm_num} in {fold}',path_to_log)
                if f"tsq{gsm_num}.xyz" not in os.listdir(f'scratch/'):
                    write_log(f"ts file not found, running {what_am_i} now",path_to_log)
                    if gsm_theory == 'xtb':
                        command = f"./gsm.orca {i+1} 30 &> out{gsm_num}.log"
                    else:
                        command = f"gsm {i+1} 30 &> out{gsm_num}.log"
                    with subprocess.Popen(command, 
                                          stdout=subprocess.PIPE, shell=True) as p:
                        p.communicate()  
                        p.wait()
                else:
                    write_log(f"ts NOT FOUND, skipping {what_am_i}",path_to_log)

            os.chdir('..')

        else: 
            write_log(f'{fold} empty folder',path_to_log)

##################################


def postproc(gsm_theory,reacs_set,prods_set,path_to_log):
    # My postproc data structure is a dictionary
    # keys = "reactant = product"
    # values = case oneTS - mTSs -> read stringfile directly to get this info
    #           list of TS enes
    #           list of TS geos
    #           energy reaction path
    #           list of geometries for reaction path

    postproc_dct,not_concerted = {},{}
    rxns = read_pickle("rxns")

    geos = read_pickle("geos")
    geo_dct = {}
    for i, geo in enumerate(geos):
        geo_dct[f'species_{i}'] = geo

    wells = read_pickle("wells")
    name_dct = {} # dct of species names
    for well, gras in wells.items():
        name_dct[well] = ' + '.join([functional_group_name(
                         automol.graph.inchi(gra)) for gra in gras])
    name_dct2 = {} # dct of smiles
    for well, gras in wells.items():
        name_dct2[well] = ' + '.join([automol.chi.smiles(
                           automol.graph.chi(automol.graph.implicit(gra))) for gra in gras])

    with open('reactions.csv', 'r') as f:
        gsm_paths = f.read()
    gsm_paths = [gsm_path.replace(' ','').split(','
                ) for gsm_path in gsm_paths.splitlines()]

    #gsm_folds = [el for el in os.listdir('GSM_FOLDS') if 'gsm_fold' in el]
    pattern = re.compile(r'(\d+)')

    # If no indication of reacs or prods, consider all of them
    if -1 in reacs_set:
        reacs_set = set(wells.keys())
    if -1 in prods_set:
        prods_set = set(wells.keys())

    edge_lst = []
    edge_ene_lst = []
    for rxn_string, rxn_info in rxns.items():
        rname, pname = rxn_string.replace(' ','').split('=')
        if rname not in reacs_set:
            continue
        if pname not in prods_set:
            continue

        prnt_str = ' + '.join(
            automol.chi.smiles(automol.graph.chi(gra)) 
            for gra in wells[rname])
        prnt_str += ' = '
        prnt_str += ' + '.join(
            automol.chi.smiles(automol.graph.chi(gra)) 
            for gra in wells[pname])
        # Go through gsm_folds and find single step reactions
        for path in gsm_paths:
            if path[0] == rname and path[1] == pname:
                print(os.path.join(path[2], f'stringfile.xyz{path[3]}'))
                if f'stringfile.xyz{path[3]}' not in os.listdir(path[2]):
                    print("strinfile does not exist, moving on")
                    continue
                with open(os.path.join(path[2], f'stringfile.xyz{path[3]}'), 'r') as f:
                    traj_str = f.read()
                traj = automol.geom.from_xyz_trajectory_string(traj_str)

                enes = [float(energy) for _,energy in traj]
                ene_max_idx = enes.index(max(enes))
                ene_max_idxs, _ = find_peaks(enes, height=0.3)
                is_bless = "TS"
                if ene_max_idx in [0,len(enes)-1]:
                    is_bless = "BLESS"

                if len(ene_max_idxs) < 2:
                    ts_geo, ene_max = traj[ene_max_idx]
                    edge_lst.append((rname, pname))
                    edge_ene_lst.append(float(ene_max))
                    postproc_dct[f'{rname}+{pname}'] = (prnt_str,is_bless,ts_geo,ene_max,traj)
                else:
                    print('not concerted!', prnt_str)
                    not_concerted[f'{rname}+{pname}'] = (prnt_str,is_bless,traj)

    write_pickle(postproc_dct,"postproc")
    write_pickle(not_concerted,"not_concerted")
    write_pickle(edge_lst,"edge_lst")
    write_pickle(edge_ene_lst,"edge_ene_lst")

    # for fold in gsm_folds:
    #     fold_split = fold.split('_')
    #     if int(fold_split[1]) in species_set:
    #         print(fold)
    # exit()
    #     if os.listdir(f'{fold}/scratch/') is not []:
    #         os.chdir(fold)
    #         print(f"In folder {fold}, Now running postproc operations...")
    #         write_log(f"In folder {fold}, Now running postproc operations...",path_to_log)

    #         inputss = [el for el in os.listdir('scratch') if el.startswith('initial')]
    #         for inp in inputss:
    #             match = pattern.search(inp)
    #             if match:
    #                 number = int(match.group())  # Convert to int to remove leading zeros
    #             gsm_num = str(number+1).zfill(4)
    #             print(f'Looking at gsm num {gsm_num} in {fold}')
    #             write_log(f'Looking at gsm num {gsm_num} in {fold}',path_to_log)
    #             if f"tsq{gsm_num}.xyz" in os.listdir(f'scratch/'):
    #                 write_log(f"ts file found",path_to_log)
    #                 print(f"ts file found")
    #             else:
    #                 write_log(f"ts NOT FOUND, skipping",path_to_log)
    #                 print(f"ts NOT FOUND, skipping ")

    #         os.chdir('..')

    #     else: 
    #         write_log(f'{fold} empty folder',path_to_log)

##################################


def main():
    message='''
Need to provide an argument to the explorer!
Possible commands are:
- runcrest -> runs preliminary CREST calculations
- unite -> reads the output of CREST calculations and collects results in a trajectory text file
- [Deprecated] bond_check -> select unique isomers and eventually bimolecular products
- selpaths -> create global_gsm.txt
- [Deprecated] filter -> creates filter_gsm.txt
- setup -> sets up GSM folders for each isomer
- rungsm -> runs all GSM calcs
- runssm -> runs SSM calculations as fallback
- postproc -> runs postprocess on GSM_FOLDS
'''
 #   commands = ["runcrest","unite","bond_check","selpaths","filter","setup","rungsm"]
    commands = ["runcrest","unite","selpaths","setup","rungsm","runssm","postproc"]
    if len(sys.argv) != 2:
        print(message)
        exit()
    elif sys.argv[1] in ["help","-help","--help","-h","--h"]:
        print(message)
        exit()
    assert sys.argv[1] in commands
    command = sys.argv[1]

    dirs = os.listdir()
    assert "input.dat" in dirs, "Missing input data file"
    # Read input information
    with open("input.dat","r") as f:
        lines = [line.strip() for line in f.readlines() if (
                 line.strip() != '' and not line.strip().startswith('#'))]
        
    config = {l.split('=')[0].strip():l.split('=')[1].strip() for l in lines}
    
    model_gsm = config['model_gsm']
    # crest_version = config['crest_version']
    crest_md_path = config['crest_md_path']
    crest_out_name = config['crest_out_name']
    charge = config['charge']
    spin = config['spin']
    gsm_theory = config['gsm_theory']
    assert gsm_theory in ['g16','xtb'], "xtb or g16 directives required for GSM"
    reacs_set, prods_set = set([-1]), set([-1])
    if "reacs" in config:
        if config["reacs"] is not "-1":
            reacs_set = species_parser("reacs",config)
    if "prods" in config:
        if config["prods"] is not "-1":
            prods_set = species_parser("prods",config)

    open('logfile.out','w').close() # Create logfile
    path_to_log = os.getcwd()+'/'
    write_log(lines)

    at_num = {'1':'H','6':'C','8':'O','7':'N'}


# #### 1 ###
    if command == 'runcrest':
        crest_calc('input-stru.xyz', crest_out_name, ['opt','msreact','ensemble'],charge,spin,path_to_log)
    elif command == 'unite':
        unite_xyz(crest_md_path,crest_out_name) # -> allcrestprods_unite.xyz
#        sort_xyz('allcrestprods_unite.xyz','crest_ensemble.xyz',charge,spin) # -> allcrestprods_sort.xyz
# Sort does not handle correclty bimol products
# #### 2 ###
    elif command == 'selpaths':
        find_reactions()
# #### 3 ###
    elif command == 'setup':
        setup_gsm('allcrestprods_sort.xyz',model_gsm,spin,charge)
# #### 4 ###
    elif command == 'rungsm':
        run_gsm(False,gsm_theory,path_to_log)
    elif command == 'runssm':
        run_gsm(True,gsm_theory,path_to_log)
    elif command == 'postproc':
        postproc(gsm_theory,reacs_set,prods_set,path_to_log)


 #   unite_molgen_xyz("C4H10_geo/",suffix="opt.xyz")
# # #### 1 ###
#     if command == 'runcrest':
#         crest_calc('input-stru.xyz', ['opt','msreact','ensemble'],charge,spin)
#     elif command == 'unite':
#         unite_xyz(crest_md_path,crest_out_name) #Skips if allcrestprods_unite is already there; protocol depends on version of crest
#         sort_xyz('allcrestprods_unite.xyz',charge,spin) # -> allcrestprods_sort.xyz
# # #### 2 ###
#     elif command == 'bond_check':
#         _,connect,_ = bond_check('allcrestprods_sort.xyz',at_num) # -> unique_bondcheck.xyz
#         #_,_,_ = remove_frags('unique_bondcheck.xyz') # -> nofrags_bondcheck.xyz
#         os.system('cp unique_bondcheck.xyz nofrags_bondcheck.xyz')
#         bond_check3('nofrags_bondcheck.xyz') # -> perm_bondcheck.xyz
#     elif command == 'selpaths':
#         sel_paths('perm_bondcheck.xyz')
# # #### 3 ###
#     elif command == 'filter':
#         filter_gsm('gsm_global.txt')
#     elif command == 'setup':
#         setup_gsm('gsm_filter.txt','perm_bondcheck.xyz',model_gsm)
# # #### 4 ###
#     elif command == 'rungsm':
#         run_gsm()

if __name__=="__main__":
    main()
