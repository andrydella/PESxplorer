# Automatic exploration of a general PES from CREST metadynamics or MSREACT workflow output
# 30 Ott 24
# Andrea Della Libera & Sarah N. Elliott

# IN PROGRESS: #
# 1. Use EStokTP grid search for TS finding as fallback to GSM
# 2. Setup AMech working dirs
# 3. Identify also minima with findpeaks(-x) to separate multistep rxns in el. steps
# 4. for Single ended, need to check the product connectivity

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

#############################

# Functions
#0. Setup CREST calculation
def crest_calc(filename,crest_out_name,calcs,charge,spin,enewin,path_to_log,log_name):
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
    write_log(f"\n####\nWorking in {crest_dir}\n####\n",log_name)

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
            write_log(f"Pre processing the CREST ensemble file {outfile}",log_name,path_to_log)
            process_ensemble_file(crest_out_name,crest_dir)
            os.system(f"cat {crest_dir}/crestopt.xyz >> {crest_dir}/{crest_out_name}")
            command = 'crest crestopt.xyz --cregen INPUT --ewin ENEWIN --notopo --T 30 &> ensemble.out'
            command = command.replace('INPUT',f'{filename}')
            command = command.replace('ENEWIN',f'{enewin}')
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
        write_log(f"Command is: {command}\nOutput is: {outfile}\n",log_name)
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
    write_log(xyz_folder+xyz_list[0],log_name)
    write_log(f'n atoms: {n_atoms}',log_name)
    with open('natoms_unite.txt','w') as f: f.write(str(n_atoms))

    for xyzfil in xyz_list:
        write_log('Working on: '+xyzfil,log_name)
        with open(xyz_folder+xyzfil,'r') as f:
            lines = [line.strip() for line in f.readlines()]
        for i,line in enumerate(lines):
            try:
                n_at = int(line)
                if n_at == n_atoms:
                    isomers.append([li for li in lines[i:i+n_atoms+2]])
                    n_isom += 1
                    write_log('New isomer! n isom '+str(n_isom),log_name)
                else:
                    pass
            except:
                pass
    write_xyzlist(isomers,'allcrestprods_unite.xyz')


# 1. Unite metadynamics results in one single text file and create list of isomers
#    Names of md folders should all have the same root name
def unite_xyz(md_path,md_name,log_name): # COmputes high enes for bimoleculars as they are considered together
    # isomers = []
    # n_isom = 0
    # filename = 'crest_products.xyz' # Output file CREST MD

    # if float(crest_vers) >=3.0:
    os.system(f"cp crest_ensemble.xyz allcrestprods_unite.xyz")
    with open("allcrestprods_unite.xyz",'r') as f:
        lines = [line.strip() for line in f.readlines()]
    n_atoms = int(lines[0])
    write_log(f"Crest output was found in {md_path}{md_name}\n",log_name)
    write_log(f'n atoms: {n_atoms}',log_name)
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
def find_reactions(log_name):
    geos,well_dct,reac_dct = finder(os.getcwd())
    write_pickle(geos,"geos")
    write_pickle(well_dct,"wells")
    write_pickle(reac_dct,"rxns")
    write_log(f"Number of species considered: {len(geos)}",log_name)
    write_log(f"Number of reactions found: {len(reac_dct.keys())}",log_name)

# 4. Convert data structures and setup GSM calculations
def setup_gsm(filinp,model_gsm,spin,charge,log_name):
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

        write_log(f"Working on reaction {reac_prod}",log_name)
        fold_name = f'GSM_FOLDS/{reac_prod[0]}_gsm_fold'
        write_log(f"Making folder {fold_name}",log_name)
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
        write_log(f"Current counter for reactions of {reac_prod[0]}: {gsm_num}",log_name)
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
def run_gsm(is_ssm,gsm_theory,reacs_set,prods_set,path_to_log,log_name):
    # os.chdir('GSM_FOLDS')
    # gsm_folds = [el for el in os.listdir() if 'gsm_fold' in el]
    # gsm_folds.sort()
    what_am_i = 'GSM'
    if is_ssm:
        what_am_i = 'SSM'
    rxns,geos,geo_dct,wells,gsm_paths,reacs_set,prods_set = setup_from_pickles(reacs_set, prods_set)
    # for fold in gsm_folds:
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
        maindir = os.getcwd()
        for path in gsm_paths:
            if path[0] == rname and path[1] == pname:
                write_log(os.path.join(path[2], f'initial{path[3]}.xyz'),log_name,path_to_log)
                if f'initial{path[3]}.xyz' not in os.listdir(path[2]):
                    write_log("initial file does not exist, moving on",log_name,path_to_log)
                    continue
        # if os.listdir(f'{fold}/scratch/') is not []:
                os.chdir(path[2])
                write_log(f"In folder {path[2]}, is ssm on? {is_ssm}",log_name,path_to_log)

                os.system("cp inpfileq.gsm inpfileq")
                if is_ssm:
                    os.system("cp inpfileq.ssm inpfileq")
            # inputss = [el for el in os.listdir('scratch') if el.startswith('initial')]
            # for i in range(len(inputss)):
                gsm_num = path[3]
                gsm_int = int(path[3])
                inputss = f'initial{gsm_num}.xyz'
                # gsm_num = str(i+1).zfill(4)
                write_log(f'Running gsm.{gsm_theory} {gsm_num} in {path[2]}',log_name,path_to_log)
                if f"stringfile.xyz{gsm_num}" not in os.listdir():
                    write_log(f"ts file NOT found, running {what_am_i} now",log_name,path_to_log)
                    if gsm_theory == 'xtb':
                        command = f"./gsm.orca {gsm_int} 10 &> out{gsm_num}.log"
                    else:
                        command = f"gsm {gsm_int} 10 &> out{gsm_num}.log"
                    with subprocess.Popen(command, 
                                          stdout=subprocess.PIPE, shell=True) as p:
                        p.communicate()  
                        p.wait()
                else:
                    write_log(f"Stringfile FOUND, skipping {what_am_i}",log_name,path_to_log)

                os.chdir(maindir)

        # else: 
        #     write_log(f'{fold} empty folder',path_to_log)

##################################
def are_gras_same(gras, gras2):
    ret = False
    if len(gras) == 2 and len(gras2) == 2:
        if automol.graph.isomorphic(gras[0], gras2[0]):
            if automol.graph.isomorphic(gras[1], gras[1]):
                ret = True
        elif automol.graph.isomorphic(gras[0], gras2[1]):
            if automol.graph.isomorphic(gras[1], gras2[0]):
                ret = True 
    elif len(gras) == 1 and len(gras2) == 1:
        if automol.graph.isomorphic(gras[0], gras2[0]):
                ret = True
    return ret


def postproc(gsm_theory,reacs_set,prods_set,path_to_log,log_name):
    # My postproc data structure is a dictionary
    postproc_dct, not_concerted, not_found = {},{},{}
    failed_gsm = []

    (
        rxns, geos, geo_dct,
        wells, gsm_paths, reacs_set, prods_set
    ) = setup_from_pickles(reacs_set, prods_set)
    check_rxns = {key.replace(' ',''): value for key, value in rxns.items()}

    edge_lst = []
    edge_ene_lst = []
    for rxn_string, rxn_info in check_rxns.items():
        rname, pname = rxn_string.split('=')
        if rname not in reacs_set:
            continue
        if pname not in prods_set:
            continue

        # Go through gsm_folds and find single step reactions
        for path in gsm_paths:
            if path[0] == rname and path[1] == pname:
                # Read stringfile to find barrier geos and enes
                if f'stringfile.xyz{path[3]}' not in os.listdir(path[2]):
                    write_log(f"{path[2]}/stringfile.xyz{path[3]} does not exist",log_name)
                    failed_gsm.append(os.path.join(path[2], f'stringfile.xyz{path[3]}'))
                    continue
                with open(os.path.join(path[2], f'stringfile.xyz{path[3]}'), 'r') as f:
                    traj_str = f.read()
                traj = automol.geom.from_xyz_trajectory_string(traj_str)
                enes = [float(energy) for _, energy in traj]
                ene_max_idx = enes.index(max(enes))
                is_bless = "TS"
                if ene_max_idx in [0,len(enes)-1]:
                    is_bless = "BLESS"

                # find well geos and enes (add first and last structure if not found)
                ene_well_idxs, _ = find_peaks([-ene for ene in enes])
                ene_well_idxs = list(ene_well_idxs)
                if len(ene_well_idxs) < 1:
                    ene_well_idxs.insert(0,0)
                elif ene_well_idxs[0] != 0:
                    ene_well_idxs.insert(0,0)
                if ene_well_idxs[-1] != len(enes) - 1:
                    ene_well_idxs.append(len(enes) - 1)

                # Look at the path from each well
                for well_num in range(len(ene_well_idxs)-1):
                    # split the trajectory to only include the well's elementary channel
                    elem_traj = traj[ene_well_idxs[well_num]:ene_well_idxs[well_num+1]+1]
                    elem_enes = [float(energy) for _, energy in elem_traj]
                    elem_ene_max_idx = elem_enes.index(max(elem_enes))
                    # read the elementary barrier and well geometries and energies
                    elem_ts_geo, elem_ene_max = elem_traj[elem_ene_max_idx]
                    rgeo = elem_traj[0][0] 
                    pgeo = elem_traj[-1][0]
                    rgras = automol.graph.connected_components(automol.geom.graph(rgeo))
                    pgras = automol.graph.connected_components(automol.geom.graph(pgeo))
                    if are_gras_same(rgras, pgras):
                        # this makes sure that the barrier between two wells is not 
                        # simply a reorientation or conformational barrier between
                        # two conformations of the same species
                        continue

                    # identify which species numbers from wells.pickle the wells correspond with
                    rspc = None
                    pspc = None
                    for spc, gras in wells.items():
                        if rspc is not None and pspc is not None:
                            break
                        if are_gras_same(gras, rgras):
                            rspc = spc
                        elif are_gras_same(gras, pgras):
                            pspc = spc
                    # if this elementary reaction was already found through
                    # another trajectory we don't want to double add it
                    if f'{rspc} = {pspc}' in rxns.keys():
                        continue

                    # add well if it was not a part of well.pickle
                    if rspc is None:
                        new_idx = 1 + max([int(key.split('_')[1]) for key in wells.keys()])
                        rspc = f'species_{new_idx}' 
                        wells[rspc] = rgras
                        geo_dct[rspc] = rgeo
                    if pspc is None:
                        new_idx = 1 + max([int(key.split('_')[1]) for key in wells.keys()])
                        pspc = f'species_{new_idx}' 
                        wells[pspc] = pgras
                        geo_dct[pspc] = pgeo

                    # add reaction to rxn.pickle if it was not originally elementary
                    if rspc != rname or pspc != pname:
                        iso_dct, frm_bnd_lst, brk_bnd_lst = automol.reac.arbitrary_reactions(
                            rgras, pgras)
                        if iso_dct:
                            rxns[rspc + ' = ' + pspc] = (iso_dct, frm_bnd_lst, brk_bnd_lst)
                        else:
                            print(f'elementary step for {rname}={pname} could not be identified')
                            not_found[f'{rname}={pname}'] = len(ene_well_idxs)
                            print(automol.geom.string(rgeo))                      
                            print(automol.geom.string(pgeo))

                    # human understandable reaction print
                    prnt_str = ' + '.join(
                        automol.chi.smiles(automol.graph.chi(gra)) 
                            for gra in rgras)
                    prnt_str += ' = '
                    prnt_str += ' + '.join(
                        automol.chi.smiles(automol.graph.chi(gra)) 
                            for gra in pgras)
                        
                    # add elementary reaction to postprocess dictionary                
                    edge_lst.append((rspc, pspc))
                    edge_ene_lst.append(float(elem_ene_max))
                    postproc_dct[f'{rname}={pname}'] = (
                        prnt_str, is_bless,
                        elem_ene_max, elem_ts_geo, elem_traj)
                    
                    write_pickle(postproc_dct,"postproc")
                    # write_pickle(not_concerted,"not_concerted")
                    write_pickle(not_found,"not_found")
                    write_pickle(edge_lst,"edge_lst")
                    write_pickle(edge_ene_lst,"edge_ene_lst")
                    write_pickle(failed_gsm,"failed_paths")
                    write_pickle(rxns, "pp_rxns")
                    write_pickle(wells, "pp_wells")
                    write_pickle(geo_dct, "pp_geos")

    with open(f"{gsm_theory}-gsm-output.csv","w") as f:
        for key,value in postproc_dct.items():
            stuff = ('\t').join(list(map(str,value[:3])))
            f.write(f"{key}\t{stuff}\n")


def prepare_amech_therm(input_str, well_chg, log_name, fs_path, path='.'):

    # get postprocessed info
    (
        _, well_geos, well_gra_dct, well_ene_dct, 
        name_dct, postproc_dct
    ) = setup_from_pp_pickles(path)
    
    # fill in info that automated_insert needs
    insert_dct = INSERT_DCT
    insert_dct['save_filesystem'] = f'{fs_path}/SAVE'
    insert_dct['input_string'] = input_str
    insert_dct['charge'] = well_chg
    insert_dct['program'] = 'gaussian09'
    insert_dct['method'] = 'b3lyp'
    insert_dct['basis'] = 'sto-3g'
    insert_dct['orb_res'] = 'RU'

    # update species filesystem
    save_wells = list(well_ene_dct.keys())
    _, spc_csv_str = prepare_amech_spcdb(
        save_wells, insert_dct, well_chg, well_geos, well_gra_dct,
        well_ene_dct, name_dct, postproc_dct, log_name)
    mech_str = "REACTIONS     CAL/MOLE     MOLES\nEND\n\n\n"

    # set up amech input files
    run_dat_str = read_template('automech_templates', 'run_thermo.dat')
    theo_dat_str = read_template('automech_templates', 'theory.dat')
    models_dat_str = read_template('automech_templates', 'models.dat')
    spc_dat_str = read_template('automech_templates', 'species.dat')
    run_dat_str = run_dat_str.replace('RUN_PREFIX', f'{fs_path}/RUN')
    run_dat_str = run_dat_str.replace('SAVE_PREFIX', f'{fs_path}/SAVE')
    run_dat_str = run_dat_str.replace('NUM_SPECIES', str(len(save_wells)))
    write_log(spc_csv_str, 'species.csv', path_to_log=f"{path}/inp/", create=True)
    write_log(mech_str, 'mechanism.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(run_dat_str, 'run.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(models_dat_str, 'models.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(theo_dat_str, 'theory.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(spc_dat_str, 'species.dat', path_to_log=f"{path}/inp/", create=True)


def prepare_amech_kin(input_str, well_chg, well_spin, log_name, fs_path, path='.'):

    # get postprocessed info
    (
        rxns, well_geos, well_gra_dct, well_ene_dct, 
        name_dct, postproc_dct
    ) = setup_from_pp_pickles(path)

    # fill in info that automated_insert needs
    insert_dct = INSERT_DCT
    insert_dct['save_filesystem'] = f'{fs_path}/SAVE'
    insert_dct['input_string'] = input_str
    insert_dct['charge'] = well_chg
    insert_dct['program'] = 'gaussian09'
    insert_dct['method'] = 'b3lyp'
    insert_dct['basis'] = 'sto-3g'
    insert_dct['orb_res'] = 'RU'

    # update species filesystem
    reached_wells = [
        well.replace(' ','') for wells in list(rxns.keys()) for well in wells.split('=')]
    reached_wells = set(reached_wells)
    well_info_dct, spc_csv_str = prepare_amech_spcdb(
        insert_dct, reached_wells, well_chg, well_geos, well_gra_dct,
        well_ene_dct, name_dct, log_name)
    
    # update reactions filesystem
    insert_dct['ts_mult'] = well_spin + 1
    insert_dct['saddle'] = True
    mech_dat_str = prepare_amech_rxnsdb(
        insert_dct, well_info_dct, rxns, well_gra_dct,
        well_ene_dct, name_dct, postproc_dct, log_name)
    
    # set up amech input files
    run_dat_str = read_template('automech_templates', 'run_rate.dat')
    theo_dat_str = read_template('automech_templates', 'theory.dat')
    models_dat_str = read_template('automech_templates', 'models.dat')
    spc_dat_str = read_template('automech_templates', 'species.dat')
    run_dat_str = run_dat_str.replace('RUN_PREFIX', f'{fs_path}/RUN')
    run_dat_str = run_dat_str.replace('SAVE_PREFIX', f'{fs_path}/SAVE')
    run_dat_str = run_dat_str.replace('NUM_CHANNELS', str(len(mech_dat_str.splitlines())-2))
    write_log(spc_csv_str, 'species.csv', path_to_log=f"{path}/inp/", create=True)
    write_log(mech_dat_str, 'mechanism.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(run_dat_str, 'run.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(models_dat_str, 'models.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(theo_dat_str, 'theory.dat', path_to_log=f"{path}/inp/", create=True)
    write_log(spc_dat_str, 'species.dat', path_to_log=f"{path}/inp/", create=True)


def prepare_amech_spcdb(
        insert_dct, save_wells, well_chg, well_geos, well_gra_dct,
        well_ene_dct, name_dct, log_name):

    # hardwired script import until i finish setting up script as callable cli thing
    current_dir = os.getcwd()
    bin_path = os.path.join(current_dir, '/home/elliott/Packages/AutoMech/mechdriver/bin')
    sys.path.insert(0, bin_path)
    import automated_insert

    well_info_dct = {}
    spc_csv = 'name,inchi,mult,charge'
    saved_chis = []

    #  loop over species and add them to the amech filesystem
    for well, ene in well_ene_dct.items():
        if well not in save_wells:
            write_log(f'not saving {well} because no rxn uses it', log_name)
            continue
        if well in well_info_dct:
            well_info = well_info_dct[well]
        else:
            well_info = ()
        for i, gra in enumerate(well_gra_dct[well]):
            if len(well_info) < i + 1:
                mults_allowed = automol.graph.possible_spin_multiplicities(gra)
                well_mult = mults_allowed[0]
                chi = automol.graph.chi(gra)
                well_info += ((chi, well_chg, well_mult),)
            well_mult = well_info[i][2]
            chi = well_info[i][0]
            if chi not in saved_chis:
                insert_dct['mult'] = well_mult
                insert_dct['inchi'] = chi
                insert_dct['output_string'] = automol.geom.xyz_string(
                    automol.geom.subgeom(well_geos[well], automol.graph.atom_keys(gra)),
                    f'{ene/len(well_gra_dct[well]):.8f}')
                automated_insert.main(insert_dct)
                name = name_dct[well].split('+')[i].replace(' ','')
                spc_csv += f"\n{name},'{automol.graph.chi(gra)}',{str(well_mult)},{str(well_chg)} ! {well}"
                saved_chis.append(chi)
            well_info_dct[well] = well_info
    return well_info_dct, spc_csv


def prepare_amech_rxnsdb(
        insert_dct, well_info_dct, rxns, well_gra_dct,
        well_ene_dct, name_dct, postproc_dct, log_name):
    
    # hardwired script import until i finish setting up script as callable cli thing
    current_dir = os.getcwd()
    bin_path = os.path.join(current_dir, '/home/elliott/Packages/AutoMech/mechdriver/bin')
    sys.path.insert(0, bin_path)
    import automated_insert

    mech_str = "REACTIONS     CAL/MOLE     MOLES"
    for name, rxn in rxns.items():

        name = name.replace(' ','')
        if name not in postproc_dct:
            continue

        rxn_dets = postproc_dct[name]
        _, cla, barrier, ts_geo, _ = rxn_dets

        reacs, prods = name.split('=')
        if name_dct[reacs] + ' = ' + name_dct[prods] + ' ' in mech_str:
            continue

        rct_chis = [info[0] for info in well_info_dct[reacs]]
        rct_mults = [info[2] for info in well_info_dct[reacs]]
        rct_chgs = [info[1] for info in well_info_dct[reacs]]
        prd_chis = [info[0] for info in well_info_dct[prods]]
        prd_mults = [info[2] for info in well_info_dct[prods]]
        prd_chgs = [info[1] for info in well_info_dct[prods]]
        rxn_chis = (rct_chis, prd_chis,)
        rxn_mults = (rct_mults, prd_mults,)
        rxn_chgs = (rct_chgs, prd_chgs,)
        rct_gra = well_gra_dct[reacs]
        rct_gra = automol.graph.union_from_sequence(rct_gra)
        prd_gra = well_gra_dct[prods]
        prd_gra = automol.graph.union_from_sequence(prd_gra)
        if cla == 'TS':
            ene = well_ene_dct[reacs] * float(barrier)*627.51
            insert_dct['mult'] = rxn_mults
            insert_dct['charge'] = rxn_chgs
            insert_dct['inchi'] = rxn_chis
            insert_dct['output_string'] = automol.geom.xyz_string(
                ts_geo, f'{ene:.8f}')
            insert_dct['forming_bonds'] = rxn[1]
            insert_dct['breaking_bonds'] = rxn[2]
            insert_dct['rct_gra'] = rct_gra
            insert_dct['prd_gra'] = prd_gra
            if len(rct_chis) > 1:
                insert_dct['rxn_class'] = "unclassified_bimol"
            else:
                insert_dct['rxn_class'] = "unclassified_unimol"
            try:
                automated_insert.main(insert_dct)
                mech_str += '\n' + name_dct[reacs] + ' = ' + name_dct[prods] + '    1.0  0  0'
            except:
                write_log('failed to save', log_name)
                write_log(name_dct[reacs] + ' = ' + name_dct[prods], log_name)
    
    mech_str += "\nEND\n"
    return mech_str

##################################


def main():
    message='''
Need to provide an argument to the explorer!
Possible commands are:
- runcrest -> runs preliminary CREST calculations
- [Deprecated] unite -> reads the output of CREST calculations and collects results in a trajectory text file
- [Deprecated] bond_check -> select unique isomers and eventually bimolecular products
- selpaths -> create global_gsm.txt
- [Deprecated] filter -> creates filter_gsm.txt
- setup -> sets up GSM folders for each isomer
- rungsm -> runs all GSM calcs
- runssm -> runs SSM calculations as fallback
- postproc -> runs postprocess on GSM_FOLDS
'''
 #   commands = ["runcrest","unite","bond_check","selpaths","filter","setup","rungsm"]
    commands = ["runcrest","selpaths","setup","rungsm","runssm","postproc","prepamech"]
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
        if "-1" not in config["reacs"]:
            reacs_set = species_parser("reacs",config)
    if "prods" in config:
        if "-1" not in config["prods"]:
            prods_set = species_parser("prods",config)
    enewin = config['confs_window']


    log_strings = [".log",command]
    dirs_lst = [fil for fil in os.listdir() if all(x in fil for x in log_strings)]
    fil_nums = []
    if not dirs_lst: 
        log_name = command + ".log_1"
    else:
        for direc in dirs_lst:
            fil_nums.append(int(direc.split("_")[1])) 
        log_name = f"{command}.log_{max(fil_nums) + 1}"
    open(log_name,'w').close() # Create logfile
    path_to_log = os.getcwd()+'/'
    write_log("Read the input file input.dat",log_name)
    # Write info on current run in decent way
    write_log(lines,log_name)

    at_num = {'1':'H','6':'C','8':'O','7':'N'}


# #### 1 ###
    if command == 'runcrest':
        crest_calc('input-stru.xyz', crest_out_name, ['opt','msreact','ensemble'],charge,spin,enewin,path_to_log,log_name)
    #elif command == 'unite':
        unite_xyz(crest_md_path,crest_out_name,log_name) # -> allcrestprods_unite.xyz
        os.system('cp allcrestprods_unite.xyz allcrestprods_sort.xyz')
#        sort_xyz('allcrestprods_unite.xyz','crest_ensemble.xyz',charge,spin) # -> allcrestprods_sort.xyz
# Sort does not handle correclty bimol products
# #### 2 ###
    elif command == 'selpaths':
        find_reactions()
# #### 3 ###
    elif command == 'setup':
        setup_gsm('allcrestprods_sort.xyz',model_gsm,spin,charge,log_name)
# #### 4 ###
    elif command == 'rungsm':
        run_gsm(False,gsm_theory,reacs_set,prods_set,path_to_log,log_name)
    elif command == 'runssm':
        run_gsm(True,gsm_theory,reacs_set,prods_set,path_to_log,log_name)
    elif command == 'postproc':
        postproc(gsm_theory,reacs_set,prods_set,path_to_log,log_name)
    elif command == 'prepamech':
        prepare_amech_kin(
            '\n'.join(lines), int(charge), int(spin), 
            log_name, '/home/elliott/projects/explorer/5_andrearun/GSM_FS/',)


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
