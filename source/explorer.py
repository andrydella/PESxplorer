# Automatic exploration of a general PES from CREST metadynamics output
# 26/05/24
# Andrea Della Libera

# IN PROGRESS: #
# 1. Interface with molgen output

# TO DO: #
# 1. Crivello di Eratostene in bond_check x scalabilità
# 2. Eliminare ridondanze calcolo permutazioni
# 3. Eliminare ridondanze calcolo matrice distanze
# 4. Formatta meglio logfile
# 5. Generalize once more CREST output reading as even with new version more iterations are required! 

# Modules to be imported
import os
from itertools import permutations
import csv
import math as m
import generalized_permutations as gen_perm

# Auxilliary functions (| in the comments means 'or')
# a. write xyz from list in the form [natom,energy|empty,coords1,...,coordsn]
def write_xyz(xyz_list,filname):
    with open(filname,'w') as f:
        for line in xyz_list:
            f.write(line+'\n')

# a.1 append xyz from list in the form [natom,energy|empty,coords1,...,coordsn]
def append_xyz(xyz_list,filname):
    with open(filname,'a') as f:
        for line in xyz_list:
            f.write(line+'\n')

# b. Write log file
def write_log(stuff_to_write):
    with open('logfile.out','a') as f:
        if type(stuff_to_write) == type('str'): f.write(stuff_to_write+'\n')
        elif type(stuff_to_write) == type(1): f.write(str(stuff_to_write)+'\n')
        else: f.writelines(stuff_to_write) # For lists

# c. Write all xyz in isomers into file
def write_xyzlist(xyz_list,filname):
    with open(filname,'w') as f:
        for isomer in xyz_list:
            for line in isomer: f.write(line+'\n')

# d. create isomers list from xyz list file
def read_xyzlist(filname,n_atoms):
    isomers = []
    with open(filname,'r') as fil:
        lines = [line.strip() for line in fil.readlines()]
    for i,line in enumerate(lines):
        try:
            n_atoms = int(line)
            isomers.append([li for li in lines[i:i+n_atoms+2]])
        except: pass
    return isomers

# e. list permutations
def list_permutations(input_list, indexes):
    perm_isom = []
    sublist = [input_list[i] for i in indexes]
    permlist = list(permutations(sublist))
    for perm in permlist:
        cp_input_list = [el for el in input_list]
        j = 0
        for i in range(len(cp_input_list)):
            if i in indexes: 
                cp_input_list[i] = perm[j]
                j += 1
        perm_isom.append(cp_input_list)
    return perm_isom

#f. compute distance, dist_mat e conn_mat
def conn_dist(isomers):
    dist_mat = []
    conn_mat = []
    coords_mat = []
    for isomer in isomers:
    # Create distance matrix and connectivity matrix (1|0 -> ye|no)
    # Both matrices are lists of lists of floats or booleans
    # [0] isomer i (list) 
    # [1] atom n of isomer i (list)
    # [2] distance atom x and n with x>n (float|bool)
        dist_line = []
        conn_line = []
        atoms = []
        coords = []
        for line in isomer[2:]:
            splitted = line.split()
            atoms.append(splitted[0])
            coords.append([float(x) for x in splitted[1:]])
        for i,line in enumerate(coords):
            distance = []
            for j in range(i+1,len(coords)):
                distance.append(((line[0]-coords[j][0])**2+(line[1]-coords[j][1])**2+(line[2]-coords[j][2])**2)**0.5)
            dist_line.append(distance)
            #print(distance)
    # Use largest caovalent radius (CC bond=0.75) *2 * 1.15 for safety
            conn_line.append(list(map(lambda x: x<1.5*1.15, distance)))
            #print(conn_line[-1])
        dist_mat.append(dist_line)
        conn_mat.append(conn_line)
        coords_mat.append(coords)
    return dist_mat,conn_mat,coords_mat

#g. filter GSM global list:
def filter_gsm(filname):
    with open(filname,'r') as f:
        matrix = [list(map(int, line.strip()[1:-1].split(', '))) for line in f.readlines()]

    seen_first_two = set()  # To keep track of the first two elements
    filtered_matrix = []  # To store the filtered list of lists

    for inner_list in matrix:
        first_two = tuple(inner_list[:2])
        if first_two not in seen_first_two:
            seen_first_two.add(first_two)
            filtered_matrix.append(inner_list)

    # Print the filtered matrix
    with open('gsm_filter.txt','w') as f:
        wr = csv.writer(f)
        wr.writerows(filtered_matrix)

#h. DFS algorithm
def dfs(graph, start, visited):
    visited.add(start)
    for neighbor in graph[start]:
        if neighbor not in visited:
            dfs(graph, neighbor, visited)

#i. Is single molecule?
def is_single_molecule(graph):
    visited = set()
    start_vertex = next(iter(graph.keys()))  # Start DFS from an arbitrary vertex
    dfs(graph, start_vertex, visited)
    return len(visited) == len(graph)

#############################

# Functions

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
def unite_xyz(md_path,md_name,crest_vers):
    isomers = []
    n_isom = 0
    filename = 'crest_products.xyz' # Output file CREST MD

    if float(crest_vers) >=3.0:
        os.system(f"cp {md_path}/{md_name} ./allcrestprods_unite.xyz")
        with open("allcrestprods_unite.xyz",'r') as f:
            lines = [line.strip() for line in f.readlines()]
        n_atoms = int(lines[0])
        # Write logfile here
        write_log("Crest version 3.0 was used..\n")
        write_log(f'n atoms: {n_atoms}')
        with open('natoms_unite.txt','w') as f: f.write(str(n_atoms))
        return

    if 'allcrestprods_unite.xyz' not in os.listdir():
    #if already run, don't recreate allprods (it takes a while...)
        md_list = [fold+'/' for fold in os.listdir(md_path) if md_name in fold]

        with open(md_path+md_list[0]+filename,'r') as f:
            lines = [line.strip() for line in f.readlines()]
        n_atoms = int(lines[0])
        # Write logfile here
        write_log(md_path+md_list[0]+filename)
        write_log(f'n atoms: {n_atoms}')
        with open('natoms_unite.txt','w') as f: f.write(str(n_atoms))

        for metadyn in md_list:
            write_log('Working on: '+metadyn)
            if filename not in os.listdir(md_path+metadyn): 
                write_log('Skipping, not found crest_products.xyz')
                continue
            else:
                with open(md_path+metadyn+filename,'r') as f:
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
        # No energy selection this time

# 2. Reorder isomers based on energy
def sort_xyz(filinp):
    isomers = []
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(filinp,n_atoms)
    # Sort by energy
    isomers.sort(key=lambda x:float(x[1].split()[0]))
    write_xyzlist(isomers,'allcrestprods_sort.xyz')

# 3. Select unique isomers using distance matrix
def bond_check(filinp,at_num):
    equal_pairs = set()
    isomers = []
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(filinp,n_atoms)

    ##CONN MAT QUI
    dist_mat,conn_mat,coords_mat = conn_dist(isomers)

    # Define equal pairs based on conn_mat 
    for i,el in enumerate(conn_mat):
        for j in range(i+1,len(conn_mat)):
            is_same = [ x for x in map(lambda x,y: x==y, el,conn_mat[j]) if x==0]
            if not is_same: equal_pairs.add((i,j))
            # Add +1 for correspondence with molden which starts from 1

    if list(equal_pairs) is not []: write_log(['%s %s\n' % x for x in equal_pairs])
    # Get rid of doubles (use set to include only uniques)
    to_be_removed = list(set([x for x in map(lambda x: x[1],equal_pairs)]))
    # Should already be good but let's sort it to make sure
    to_be_removed.sort(key=lambda x: x,reverse=1)
    write_log(['%s ' % x for x in to_be_removed])

    for el in to_be_removed: isomers.pop(el)
    write_xyzlist(isomers,'unique_bondcheck.xyz')
    with open('numisoms.txt','w') as f:
        f.write(str(len(isomers))+'\n')

    return dist_mat,conn_mat,coords_mat

# 3.5 Remove unique isomers that are actually fragments exploiting DFS algorithm
def remove_frags(filinp):
    write_log("\nEntering remove_frags function..\n")
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(filinp,n_atoms)
    dist_mat,conn_mat,coords_mat = conn_dist(isomers)
    to_be_removed=[]

    for i in range(len(conn_mat)):
        graph = {num:[] for num in range(n_atoms)}
        for counter,item in enumerate(conn_mat[i]):
            for ccounter,it in enumerate(item): 
                if it: 
                    graph[counter].append(counter+ccounter+1)
                    graph[counter+ccounter+1].append(counter)

        result = is_single_molecule(graph)
        if not result: to_be_removed.append(i)
        write_log(f"The given structure {i} is a single molecule: {result}\n")

    to_be_removed.sort(key=lambda x: x,reverse=1)
    print(f"{to_be_removed} - {len(isomers)}")
    for el in to_be_removed: isomers.pop(el)
    write_xyzlist(isomers,'nofrags_bondcheck.xyz')
    dist_mat,conn_mat,coords_mat = conn_dist(isomers)
    return dist_mat,conn_mat,coords_mat


# 4.V2 Select unique isomers using distance matrix AND permutations - GENERALIZED
def bond_check3(filinp):
    open('perm_equalpairs.txt','w').close()
    equal_pairs = set()
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(filinp,n_atoms)

    write_log('\nEntering permutations loop..:\n')
    for k,isomer in enumerate(isomers):
        prefix = [el for el in isomer[:2]]
        # Crea permutations of isomer K
        isom_k_permuts = gen_perm.generate_molecule_permutations(isomer[2:])

        for h,permutation in enumerate(isom_k_permuts):
            # Define useful variables
            isomers_copy = [el for el in isomers]
            isomers_copy[k] = prefix + permutation

            ##CONN MAT QUI
            _,conn_mat_copy,_ = conn_dist(isomers_copy)

            # Define equal pairs based on conn_mat
            for j in range(k+1,len(conn_mat_copy)):
                is_same = [ x for x in map(lambda x,y: x==y, conn_mat_copy[k],conn_mat_copy[j]) if x==0]
                if not is_same: 
                    equal_pairs.add((k,j))
                    with open('perm_equalpairs.txt','a') as f:
                        f.write(f'\nPermutation {h} - Eq_pair: {k} {j}')
                    write_log(f'\nPermutation {h} - Eq_pair: {k} {j}')
                # Add +1 for correspondence with molden which starts from 1

            #if list(equal_pairs) is not []: write_log(['%s %s\n' % x for x in equal_pairs])
            # Get rid of doubles (use set to include only uniques)
            to_be_removed = list(set([x for x in map(lambda x: x[1],equal_pairs)]))
            # Should already be good but let's sort it to make sure
            to_be_removed.sort(key=lambda x: x,reverse=1)
            #print(to_be_removed)

    write_log('Elements to be removed:\n')
    write_log(['%s ' % x for x in to_be_removed])
    for el in to_be_removed: 
        isomers.pop(el)
        write_log(f'\nIsomer {el} was successfully removed\n')
    write_log(f"\nCurrent number of isomers: {len(isomers)}\n")
    write_xyzlist(isomers,'perm_bondcheck.xyz')
    with open('numisoms.txt','w') as f:
        f.write(str(len(isomers))+'\n')

# 5. Select paths with permutations
def sel_paths(filinp):
    # Creates folder to store permutations of isomers
    if 'perm_folder' not in os.listdir(): os.mkdir('perm_folder')

    # Contains [isom K - isom I - permutC Z - permutH H]
    gsm_list = []
    gsm_included = set()
    open('selectedpaths.txt','w').close()
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(filinp,n_atoms)

    write_log('\nEntering permutations loop..:\n')
    for k,isomer in enumerate(isomers):
        prefix = [el for el in isomer[:2]]
        # Crea permutations of isomer K
        isom_k_permuts = gen_perm.generate_molecule_permutations(isomer[2:])

        for h,permutation in enumerate(isom_k_permuts):
            # Define useful variables
            isomers_copy = [el for el in isomers]
            isomers_copy[k] = prefix + permutation

            ##CONN MAT QUI
            _,conn_mat_copy,_ = conn_dist(isomers_copy)

            # Define number of broken / added bonds
            # "leaf from tree in forest for leaf in tree"
            k_bonds = [item for sublist in conn_mat_copy[k] for item in sublist]
            for i in range(k+1,len(conn_mat_copy)):
                curr_bonds = [item for sublist in conn_mat_copy[i] for item in sublist]
                is_same = [ x for x in map(lambda x,y: x==y, k_bonds,curr_bonds)]
                #print(i+k,is_same.count(False))
                if is_same.count(False)<3: 
                    gsm_list.append([k, i, h])
                    gsm_included.add(k)
                    gsm_included.add(i)

        # Writes permutations in file for isomer k
        write_log(f"Writing permutations of isomer {k} in perm_folder/permutations_{k}.xyz..\n")
        write_xyzlist([prefix+iso for iso in isom_k_permuts],f'perm_folder/permutations_{k}.xyz')

    with open(f'gsm_global.txt','w') as g:
        for lst in gsm_list: g.write(f'{lst}\n')

    write_log(str(list(gsm_included)))

# 6. setup_gsm from gsm_filter file and perm_bondcheck list of isomers
# NON TIENE LO STESSO ORDINE DELLE PERMUTAZIONI CAZZO - salvo gli xyz in bond_check2
def setup_gsm(gsminp,isomlist,model_gsm):

    # Get isomers list and n atoms
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(isomlist,n_atoms)

    # Get reaction list from gsm_filter.xyz
    with open(gsminp,'r') as f:
         reaction_list = list(csv.reader(f, delimiter=","))

    # create gsm folders
    for k,isomer in enumerate(isomers):
        # Setup GSM folder
        fold_name = f'{k}_gsm_fold'
        if not os.path.exists(fold_name): os.mkdir(fold_name)
        if not os.path.exists(f'{fold_name}/scratch'): os.mkdir(f'{fold_name}/scratch')
        gsm_path = os.getcwd()+f'/{fold_name}/'

    for k,isomer in enumerate(isomers):
        isom_k_perm = read_xyzlist(f'perm_folder/permutations_{k}.xyz',n_atoms)
        isomers_copy = [el for el in isomers] # Duplicate list for modfying it
        # Get GSM folder name and path
        fold_name = f'{k}_gsm_fold'
        gsm_path = os.getcwd()+f'/{fold_name}/'

        # Get list of reactions to study
        reaction_list_k = [el for el in reaction_list if el[0]==f'{k}']

        # Skip folder if reaction list empty
        if reaction_list_k == []: continue
        else: pass

        # Copy model data folder
        os.system(f'cp -r {model_gsm} {gsm_path}')

        # Write initial file for each gsm calculation
        for i,ind in enumerate(reaction_list_k):
            gsm_count = str(i+1).zfill(4)
            isom_k_num = int(ind[2])
            isom_k_perm[isom_k_num][1] = f'isom num {i} ' + isom_k_perm[isom_k_num][1]
            write_xyz(isom_k_perm[isom_k_num],f'{gsm_path}scratch/initial{gsm_count}.xyz')
            isomers_copy[int(ind[1])][1] = f'from perm {isom_k_num} of prev isomer ' + isomers_copy[int(ind[1])][1]
            append_xyz(isomers_copy[int(ind[1])],f'{gsm_path}scratch/initial{gsm_count}.xyz')

# 7. run_gsm()
def run_gsm():
    gsm_folds = [el for el in os.listdir() if 'gsm_fold' in el]
    gsm_folds.sort()

    # Alternative for cycles for running sets of subfolders
    #for fold in gsm_folds[-2:]:
    #for fold in gsm_folds[-4:-2]:
    #for fold in gsm_folds[-6:-4]:
    #for fold in gsm_folds[-7:-6]:
    #for fold in gsm_folds[-13:-11]:
    for fold in gsm_folds:
        if os.listdir(f'{fold}/scratch/') is not []:
            os.chdir(fold)
            inputss = [el for el in os.listdir('scratch') if el.startswith('initial')]
            for i in range(len(inputss)):
                print(f'rungsm {i+1} in {fold}')
                os.system(f'rungsm {i+1}')
            os.chdir('..')
        else: print(f'{fold} empty folder')


##################################


def main():
    at_num = {'1':'H','6':'C','8':'O','7':'N'}

    # Read input information
    with open("input.dat","r") as f:
        lines = [line.strip() for line in f.readlines() if (
                 line and not line.strip().startswith('#'))]
        
    config = {l.split('=')[0].strip():l.split('=')[1].strip() for l in lines}
    
    model_gsm = config['model_gsm']
    crest_version = config['crest_version']
    crest_md_path = config['crest_md_path']
    crest_out_name = config['crest_out_name']

    open('logfile.out','w').close() # Create logfile

 #   unite_molgen_xyz("C4H10_geo/",suffix="opt.xyz")

# #### 1 ###
    unite_xyz(crest_md_path,crest_out_name,crest_version) #Skips if allcrestprods_unite is already there; protocol depends on version of crest
#     sort_xyz('allcrestprods_unite.xyz') # -> allcrestprods_sort.xyz
# #### 2 ###
#     _,_,_ = bond_check('allcrestprods_sort.xyz',at_num) # -> unique_bondcheck.xyz
#     _,_,_ = remove_frags('unique_bondcheck.xyz') # -> nofrags_bondcheck.xyz
#     bond_check3('nofrags_bondcheck.xyz') # -> perm_bondcheck.xyz
#     sel_paths('perm_bondcheck.xyz')
# #### 3 ###
#     filter_gsm('gsm_global.txt')
#     setup_gsm('gsm_filter.txt','perm_bondcheck.xyz',model_gsm)
# #### 4 ###
#     run_gsm()

if __name__=="__main__":
    main()
