from auxiliary_funcs import *
import generalized_permutations as gen_perm
import os

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
            # Mentally add +1 for correspondence with molden which starts from 1

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
    _,conn_mat,_ = conn_dist(isomers)

    write_log('\nEntering permutations loop..:\n')
    for k,isomer in enumerate(isomers):
        write_log(f'Working on isomer {k}\n')
        prefix = [el for el in isomer[:2]]
        # Crea permutations of isomer K
        isom_k_permuts = gen_perm.generate_molecule_permutations(isomer[2:]) #,conn_mat[k]
        write_log(f'Number of permutations {len(isom_k_permuts)}\n')

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
                # Mentally add +1 for correspondence with molden which starts from 1

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
    _,conn_mat,_ = conn_dist(isomers)

    write_log('\nEntering permutations loop..:\n')
    for k,isomer in enumerate(isomers):
        write_log(f'Working on isomer {k}\n')
        prefix = [el for el in isomer[:2]]
        # Crea permutations of isomer K
        isom_k_permuts = gen_perm.generate_molecule_permutations(isomer[2:]) #,conn_mat[k]
        write_log(f'Number of permutations {len(isom_k_permuts)}\n')

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
                if is_same.count(False)<5: # Up to break2form2
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
# NON TIENE LO STESSO ORDINE DELLE PERMUTAZIONI - salvo gli xyz in bond_check3
def setup_gsm(gsminp,isomlist,model_gsm):

    # Get isomers list and n atoms
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(isomlist,n_atoms)

    # Get reaction list from gsm_filter.xyz
    with open(gsminp,'r') as f:
         reaction_list = list(csv.reader(f, delimiter=","))

    # create gsm folders
    os.system('mkdir -p GSM_FOLDS')
    for k,isomer in enumerate(isomers):
        # Setup GSM folder
        fold_name = f'GSM_FOLDS/{k}_gsm_fold'
        if not os.path.exists(fold_name): os.mkdir(fold_name)
        if not os.path.exists(f'{fold_name}/scratch'): os.mkdir(f'{fold_name}/scratch')
        gsm_path = os.getcwd()+f'/{fold_name}/'

    for k,isomer in enumerate(isomers):
        isom_k_perm = read_xyzlist(f'perm_folder/permutations_{k}.xyz',n_atoms)
        isomers_copy = [el for el in isomers] # Duplicate list for modfying it
        # Get GSM folder name and path
        fold_name = f'GSM_FOLDS/{k}_gsm_fold'
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
    
# 6. setup_gsm from gsm_filter file and perm_bondcheck list of isomers
# NON TIENE LO STESSO ORDINE DELLE PERMUTAZIONI - salvo gli xyz in bond_check3
def setup_gsm(gsminp,isomlist,model_gsm):

    # Get isomers list and n atoms
    with open('natoms_unite.txt','r') as f:
        n_atoms = int(f.read())
    isomers = read_xyzlist(isomlist,n_atoms)

    # Get reaction list from gsm_filter.xyz
    with open(gsminp,'r') as f:
         reaction_list = list(csv.reader(f, delimiter=","))

    # create gsm folders
    os.system('mkdir -p GSM_FOLDS')
    for k,isomer in enumerate(isomers):
        # Setup GSM folder
        fold_name = f'GSM_FOLDS/{k}_gsm_fold'
        if not os.path.exists(fold_name): os.mkdir(fold_name)
        if not os.path.exists(f'{fold_name}/scratch'): os.mkdir(f'{fold_name}/scratch')
        gsm_path = os.getcwd()+f'/{fold_name}/'

    for k,isomer in enumerate(isomers):
        isom_k_perm = read_xyzlist(f'perm_folder/permutations_{k}.xyz',n_atoms)
        isomers_copy = [el for el in isomers] # Duplicate list for modfying it
        # Get GSM folder name and path
        fold_name = f'GSM_FOLDS/{k}_gsm_fold'
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