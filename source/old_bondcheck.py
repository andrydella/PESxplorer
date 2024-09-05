from auxiliary_funcs import *
import generalized_permutations as gen_perm

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