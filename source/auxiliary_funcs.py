from itertools import permutations
import pickle
import csv

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

#j. Interface with automol
def amech_data_structure(filinp):
    import automol
    with open(filinp,'r') as f:
        trajec_string = f.read()
    # Probably it will be the opposite, I will start from well_dct from Sarah 
    # and create geos list
    geos = [ geoi for geoi,_ in automol.geom.from_xyz_trajectory_string(trajec_string)]
    well_dct = {f"species_{i}":automol.geom.graph(geoi) for i,geoi in enumerate(geos)}
    print("Here Sarah looks for possible reactivity (assume b2f2)")
    print("Make fake reaction dct")
    isomorph = {i:i for i in range(len(geos[0]))}
    isomorph[2] = 3
    isomorph[3] = 2
    reac_dct = {"species_0 = species_2":(isomorph,frozenset([(3,8),(2,4)]),frozenset([(2,8),(3,4)]))}
    print(reac_dct)
    geos.append(automol.geom.swap_coordinates(geos[1], 2, 3))
    well_dct["species_2"] = automol.geom.graph(geos[2])
    print("I inverted position of atoms 2 and 3 to try switching them back again, so geos2 is same as 1 but inverted atoms")

    return geos,well_dct,reac_dct

# k. Swap atoms in automech geometry
def swap_atoms(atom_order,geo_p):
    import automol

    visited_ats = set()
    for at1,at2 in atom_order.items():
        if at1 not in visited_ats:
            if at1 != at2:
                print(f"Swapping ats {at1} {at2}")
                geo_p = automol.geom.swap_coordinates(geo_p, at1, at2)
            visited_ats.update([at1,at2])

    return geo_p

# l.
def write_pickle(data,name):
    with open(f"{name}.pickle","wb") as f:
        pickle.dump(data,f,pickle.HIGHEST_PROTOCOL)

# m.
def read_pickle(name):
    with open(name + '.pickle', 'rb') as f:
        return pickle.load(f)