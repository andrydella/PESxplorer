import itertools
import numpy as np
from collections import defaultdict
from typing import List, Tuple
from copy import deepcopy

# Step 1: Parse the input string
def parse_molden_input(molden_str: list) -> List[Tuple[int, str, List[str]]]:

    # Parse the atom type and coordinates, preserving the original line index
    atom_lines = []
    for idx, line in enumerate(molden_str):
        parts = line.split()
        atom_type = parts[0]
        coordinates = parts[1:]
        atom_lines.append((idx, atom_type, coordinates))
    
    return atom_lines

# Step 2: Group atom lines by type, but keep the original indices
def group_by_atom_type(atom_lines: List[Tuple[int, str, List[str]]]) -> dict:
    grouped = defaultdict(list)
    for idx, atom_type, coordinates in atom_lines:
        grouped[atom_type].append((idx, atom_type, coordinates))
    return grouped

# Step 3:
def prune_groups(upper_triangular,grouped):
    def expand_conn_mat(upper_triangular):
        # Determine the size of the matrix
        n = len(upper_triangular)
        # Create an n x n matrix initialized with False
        matrix = np.full((n, n), False, dtype=bool)
        # Populate the upper triangular part (excluding the diagonal)
        for i, row in enumerate(upper_triangular):
            matrix[i, i+1:] = row
        # Mirror the upper triangular part into the lower triangular part
        matrix = matrix | matrix.T
        return matrix
    
    conn_mat = expand_conn_mat(upper_triangular)
    groups = {}
    for atom_type, atoms in grouped.items():
        groups[atom_type] = []
        indices = [atom[0] for atom in atoms]
        not_yet_compared = set(indices)

        while not_yet_compared:
            idx = not_yet_compared.pop()
            group = {idx}

            for i in list(not_yet_compared):
                if np.array_equal(conn_mat[idx],conn_mat[i]):
                    group.add(i)
                    not_yet_compared.remove(i)
            groups[atom_type].append(group)

    grouped_copy = defaultdict(list)
    removed_count = defaultdict(int)
    for atom_type, elements in groups.items():
        for group in elements:
            min_index = min(group)
            removed_count[atom_type] += len(group)-1
            for item in grouped[atom_type]:
                if item[0] == min_index:
                    grouped_copy[atom_type].append(item)
                    break

    return grouped_copy,removed_count,groups

# Step 4: Generate permutations for each group of atoms
def generate_permutations(grouped: dict) -> dict:
    perms = {}
    for atom_type, atoms in grouped.items():
        indices = [atom[0] for atom in atoms]
        # Generate all permutations of the indices
        perms[atom_type] = list(itertools.permutations(indices))

    return perms

# Step 5: Combine permutations to create all possible configurations
def combine_permutations(perms: dict, atom_lines: List[Tuple[int, str, List[str]]], added_perms) -> List[List[Tuple[int, str, List[str]]]]:
    
    combined_permutations = list(itertools.product(*perms.values()))

    # get reference atom positions
    ref_lines = combined_permutations[0]
    all_permutations = []

    for perm in combined_permutations:
        current_perm_on_atom = [None]*len(atom_lines)
        for i,indices in enumerate(perm):
            for j,idx in enumerate(indices):
                current_perm_on_atom[ref_lines[i][j]] = (
                    ref_lines[i][j],atom_lines[idx][1],atom_lines[idx][2])
        for i,element in enumerate(current_perm_on_atom):
            if element is None:
                current_perm_on_atom[i] = atom_lines[i]
        
        #Sort by new index to ensure correct ordering
        new_lines = sorted(current_perm_on_atom, key=lambda x: x[0])
        all_permutations.append(new_lines)

    ref_stru = all_permutations[0]
    for perm in added_perms:
        for idx1,idx2 in perm:
            new_stru = deepcopy(ref_stru)
            new_stru[idx1],new_stru[idx2] = ref_stru[idx2],ref_stru[idx1]
        all_permutations.append(new_stru)

    return all_permutations

# Format output into a list of strings with "Atom coord1 coord2 coord3"
def format_permutations(permutations: List[List[Tuple[int, str, List[str]]]]) -> List[List[str]]:
    formatted_results = []
    for perm in permutations:
        formatted = [f"{atom_type} {' '.join(coords)}" for _,atom_type,coords in perm]
        formatted_results.append(formatted)
    return formatted_results


# Step 6: Generate all molecule permutations with the given input
def generate_molecule_permutations(molden_str: list, up_tr: List[List[bool]]) -> List[List[Tuple[int, str, List[str]]]]:
    atom_lines = parse_molden_input(molden_str)
    
    # Group by atom type, retaining the original indices
    grouped = group_by_atom_type(atom_lines)

    grouped,removed_count,groups = prune_groups(up_tr,grouped)
    
    extra_perms = add_structures(groups)

    # Generate permutations for the indices within each group
    perms = generate_permutations(grouped)
    
    # Combine these permutations to create all possible configurations
    all_permutations = combine_permutations(perms, atom_lines,extra_perms)
    
    return format_permutations(all_permutations)


#Add one permutation of excluded equal atoms
def add_structures(data):
    output = defaultdict(list)
    # Iterate over each key in the data
    for key, groupss in data.items():
        # Identify groups with more than one element
        for i, group in enumerate(groupss):
            if len(group) > 1:              
                for j, other_group in enumerate(groupss):
                    if i != j:  
                        for elem in sorted(list(group))[1:]:
                            output[key].append((elem,min(other_group))) 
    extra_perms = list(itertools.product(*output.values()))
    return extra_perms

# def add_structures(data):
#     output = []
#     # Iterate over each key in the data
#     for key, groupss in data.items():
#         # Identify groups with more than one element
#         for i, group in enumerate(groupss):
#             if len(group) > 1:
#                 first_element = min(group)
                
#                 # Iterate over the rest of the groups to create substitutions
#                 for j, other_group in enumerate(groupss):
#                     if i != j:  # Avoid substituting within the same group
#                         for elem in sorted(list(group))[1:]:  # Skip the first element
#                             new_data = deepcopy(data)
#                             # Substitute the element with the first element of the other group
#                             new_data[key][i].remove(elem)
#                             new_data[key][i].add(min(other_group))
#                             new_data[key][j].remove(min(other_group))
#                             new_data[key][j].add(elem)
#                             output.append(new_data) 
#     more_perms=defaultdict(list)
#     for variation in output:
#         for atom_type,group in variation.items():
#             new_permut = []
#             for subset in group:
#                 new_permut.extend(sorted(list(subset)))
#             more_perms[atom_type].append(new_permut)

#     added_perms = []
#     # Soluzione corrente è comoda per C4H10 perché i primi 4 atomi sono C e poi ho gli H. Se
#     # gli indici sono mischiati la ricostruzione non è cosi semplice e va ripensata
#     for i in range(len(output)):
#         pass
#     return more_perms




def main():
    # Example input
    molden_str = [
        "C 0. 0. 0.",
        "N 5. 0. 0.",
        "C 2. 0. 0.",
        "N 7. 1. 0.",
        "N 1. 1. 0.",
        "C 9. 1. 0.",
    ]

    # Generate all permutations for the molecule with new line arrangements
    permutations = generate_molecule_permutations(molden_str)
    # Comment prune_groups to run

    # Display the results
    for i, perm in enumerate(permutations):
        print(f"Permutation {i + 1}:")
        print(perm)

if __name__=="__main__":
    main()

