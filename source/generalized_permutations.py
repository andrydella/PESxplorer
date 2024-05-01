import itertools
from collections import defaultdict
from typing import List, Tuple

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

# Step 3: Generate permutations for each group of atoms
def generate_permutations(grouped: dict) -> dict:
    perms = {}
    for atom_type, atoms in grouped.items():
        indices = [atom[0] for atom in atoms]
        # Generate all permutations of the indices
        perms[atom_type] = list(itertools.permutations(indices))
    
    return perms

# Step 4: Combine permutations to create all possible configurations
def combine_permutations(perms: dict, atom_lines: List[Tuple[int, str, List[str]]]) -> List[List[Tuple[int, str, List[str]]]]:
    combined_permutations = list(itertools.product(*perms.values()))

# get reference atom positions
    ref_lines = combined_permutations[0]
    #print(ref_lines)

    #print(combined_permutations)

#    print(atom_lines)
    all_permutations = []
    for perm in combined_permutations:
        current_perm_on_atom = [None]*len(atom_lines)
        for i,indices in enumerate(perm):
            for j,idx in enumerate(indices):
                current_perm_on_atom[ref_lines[i][j]] = (
                    ref_lines[i][j],atom_lines[idx][1],atom_lines[idx][2])
        
        #Sort by new index to ensure correct ordering
        new_lines = sorted(current_perm_on_atom, key=lambda x: x[0])
        all_permutations.append(new_lines)

    
    return all_permutations


# Format output into a list of strings with "Atom coord1 coord2 coord3"
def format_permutations(permutations: List[List[Tuple[int, str, List[str]]]]) -> List[List[str]]:
    formatted_results = []
    for perm in permutations:
        formatted = [f"{atom_type} {' '.join(coordinates)}" for _, atom_type, coordinates in perm]
        formatted_results.append(formatted)
    return formatted_results

# Step 5: Generate all molecule permutations with the given input
def generate_molecule_permutations(molden_str: list) -> List[List[Tuple[int, str, List[str]]]]:
    atom_lines = parse_molden_input(molden_str)
    
    # Group by atom type, retaining the original indices
    grouped = group_by_atom_type(atom_lines)
    
    # Generate permutations for the indices within each group
    perms = generate_permutations(grouped)
    
    # Combine these permutations to create all possible configurations
    all_permutations = combine_permutations(perms, atom_lines)
    
    return format_permutations(all_permutations)



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
    #print(permutations)


    # Display the results
    for i, perm in enumerate(permutations):
        print(f"Permutation {i + 1}:")
        print(perm)

if __name__=="__main__":
    main()

