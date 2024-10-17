from multiprocessing import Pool
import itertools

import ioformat.pathtools as pa
import automol


def _geom_to_graph(args):
    geo, ene, i = args
    gra = automol.geom.graph(geo)
    gras =  automol.graph.connected_components(gra)
    if len(ene.split('=')) > 1:
        ene = ene.split('=')[1].split()[0]
    name = f'species_{i}'
    #name = automol.graph.smiles(gras[0])
    well_dct = {name: gras}
    return well_dct


def _names_to_reaction(args): #adl
    rxn_dct = {}
    rname, rgras, pname, pgras = args
    if len(rgras) + len(pgras) < 4:
        iso_dct, frm_bnd_lst, brk_bnd_lst = automol.reac.arbitrary_reactions(
            rgras, pgras)
        if iso_dct:
            rxn_dct[rname + ' = ' + pname] = (iso_dct, frm_bnd_lst, brk_bnd_lst)
    # isomorph = {i:i for i in range(14)}
    # rxn_dct = {"species_0 = species_1":(isomorph,frozenset([(3,8),(2,4)]),frozenset([(2,8),(3,4)]))}
    return rxn_dct


def finder(path): #adl
    traj_str = pa.read_file(path, 'allcrestprods_sort.xyz') #adl
    # input_str = pa.read_file('/lcrc/project/PACC/adl/qm/EXPLORER/test-sarah', 'input-stru.xyz') #adl
    geo_lst = automol.geom.from_xyz_trajectory_string(traj_str)
    idx_geo_lst = [(geo, ene, idx) for idx, (geo, ene) in enumerate(geo_lst)] #adl

    with Pool() as pool:
        all_well_dct = pool.map(_geom_to_graph, idx_geo_lst) # adl

    wells_dct = {}
    for well_dct in all_well_dct:
        wells_dct.update(well_dct)

    well_names = list(wells_dct.keys())
    rct_prd_pairs = itertools.combinations(well_names, 2)
    pair_checks = [
        (rname, wells_dct[rname], pname, wells_dct[pname])
        for rname, pname in rct_prd_pairs]
    
    with Pool() as pool:
        all_rxn_dct = pool.map(_names_to_reaction, pair_checks)

    rxns_dct = {}
    for rxn_dct in all_rxn_dct:
        rxns_dct.update(rxn_dct)
    
    geo_lst = [ geoi for geoi,_ in automol.geom.from_xyz_trajectory_string(traj_str)]

    return geo_lst, wells_dct, rxns_dct


if __name__ == "main":
    finder()
