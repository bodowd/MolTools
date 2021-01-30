from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List
import itertools

def enumerate_mols(core: List[List[str]], chain: List[List[str]], rxn_smarts: str = '[c:1][#0].[#0][*:2]>>[c:1]-[*:2]' ):
    """

    :param core: List of lists [[ID, SMILES], [ID, SMILES],...]
    :param chain: List of lists [[ID, SMILES], [ID, SMILES],...]
    :param rxn_smarts:
    :return:
    """
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    res = []
    combos = itertools.product()
    for t in combos:
        cur_core_id = t[0][0]
        cur_core_smi = t[0][1]
        cur_core_mol = Chem.MolFromSmiles(cur_core_smi)

        cur_chain_id = t[1][0]
        cur_chain_smi = t[1][1]
        cur_chain_mol = Chem.MolFromSmiles(cur_chain_smi)
        p = rxn.RunReactants((cur_core_mol, cur_core_smi))
        product = Chem.MolToSmiles(p[0][0])

    res.append([cur_core_id, cur_core_mol, cur_chain_id, cur_chain_mol, product])
    return res

def replace_dummy_with_star(smi_list: List = None, mol_list: List = None, dummy: str = 'U'):
    if mol_list is not None:
        return [Chem.MolToSmiles(i).replace(dummy, '*') for i in mol_list]
    elif smi_list is not None:
        return [i.replace(dummy, '*') for i in smi_list]
    else:
        return 'Please input a list of smiles or list of mol objects'

