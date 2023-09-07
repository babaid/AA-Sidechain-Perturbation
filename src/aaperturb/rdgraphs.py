from typing import Union

import rdkit
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from .constants import *
from loguru import logger as log


def create_edge_index(mol: rdkit.Chem.rdchem.Mol) -> np.ndarray:
    begin_atoms = [a.GetBeginAtomIdx() for a in mol.GetBonds()]
    end_atoms = [a.GetEndAtomIdx() for a in mol.GetBonds()]
    return np.array([begin_atoms, end_atoms])

def get_coords(mol: rdkit.Chem.rdchem.Mol) -> np.ndarray:
    for c in mol.GetConformers():
        positions = np.array(c.GetPositions())
    return positions

def set_coords(mol: rdkit.Chem.rdchem.Mol, coords, idxs):
    conf = mol.GetConformer()
    for i, idx in enumerate(idxs):
        x,y,z = coords[idx, :].tolist()
        conf.SetAtomPosition(idx, rdkit.Geometry.rdGeometry.Point3D(x, y, z))

def get_bond_labels(mol: rdkit.Chem.rdchem.Mol, onehot: bool = False) -> np.ndarray:
    bond_labels = []
    for bond in mol.GetBonds():
        bond_labels.append(BOND_TYPES_RDKIT[str(bond.GetBondType())])
    return np.array(bond_labels)


def get_residue(mol: rdkit.Chem.rdchem.Mol, chain: str, resnum: int, exclude_atom_types: list[str] = ["N", "CA", "C", "O"]):
    residue = []
    atom_idxs = []
    if exclude_atom_types == []:
        for i, atom in enumerate(mol.GetAtoms()):
            resinfo = atom.GetPDBResidueInfo()
            # find residue
            if resinfo.GetChainId() == chain and resinfo.GetResidueNumber() == resnum:
                resname = resinfo.GetResidueName()
                if resinfo.GetName().strip() not in exclude_atom_types:
                    residue.append(atom)
                    atom_idxs.append(i)
    else:
        for i, atom in enumerate(mol.GetAtoms()):
            resinfo = atom.GetPDBResidueInfo()
            # find residue
            if resinfo.GetChainId() == chain and resinfo.GetResidueNumber() == resnum:
                    #residue.append(atom)
                    atom_idxs.append(i)
    return atom_idxs, resname

def get_residue_for_rotation(mol: rdkit.Chem.rdchem.Mol, chain: str, resnum: int, include_lst):
    residue = {}
    atom_idxs = []
    idx_dict = {}
    for i, atom in enumerate(mol.GetAtoms()):
            resinfo = atom.GetPDBResidueInfo()
            # find residue
            if resinfo.GetChainId() == chain and resinfo.GetResidueNumber() == resnum:
                resname = resinfo.GetResidueName()
                if resinfo.GetName().strip() in include_lst:      
                    #residue[resinfo.GetName().strip()] = atom
                    atom_idxs.append(i)
                idx_dict[resinfo.GetName().strip()] = i
    return idx_dict, atom_idxs, resname


def no_clashes(coords: np.ndarray, tolerance: int = 0):
    distmat = cdist(coords, coords)
    np.fill_diagonal(distmat, 1.0)
    n_clashes =  np.where(distmat < 0.5, 1, 0).sum()
    if n_clashes <= tolerance:
        return True
    return False
