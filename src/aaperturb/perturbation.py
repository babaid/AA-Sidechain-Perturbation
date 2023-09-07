import rdkit
import numpy as np
from scipy.spatial.transform import Rotation
from .rdgraphs import *
from .constants import *

from loguru import logger as log


class AAPerturbator(object):

    perturbationcounter = 0
    def __init__(self, pdb_path) -> None:
        self.mol = rdkit.Chem.rdmolfiles.MolFromPDBFile(pdb_path)
        self.coords = get_coords(self.mol)
        self.ref_coords = self.coords.copy()
        self.rmsd = 0.0

    def save(self, path: str):
        if self.perturbationcounter > 0:
            rdkit.Chem.rdmolfiles.MolToPDBFile(self.mol, path)
        else:
            print("It seems there was no perturbation, so we wont save this.")


    def gradual_sidechain_perturbation(self, chain: str, resnum: str, perturbation_rmsd: float, angle_step_size: float = 1e-1):
        self.perturbationcounter+=1
        temp_coords = self.coords.copy()
        full_rot_idxs, resname = get_residue(self.mol, chain, resnum, [])
        rmsd = 0.0
        old_rmsd = 0.0
        rmsd_period_cntr = 0
        rmsd_patience = 3
        if resname not in ["GLY", "PRO"]:
            while rmsd<perturbation_rmsd:
                for ax in ROT_AXES[resname]:
                    log.info(f"Rotating around axis {ax}")
                    axis_pivot_idx = AMINO_LST[resname].index(ax[1])
                    connected_atoms = AMINO_LST[resname][axis_pivot_idx:]
                    log.info(f"The substructure contains following atoms: {connected_atoms}")
                    idx_dict, rot_idxs, resname= get_residue_for_rotation(self.mol, chain, resnum, include_lst = connected_atoms)
                    axis = self.coords[idx_dict[ax[1]], :] #- coords[idx_dict["CB"], :]
                    axis = axis/np.linalg.norm(axis)
                    R_CA_CB = Rotation.from_rotvec(angle_step_size * axis, degrees=True)
                    temp_coords[rot_idxs, :] = R_CA_CB.apply(self.coords[rot_idxs, :]).copy()
                    if no_clashes(temp_coords):
                        self.coords[rot_idxs, :] = temp_coords[rot_idxs, :].copy()
                        log.info("There were no clashes between atoms.")
                    else:
                        log.info("Some atoms clashed. Chill. Will be resolved.")
                rmsd = np.sqrt(np.mean((self.ref_coords[full_rot_idxs, :]-self.coords[full_rot_idxs, :])**2))
                if rmsd < old_rmsd:
                    log.info("Reached maximum possible RMSD on this rotation axis.")
                    rmsd_period_cntr+=1
                else:
                    old_rmsd = rmsd
                if rmsd_period_cntr > rmsd_patience:
                    log.info("Termination due to RMSD periodicity patience timeout.")
                    #return coords, rot_idxs
                log.info(f"Residue RMSD: {rmsd}")
            set_coords(self.mol, self.coords, full_rot_idxs) 
        else:
            log.info(f"No rotation possible for {resname}, consider mutating it or just dont do anything!")
            return
        
    def random_perturbaton(self, chain: str, resnum: str, perturbation_rmsd: float):
        self.perturbationcounter+=1
        temp_coords = self.coords.copy()
        full_rot_idxs, resname = get_residue(self.mol, chain, resnum, [])
        rmsd = 0.0
        old_rmsd = 0.0
        rmsd_period_cntr = 0
        rmsd_patience = 3
        if resname not in ["GLY", "PRO"]:
            while rmsd<perturbation_rmsd:
                for ax in ROT_AXES[resname]:
                    log.info(f"Rotating around axis {ax}")
                    axis_pivot_idx = AMINO_LST[resname].index(ax[1])
                    connected_atoms = AMINO_LST[resname][axis_pivot_idx:]
                    log.info(f"The substructure contains following atoms: {connected_atoms}")
                    idx_dict, rot_idxs, resname= get_residue_for_rotation(self.mol, chain, resnum, include_lst = connected_atoms)
                    axis = self.coords[idx_dict[ax[1]], :] #- coords[idx_dict["CB"], :]
                    axis = axis/np.linalg.norm(axis)
                    R_CA_CB = Rotation.from_rotvec(np.random.uniform(0, 180)*axis, degrees=True)
                    temp_coords[rot_idxs, :] = R_CA_CB.apply(self.coords[rot_idxs, :]).copy()
                    if no_clashes(temp_coords):
                        self.coords[rot_idxs, :] = temp_coords[rot_idxs, :].copy()
                        log.info("There were no clashes between atoms.")
                    else:
                        log.info("Some atoms clashed. Chill. Will be resolved.")
                rmsd = np.sqrt(np.mean((self.ref_coords[full_rot_idxs, :]-self.coords[full_rot_idxs, :])**2))
                if rmsd < old_rmsd:
                    log.info("Reached maximum possible RMSD on this rotation axis.")
                    rmsd_period_cntr+=1
                else:
                    old_rmsd = rmsd
                if rmsd_period_cntr > rmsd_patience:
                    log.info("Termination due to RMSD periodicity patience timeout.")
                    #return coords, rot_idxs
                log.info(f"Residue RMSD: {rmsd}")
            set_coords(self.mol, self.coords, full_rot_idxs) 
        else:
            log.info(f"No rotation possible for {resname}, consider mutating it or just dont do anything!")
            return

    
"""
def gradual_sidechain_perturbation(mol: rdkit.Chem.rdchem.Mol, chain: str, resnum: int, perturbation_rmsd: float, angle_step_size: float):
    
    Creates gradual rotation around the possible angles of an amino acid sidechain until a certain rmsd is reached.
    
    coords = get_coords(mol).numpy()
    ref_coords = coords.copy()
    temp_coords = coords.copy()
    full_rot_idxs, resname = get_residue(mol, chain, resnum, [])
    rmsd = 0.0
    old_rmsd = 0.0
    rmsd_period_cntr = 0
    rmsd_patience = 3
    
    if resname not in ["GLY", "PRO"]:
        while rmsd<perturbation_rmsd:
            for ax in ROT_AXES[resname]:
                log.info(f"Rotating around axis {ax}")
                axis_pivot_idx = AMINO_LST[resname].index(ax[1])
                connected_atoms = AMINO_LST[resname][axis_pivot_idx:]
                log.info(f"The substructure contains following atoms: {connected_atoms}")
                idx_dict, rot_idxs, resname= get_residue_for_rotation(mol, chain, resnum, include_lst = connected_atoms)
                axis = coords[idx_dict[ax[1]], :] #- coords[idx_dict["CB"], :]
                axis = axis/np.linalg.norm(axis)
                R_CA_CB = Rotation.from_rotvec(axis, angle_step_size)
                temp_coords[rot_idxs, :] = R_CA_CB.apply(coords[rot_idxs, :]).copy()
                if no_clashes(temp_coords):
                    coords[rot_idxs, :] = temp_coords[rot_idxs, :].copy()
                    log.info("There were no clashes between atoms.")
                else:
                    log.info("Some atoms clashed. Chill. Will be resolved.")
            rmsd = np.sqrt(np.mean((ref_coords[full_rot_idxs, :]-coords[full_rot_idxs, :])**2))
            if rmsd < old_rmsd:
                log.info("Reached maximum possible RMSD on this rotation axis.")
                rmsd_period_cntr+=1
            else:
                old_rmsd = rmsd
            if rmsd_period_cntr > rmsd_patience:
                log.info("Termination due to RMSD periodicity patience timeout.")
                #return coords, rot_idxs
            log.info(f"Residue RMSD: {rmsd}")
        return coords, rot_idxs
    else:
        log.info(f"No rotation possible for {resname}")
        return coords, rot_idxs

"""