"""
This module contains functions to parse a CIF file and extract the structure, models, chains, residues, and atoms.
"""

import os

import Bio.PDB.Structure
import gradio
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure

from .output import save_csv, save_gltf


def load_cif_file(file_name: str) -> tuple[Structure, str, list]:
    """
    Load a CIF file and extract the structure, structureID, and models.

    Args:
        file_name (str): The name of the CIF file to load.

    Returns:
        structure (Bio.PDB.Structure.Structure): The structure object.
        structureID (str): The structure ID.
        models (list): A list of model IDs.
    """
    # Create an MMCIFParser Object
    parser = MMCIFParser()

    # Extract the overall structure
    structure = parser.get_structure(
        os.path.splitext(os.path.basename(file_name))[0], file_name
    )  # Using inputs

    # read in the CIF file.
    pdb_info = MMCIF2Dict(file_name)

    # structureID = pdb_info["_pdbx_database_status.entry_id"][0]
    structure_id = pdb_info["_entry.id"][0]

    # Extract models
    models = sorted([model.id for model in structure])

    # Return stuff
    return structure, structure_id, models


def extract_chains_from_model(structure: Structure, selected_model: str) -> list:
    """
    Extract the chains from a selected model.

    Args:
        structure (Bio.PDB.Structure.Structure): The structure object.
        selected_model (str): The selected model ID.

    Returns:
        chains (list): A list of chain IDs.

    """
    chains = []
    for model in structure:
        if model.id == selected_model:
            chains = sorted([chain.id for chain in model])
    return chains


def extract_residues_from_chain(
    structure: Structure, selected_model: str, selected_chain: str
) -> list:
    """
    Extract the residues from a selected chain.

    Args:
        structure (Bio.PDB.Structure.Structure): The structure object.
        selected_model (str): The selected model ID.
        selected_chain (str): The selected chain ID.

    Returns:
        residues (list): A list of residue IDs.

    """
    residues = []
    for model in structure:
        if model.id == selected_model:
            for chain in model:
                if chain.id == selected_chain:
                    residues = sorted(
                        set(
                            (residue.resname, residue.id[1], residue.id[2])
                            for residue in chain.get_residues()
                        )
                    )
    return [f"{res[1]}: {res[0]}" for res in residues] if residues else []


def extract_chain_atoms(
    structure: Structure,
    selected_model: str,
    selected_chain: str,
    selected_residues: list[str],
) -> list[list[float]]:
    """
    Extract the atoms from a selected chain.

    Args:
        structure (Bio.PDB.Structure.Structure): The structure object.
        selected_model (str): The selected model ID.
        selected_chain (str): The selected chain ID.
        selected_residues (list): The selected residue IDs.

    Returns:
        atom_coords (list): A list of atom coordinates.
    """
    atom_coords = []
    selected_residue_numbers = {int(res.split(":")[0]) for res in selected_residues}
    for model in structure:
        if model.id == selected_model:
            for chain in model:
                if chain.id == selected_chain:
                    for residue in chain.get_residues():
                        if residue.id[1] in selected_residue_numbers:
                            for atom in residue.get_atoms():
                                atom_coords.append(
                                    [
                                        atom.element,
                                        atom.coord[0],
                                        atom.coord[1],
                                        atom.coord[2],
                                    ]
                                )
    return atom_coords


def rescale_dataset(atom_coords: list[list[float]]) -> list[list[float]]:
    """
    Rescale the atom coordinates to a smaller scale.

    Args:
        atom_coords (list): A list of atom coordinates.

    Returns:
        scaled_coords (list): A list of scaled atom coordinates.
    """

    scale_down_factor = 0.68

    # Do the operation for the desired residue, entity, model, etc...
    new_dataset = atom_coords
    for atom in new_dataset:
        if atom[0] == "N":
            atom.append(round(1.55 * scale_down_factor, 1))  # add radius
            atom.pop(0)  # delete first item in the list (atom symbol)
        elif atom[0] == "C":
            atom.append(round(1.7 * scale_down_factor, 1))
            atom.pop(0)
        elif atom[0] == "O":
            atom.append(round(1.52 * scale_down_factor, 1))
            atom.pop(0)
        elif atom[0] == "H":
            atom.append(round(1.10 * scale_down_factor, 2))
            atom.pop(0)
        elif atom[0] == "P":
            atom.append(round(1.8 * scale_down_factor, 1))
            atom.pop(0)
        else:  # By default all other atoms not listed above are given a radius equal to a carbon atom
            atom.append(round(1.55 * scale_down_factor, 1))
            atom.pop(0)

    # Scale the values.
    scaled_coords = [[round(entry * 10, 3) for entry in atom] for atom in new_dataset]

    return scaled_coords


def process_cif(
    file: gradio.File,
    selected_model: str,
    selected_chain: str,
    selected_residues: list[str],
) -> tuple[gradio.File, gradio.Model3D]:
    """
    Process the CIF file and save the scaled dataset as a CSV file.

    Args:
        file (gradio.File): The CIF file to process.
        selected_model (str): The selected model ID.
        selected_chain (str): The selected chain ID.
        selected_residues (list): The selected residue IDs.

    Returns:
        file_output (gradio.File): The CSV file to download.
        model_preview (gradio.Model3D): The model preview.
    """
    if file is None:
        return gradio.File(
            None, label="Download CSV", interactive=False
        ), gradio.Model3D(None, label="Model Preview", interactive=False)
    else:
        structure, structure_id, _ = load_cif_file(file.name)
        atom_coords = extract_chain_atoms(
            structure, selected_model, selected_chain, selected_residues
        )
        if not atom_coords:
            return gradio.File(
                None, label="Download CSV", interactive=False
            ), gradio.Model3D(None, label="Model Preview", interactive=False)

        scaled_dataset = rescale_dataset(atom_coords)
        csv_name = (
            f"{structure_id}_Model_{selected_model}_Chain_{selected_chain}_Residues.csv"
        )
        gltf_name = csv_name.replace(".csv", ".obj")
        save_csv(csv_name, scaled_dataset)
        save_gltf(gltf_name, scaled_dataset)

        return gradio.File(
            csv_name, label=f"Download {csv_name}", interactive=True
        ), gradio.Model3D(gltf_name, label="Model Preview")
