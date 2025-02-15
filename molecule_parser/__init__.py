import csv
import os

import Bio.PDB.Structure
import gradio
import trimesh
import trimesh.visual
from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure


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


def save_csv(csv_name: str, scaled_dataset: list[list[float]]):
    """
    Save the scaled dataset to a CSV file.

    Args:
        csv_name (str): The name of the CSV file to save.
        scaled_dataset (list): The scaled dataset to save.

    Return:
        None
    """
    header = ["X_Coord", "Y_Coord", "Z_Coord", "Atomic_Radius"]

    with open(csv_name, "w") as f:
        write = csv.writer(f)

        write.writerow(header)
        write.writerows(scaled_dataset)


def save_gltf(model_name: str, atom_coordiantes: list[list[float]]):
    """
    Save the scaled dataset as a glTF file.

    Args:
        model_name (str): The name of the glTF file to save.
        atom_coordiantes (list): The scaled dataset to save.

    Return:
        None

    """
    # Create an empty scene
    scene = trimesh.Scene()

    # For each atom, create a sphere and translate it to the correct position.
    for i, atom in enumerate(atom_coordiantes):
        # Get info for the atom
        x, y, z, radius = atom

        # Create a sphere using an icosphere
        sphere = trimesh.creation.icosphere(subdivisions=1, radius=radius)

        # Move the sphere to the atom's location
        sphere.apply_translation([x, y, z])

        # Add the sphere to the scene with a unique node name
        scene.add_geometry(sphere, node_name=f"sphere_{i}")

    # Export the entire scene as a glTF file.
    mesh = scene.to_geometry()
    mesh.export(model_name)


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


def update_model_dropdown_from_file(file: gradio.File) -> dict[str, any]:
    """
    Update the model dropdown from a CIF file.

    Args:
        file (gradio.File): The CIF file to process.

    Returns:
        model_dropdown (gradio.Dropdown): The model dropdown.

    """
    if file is None:
        models = []
    else:
        _, _, models = load_cif_file(file.name)
    return gradio.update(choices=models, value=models[0] if models else None)


def update_model_dropdown_from_code(code: str) -> tuple[dict[str, any], dict[str, any]]:
    """
    Update the model dropdown from a PDB code.

    Args:
        code (str): The PDB code to process.

    Returns:
        model_dropdown (gradio.Dropdown): The model dropdown.
        file_input (gradio.File): The CIF file to process.
    """
    if len(code) < 4:
        pdb_filename = None
        models = []
    else:

        # Create an instance of the PDBList class
        pdb_list = PDBList()

        # Download the MMCIF file using the retrieve_pdb_file method
        pdb_filename = pdb_list.retrieve_pdb_file(
            code, pdir="data", file_format="mmCif"
        )

        _, _, models = load_cif_file(pdb_filename)

    return gradio.update(
        choices=models, value=models[0] if models else None
    ), gradio.update(value=pdb_filename)


def update_chain_dropdown(file: gradio.File, selected_model: str) -> dict[str, any]:
    """
    Update the chain dropdown from a CIF file.

    Args:
        file (gradio.File): The CIF file to process.
        selected_model (str): The selected model ID.

    Returns:
        chain_dropdown (gradio.Dropdown): The chain dropdown

    """
    if file is None:
        chains = []
    else:
        structure, _, _ = load_cif_file(file.name)
        chains = extract_chains_from_model(structure, selected_model)
    return gradio.update(choices=chains, value=chains[0] if chains else None)


def update_residue_dropdown(
    file: gradio.File, selected_model: str, selected_chain: str
) -> dict[str, any]:
    """
    Update the residue dropdown from a CIF file.

    Args:
        file (gradio.File): The CIF file to process.
        selected_model (str): The selected model ID.
        selected_chain (str): The selected chain ID.

    Returns:
        residue_dropdown (gradio.Dropdown): The residue dropdown.
    """

    if file is None:
        residues = []
    else:
        structure, _, _ = load_cif_file(file.name)
        residues = extract_residues_from_chain(
            structure, selected_model, selected_chain
        )
    return gradio.update(choices=residues, value=residues if residues else [])


def gui() -> gradio.Blocks:
    """
    Create the GUI for the molecule parser.

    Returns:
        demo (gradio.Blocks): The GUI for the molecule parser.
    """
    with gradio.Blocks() as demo:
        with gradio.Row():
            with gradio.Column():
                with gradio.Row():
                    file_input = gradio.File(label="Upload CIF File")
                    code_dropdown = gradio.Textbox(label="Enter PDB Code")

                model_dropdown = gradio.Dropdown(
                    label="Select Model", choices=[], interactive=True
                )
                chain_dropdown = gradio.Dropdown(
                    label="Select Chain", choices=[], interactive=True
                )
                residue_dropdown = gradio.Dropdown(
                    label="Select Residues",
                    choices=[],
                    interactive=True,
                    multiselect=True,
                )
            with gradio.Column():
                model_preview = gradio.Model3D(label="Model Preview", interactive=False)
                file_output = gradio.File(label="Download CSV", interactive=False)

        file_input.change(
            update_model_dropdown_from_file,
            inputs=[file_input],
            outputs=[model_dropdown],
        )
        code_dropdown.change(
            update_model_dropdown_from_code,
            inputs=[code_dropdown],
            outputs=[model_dropdown, file_input],
        )
        model_dropdown.change(
            update_chain_dropdown,
            inputs=[file_input, model_dropdown],
            outputs=[chain_dropdown],
        )
        chain_dropdown.change(
            update_residue_dropdown,
            inputs=[file_input, model_dropdown, chain_dropdown],
            outputs=[residue_dropdown],
        )
        residue_dropdown.change(
            process_cif,
            inputs=[
                file_input,
                model_dropdown,
                chain_dropdown,
                residue_dropdown,
            ],
            outputs=[file_output, model_preview],
        )

    return demo
