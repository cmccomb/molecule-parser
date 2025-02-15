"""
This module contains the functions to launch the molecule parser GUI.
"""

import os

import gradio
from Bio.PDB import PDBList

from .parsing import (
    process_cif,
    load_cif_file,
    extract_chains_from_model,
    extract_residues_from_chain,
)


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
        pdb_filename = pdb_list.retrieve_pdb_file(code, file_format="mmCif")

        # Check if file at pdb_filename exists us os.path.exists
        if os.path.exists(pdb_filename) is False:
            pdb_filename = None
            models = []
            gradio.Error(f"{code} is not a valid PDB code.")
        else:
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


def launch(**kwargs) -> tuple[gradio.routes.App, str, str]:
    """
    Launch the molecule parser GUI.

    Args:
        **kwargs: Any keyword argument that can be passed to the gradio.App launch method.

    Returns:
        demo (gradio.App): The Gradio App object.

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

    return demo.launch(**kwargs)
