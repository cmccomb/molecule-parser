import csv
import os

import gradio
import numpy
import pyvista
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
import Bio.PDB


def load_cif_file(file_name: str):
    # Create an MMCIFParser Object
    parser = MMCIFParser()

    # Extract the overall structure
    structure = parser.get_structure(
        os.path.splitext(os.path.basename(file_name))[0], file_name
    )  # Using inputs

    # read in the CIF file.
    pdb_info = MMCIF2Dict(file_name)

    # structureID = pdb_info["_pdbx_database_status.entry_id"][0]
    structureID = pdb_info["_entry.id"][0]

    # Extract models
    models = sorted([model.id for model in structure])

    # Return stuff
    return structure, structureID, models


def extract_chains_from_model(structure, selected_model):
    for model in structure:
        if model.id == selected_model:
            chains = sorted([chain.id for chain in model])
            return chains
    return []


def extract_residues_from_chain(structure, selected_model, selected_chain):
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


def extract_chain_atoms(structure, selected_model, selected_chain, selected_residues):
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


def make_data_the_proper_format(atom_coords):

    scaleDown = 0.68

    # Do the operation for the desired residue, entity, model, etc...
    newDataSet = atom_coords
    for atom in newDataSet:
        if atom[0] == "N":
            atom.append(round(1.55 * scaleDown, 1))  # add radius
            atom.pop(0)  # delete first item in the list (atom symbol)
        elif atom[0] == "C":
            atom.append(round(1.7 * scaleDown, 1))
            atom.pop(0)
        elif atom[0] == "O":
            atom.append(round(1.52 * scaleDown, 1))
            atom.pop(0)
        elif atom[0] == "H":
            atom.append(round(1.10 * scaleDown, 2))
            atom.pop(0)
        elif atom[0] == "P":
            atom.append(round(1.8 * scaleDown, 1))
            atom.pop(0)
        else:  # By default all other atoms not listed above are given a radius equal to a carbon atom
            atom.append(round(1.55 * scaleDown, 1))
            atom.pop(0)

    # Scale the values.
    scaledDataSet = [[round(entry * 10, 3) for entry in atom] for atom in newDataSet]

    return scaledDataSet


def save_csv(csvName, scaledDataSet):

    # Non-poly residues don't use neighbor search to export extra info... so just use this...

    header = ["X_Coord", "Y_Coord", "Z_Coord", "Atomic_Radius"]

    with open(csvName, "w") as f:
        write = csv.writer(f)

        write.writerow(header)
        write.writerows(scaledDataSet)


def process_cif(file, selected_model, selected_chain, selected_residues):
    if file is None:
        return gradio.File(
            None, label="Download CSV", interactive=False
        ), gradio.Model3D(None, label="Model Preview", interactive=False)
    else:
        structure, structureID, _ = load_cif_file(file.name)
        atom_coords = extract_chain_atoms(
            structure, selected_model, selected_chain, selected_residues
        )
        if not atom_coords:
            return "Error: No atoms found for the selected criteria."
        scaled_dataset = make_data_the_proper_format(atom_coords)
        csvName = (
            f"{structureID}_Model_{selected_model}_Chain_{selected_chain}_Residues.csv"
        )
        save_csv(csvName, scaled_dataset)

        # Create a pyvista plotter
        pyvista.start_xvfb()
        plotter = pyvista.Plotter()

        # Create a point cloud from your scaled dataset
        scaledData = numpy.array(scaled_dataset)
        points = scaledData[:, 0:3]
        radii = scaledData[:, 3]

        # Create a PolyData object with the points
        point_cloud = pyvista.PolyData(points)
        # Add the radii as point data
        point_cloud["radius"] = radii

        # Use glyphing to place a sphere at each point scaled by its radius
        glyph = point_cloud.glyph(
            scale="radius",
            geom=pyvista.Sphere(radius=1.0, phi_resolution=5, theta_resolution=5),
        )
        plotter.add_mesh(glyph)

        # Save the plot to a file
        plotter.export_gltf("model.gltf", inline_data=True)
        return gradio.File(
            csvName, label=f"Download {csvName}", interactive=True
        ), gradio.Model3D("model.gltf", label="Model Preview")


def update_model_dropdown_from_file(file):
    if file is None:
        models = []
    else:
        _, _, models = load_cif_file(file.name)
    return gradio.update(choices=models, value=models[0] if models else None)


def update_model_dropdown_from_code(code: str):
    if len(code) < 4:
        pdb_filename = None
        models = []
    else:

        # Create an instance of the PDBList class
        pdb_list = Bio.PDB.PDBList()

        # Download the MMCIF file using the retrieve_pdb_file method
        pdb_filename = pdb_list.retrieve_pdb_file(
            code, pdir="data", file_format="mmCif"
        )

        _, _, models = load_cif_file(pdb_filename)

    return gradio.update(
        choices=models, value=models[0] if models else None
    ), gradio.update(value=pdb_filename)


def update_chain_dropdown(file, selected_model):
    if file is None:
        chains = []
    else:
        structure, _, _ = load_cif_file(file.name)
        chains = extract_chains_from_model(structure, selected_model)
    return gradio.update(choices=chains, value=chains[0] if chains else None)


def update_residue_dropdown(file, selected_model, selected_chain):
    if file is None:
        residues = []
    else:
        structure, _, _ = load_cif_file(file.name)
        residues = extract_residues_from_chain(
            structure, selected_model, selected_chain
        )
    return gradio.update(choices=residues, value=residues if residues else [])


def gui():

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
