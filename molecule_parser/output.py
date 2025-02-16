"""
This module contains functions to save the scaled dataset to a CSV file and a glTF file.
"""

import csv

import matplotlib.pyplot
import trimesh


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

    # Get unique radii
    radii = set([atom[3] for atom in atom_coordiantes])

    # Create a color map based on the number of unique radii
    colors = matplotlib.pyplot.cm.get_cmap("jet", len(radii))

    # For each atom, create a sphere and translate it to the correct position.
    for i, atom in enumerate(atom_coordiantes):
        # Get info for the atom
        x, y, z, radius = atom

        # Create a sphere using an icosphere
        sphere = trimesh.creation.icosphere(subdivisions=0, radius=radius)

        # Move the sphere to the atom's location
        sphere.apply_translation([x, y, z])

        # Assign color based on radius
        color = colors(list(radii).index(radius))

        sphere.visual.face_colors = color

        # Add the sphere to the scene with a unique node name and color
        scene.add_geometry(sphere, node_name=f"sphere_{i}")

    # Export the entire scene as a glTF file.
    mesh = scene.to_geometry()
    mesh.export(model_name)
