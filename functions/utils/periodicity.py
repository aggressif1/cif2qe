import numpy as np
from itertools import product

def generate_periodic_images(atoms, cell_vectors, radius=1):
    """
    Generate periodic images of atoms in neighboring unit cells.
    
    Args:
        atoms (list): List of atomic information dictionaries
        cell_vectors (dict): Dictionary containing cell vectors a, b, c
        radius (int): Number of unit cells to consider in each direction (+/-)
    
    Returns:
        list: List of dictionaries containing atomic information for all atoms 
              (original + periodic images)
    """
    # Convert cell vectors to numpy arrays for easier calculation
    a_vec = np.array(cell_vectors['a'])
    b_vec = np.array(cell_vectors['b'])
    c_vec = np.array(cell_vectors['c'])
    
    # Generate all possible combinations of translations
    translations = list(product(range(-radius, radius+1), repeat=3))
    
    # List to store all atoms including periodic images
    all_atoms = []
    
    # For each atom in the original unit cell
    for atom in atoms:
        pos = np.array([atom['cart_x'], atom['cart_y'], atom['cart_z']])
        
        # For each possible translation
        for tx, ty, tz in translations:
            # Calculate the translated position
            translation = tx * a_vec + ty * b_vec + tz * c_vec
            new_pos = pos + translation
            
            # Create new atom dictionary with translated position
            new_atom = atom.copy()
            new_atom.update({
                'cart_x': float(new_pos[0]),
                'cart_y': float(new_pos[1]),
                'cart_z': float(new_pos[2]),
                'cell_x': tx,
                'cell_y': ty,
                'cell_z': tz,
                'is_periodic_image': not (tx == ty == tz == 0)
            })
            
            all_atoms.append(new_atom)
    
    return all_atoms

def is_periodic_bond(atom1, atom2):
    """
    Check if the bond between two atoms crosses periodic boundaries.
    
    Args:
        atom1 (dict): First atom information
        atom2 (dict): Second atom information
    
    Returns:
        bool: True if the bond crosses periodic boundaries
    """
    return (atom1.get('is_periodic_image', False) != 
            atom2.get('is_periodic_image', False))

def get_periodic_bond_description(atom1, atom2):
    """
    Generate a description of the periodic bond between two atoms.
    
    Args:
        atom1 (dict): First atom information
        atom2 (dict): Second atom information
    
    Returns:
        str: Description of the periodic bond
    """
    if not is_periodic_bond(atom1, atom2):
        return "same cell"
    
    # Determine which atom is in the original cell
    orig_atom = atom1 if not atom1.get('is_periodic_image', False) else atom2
    per_atom = atom2 if orig_atom is atom1 else atom1
    
    # Get the cell indices of the periodic image
    cell_x = per_atom.get('cell_x', 0)
    cell_y = per_atom.get('cell_y', 0)
    cell_z = per_atom.get('cell_z', 0)
    
    # Create direction description
    directions = []
    if cell_x != 0:
        directions.append("+" + "a" if cell_x > 0 else "-" + "a")
    if cell_y != 0:
        directions.append("+" + "b" if cell_y > 0 else "-" + "b")
    if cell_z != 0:
        directions.append("+" + "c" if cell_z > 0 else "-" + "c")
    
    return "across " + ",".join(directions)
