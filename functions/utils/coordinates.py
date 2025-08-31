import numpy as np

def convert_to_cartesian(frac_coords, cell_vectors):
    """
    Convert fractional coordinates to cartesian coordinates (angstrom).
    
    Args:
        frac_coords (list): [x, y, z] fractional coordinates
        cell_vectors (dict): Dictionary containing cell vectors a, b, c
    
    Returns:
        list: [x, y, z] cartesian coordinates in angstrom
    """
    x = (frac_coords[0] * cell_vectors['a'][0] + 
         frac_coords[1] * cell_vectors['b'][0] + 
         frac_coords[2] * cell_vectors['c'][0])
    
    y = (frac_coords[0] * cell_vectors['a'][1] + 
         frac_coords[1] * cell_vectors['b'][1] + 
         frac_coords[2] * cell_vectors['c'][1])
    
    z = (frac_coords[0] * cell_vectors['a'][2] + 
         frac_coords[1] * cell_vectors['b'][2] + 
         frac_coords[2] * cell_vectors['c'][2])
    
    return [x, y, z]

def convert_atoms_to_cartesian(atoms, cell_vectors):
    """
    Convert atomic positions from fractional to cartesian coordinates.
    
    Args:
        atoms (list): List of atomic information dictionaries
        cell_vectors (dict): Dictionary containing cell vectors
    
    Returns:
        list: List of atoms with cartesian coordinates added
    """
    if not atoms or cell_vectors is None:
        return atoms
    
    cartesian_atoms = []
    for atom in atoms:
        # Extract fractional coordinates using consistent key names
        frac_x = float(atom.get('x', 0.0))
        frac_y = float(atom.get('y', 0.0))
        frac_z = float(atom.get('z', 0.0))
        frac_coords = [frac_x, frac_y, frac_z]
        
        cart_coords = convert_to_cartesian(frac_coords, cell_vectors)
        
        # Create new atom dictionary with both fractional and cartesian coordinates
        new_atom = atom.copy()  # Preserve existing information
        new_atom.update({
            'x': frac_x,
            'y': frac_y,
            'z': frac_z,
            'cart_x': cart_coords[0],
            'cart_y': cart_coords[1],
            'cart_z': cart_coords[2]
        })
        cartesian_atoms.append(new_atom)
    
    return cartesian_atoms

def calculate_distance(pos1, pos2):
    """
    Calculate distance between two points in 3D space.
    
    Args:
        pos1 (list/array): [x, y, z] coordinates of first point
        pos2 (list/array): [x, y, z] coordinates of second point
    
    Returns:
        float: Distance between the points in Angstroms
    """
    return np.sqrt(np.sum((np.array(pos1) - np.array(pos2)) ** 2))

def calculate_vector(pos1, pos2):
    """
    Calculate vector from pos1 to pos2.
    
    Args:
        pos1 (list/array): [x, y, z] coordinates of first point
        pos2 (list/array): [x, y, z] coordinates of second point
    
    Returns:
        array: Vector from pos1 to pos2
    """
    return np.array(pos2) - np.array(pos1)

def normalize_vector(vector):
    """
    Normalize a vector.
    
    Args:
        vector (array): Vector to normalize
    
    Returns:
        array: Normalized vector
    """
    norm = np.linalg.norm(vector)
    return vector / (norm + 1e-10)  # Add small value for numerical stability

def calculate_angle(v1, v2):
    """
    Calculate angle between two vectors.
    
    Args:
        v1 (array): First vector
        v2 (array): Second vector
    
    Returns:
        float: Angle in degrees
    """
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    
    # Handle case when vector is too short
    if norm1 < 1e-10 or norm2 < 1e-10:
        return 0.0
    
    cos_angle = np.dot(v1, v2) / (norm1 * norm2)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical errors
    
    return np.degrees(np.arccos(cos_angle)) 