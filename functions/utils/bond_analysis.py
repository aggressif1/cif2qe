import numpy as np
from itertools import combinations
from .periodicity import generate_periodic_images
from .coordinates import (
    calculate_distance,
    calculate_vector,
    normalize_vector,
    calculate_angle
)

def calculate_bond_vector(pos1, pos2, cell_vectors):
    """
    Calculate the bond vector from pos1 to pos2 considering periodic boundary conditions.
    
    Args:
        pos1 (list): Position of the first atom [x, y, z]
        pos2 (list): Position of the second atom [x, y, z]
        cell_vectors (dict): Unit cell vectors
        
    Returns:
        tuple: (bond_vector, distance)
    """
    # Convert positions to numpy arrays for easier calculation
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    
    # Calculate direct vector
    direct_vector = pos2 - pos1
    
    # Consider periodic images
    cell_matrix = np.array([
        cell_vectors['a'],
        cell_vectors['b'],
        cell_vectors['c']
    ])
    
    # Check all neighboring periodic images
    min_distance = float('inf')
    best_vector = None
    
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                translation = np.dot(cell_matrix.T, [i, j, k])
                periodic_vector = direct_vector + translation
                distance = np.linalg.norm(periodic_vector)
                
                if distance < min_distance:
                    min_distance = distance
                    best_vector = periodic_vector
    
    # Normalize the vector
    normalized_vector = best_vector / min_distance
    
    return best_vector, min_distance

def calculate_bond_angle(pos1, pos2, pos3):
    """
    Calculate angle between three atoms (pos1-pos2-pos3).
    
    Args:
        pos1 (list/array): [x, y, z] coordinates of first atom
        pos2 (list/array): [x, y, z] coordinates of central atom
        pos3 (list/array): [x, y, z] coordinates of third atom
    
    Returns:
        float: Angle in degrees
    """
    v1 = calculate_vector(pos2, pos1)
    v2 = calculate_vector(pos2, pos3)
    return calculate_angle(v1, v2)

def find_atom_bonds(central_atom, all_atoms, max_distance=3.0, min_distance=0.1, cell_vectors=None):
    """
    Find all bonds for a specific atom.
    
    Args:
        central_atom (dict): The atom to analyze
        all_atoms (list): List of all atoms including periodic images
        max_distance (float): Maximum distance to consider for bonds
        min_distance (float): Minimum distance to consider for bonds
        cell_vectors (dict): Unit cell vectors for periodic boundary conditions
    
    Returns:
        list: List of dictionaries containing bond information
    """
    bonds = []
    central_pos = [central_atom['cart_x'], central_atom['cart_y'], central_atom['cart_z']]
    
    for other_atom in all_atoms:
        if other_atom is central_atom:
            continue
            
        other_pos = [other_atom['cart_x'], other_atom['cart_y'], other_atom['cart_z']]
        
        if cell_vectors:
            # Calculate vector and distance considering periodic boundary conditions
            bond_vector, distance = calculate_bond_vector(central_pos, other_pos, cell_vectors)
        else:
            # General distance calculation
        distance = calculate_distance(central_pos, other_pos)
            bond_vector = calculate_vector(central_pos, other_pos)
        
        if min_distance <= distance <= max_distance:
            bond_info = {
                'partner_element': other_atom['element'],
                'distance': distance,
                'vector': bond_vector,
                'partner_atom': other_atom,
                'periodicity': other_atom.get('periodicity', 'same cell')
            }
            bonds.append(bond_info)
    
    return bonds

def analyze_atom_environment(atom, all_atoms, max_distance=3.0):
    """
    Analyze the bonding environment of a specific atom.
    
    Args:
        atom (dict): The atom to analyze
        all_atoms (list): List of all atoms including periodic images
        max_distance (float): Maximum distance to consider for bonds
    
    Returns:
        dict: Dictionary containing bonding environment information
    """
    bonds = find_atom_bonds(atom, all_atoms, max_distance)
    
    # Group bonds by element type
    bonds_by_element = {}
    for bond in bonds:
        element = bond['partner_element']
        if element not in bonds_by_element:
            bonds_by_element[element] = []
        bonds_by_element[element].append(bond)
    
    # Calculate bond angles
    bond_angles = []
    if len(bonds) >= 2:
        for i, j in combinations(range(len(bonds)), 2):
            pos1 = [bonds[i]['partner_atom']['cart_x'],
                   bonds[i]['partner_atom']['cart_y'],
                   bonds[i]['partner_atom']['cart_z']]
            pos2 = [atom['cart_x'], atom['cart_y'], atom['cart_z']]
            pos3 = [bonds[j]['partner_atom']['cart_x'],
                   bonds[j]['partner_atom']['cart_y'],
                   bonds[j]['partner_atom']['cart_z']]
            
            angle = calculate_bond_angle(pos1, pos2, pos3)
            angle_info = {
                'elements': (bonds[i]['partner_element'], 
                           atom['element'],
                           bonds[j]['partner_element']),
                'angle': angle,
                'periodicity': (bonds[i]['periodicity'],
                              bonds[j]['periodicity'])
            }
            bond_angles.append(angle_info)
    
    return {
        'element': atom['element'],
        'position': [atom['cart_x'], atom['cart_y'], atom['cart_z']],
        'bonds_by_element': bonds_by_element,
        'bond_angles': bond_angles
    }

def analyze_structure_bonds(atoms, cell_vectors, cutoff=3.0):
    """
    Analyze bonds in the structure considering periodic boundary conditions.
    """
    atom_environments = []
    
    for i, atom1 in enumerate(atoms):
        env = {
            'element': atom1['element'],
            'position': [atom1['cart_x'], atom1['cart_y'], atom1['cart_z']],
            'bonds_by_element': {},
            'bond_angles': []
        }
        
        # Find bonds
        bonds = []
        for j, atom2 in enumerate(atoms):
            if i == j:
                continue
                
            pos1 = [atom1['cart_x'], atom1['cart_y'], atom1['cart_z']]
            pos2 = [atom2['cart_x'], atom2['cart_y'], atom2['cart_z']]
            
            bond_vector, distance = calculate_bond_vector(pos1, pos2, cell_vectors)
            
            if distance <= cutoff:
                bond_info = {
                    'index': j,
                    'element': atom2['element'],
                    'distance': distance,
                    'vector': bond_vector.tolist()
                }
                bonds.append(bond_info)
                
                if atom2['element'] not in env['bonds_by_element']:
                    env['bonds_by_element'][atom2['element']] = []
                env['bonds_by_element'][atom2['element']].append(bond_info)
        
        # Calculate bond angles
        for b1 in range(len(bonds)):
            for b2 in range(b1 + 1, len(bonds)):
                v1 = np.array(bonds[b1]['vector'])
                v2 = np.array(bonds[b2]['vector'])
                
                cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle = np.degrees(np.arccos(cos_angle))
                
                angle_info = {
                    'elements': [bonds[b1]['element'], atom1['element'], bonds[b2]['element']],
                    'angle': angle,
                    'vectors': [bonds[b1]['vector'], bonds[b2]['vector']]
                }
                env['bond_angles'].append(angle_info)
        
        atom_environments.append(env)
    
    return atom_environments

def format_atom_bonds(atom_index, atom_env):
    """
    Format bond information for a single atom.
    """
    output = []
    output.append(f"Atom {atom_index} ({atom_env['element']}):")
    
    if not atom_env['bonds_by_element']:
        output.append("  No bonds found within cutoff distance")
        return "\n".join(output)
    
    total_bonds = sum(len(bonds) for bonds in atom_env['bonds_by_element'].values())
    output.append(f"  Total bonds: {total_bonds}")
    
    for element, bonds in sorted(atom_env['bonds_by_element'].items()):
        output.append(f"  Bonds with {element}:")
        for i, bond in enumerate(sorted(bonds, key=lambda x: x['distance']), 1):
            vector = bond['vector']
            output.append(f"    Bond {i}: Distance = {bond['distance']:.4f} Å  |  Bond Vector = [{vector[0]:7.4f}, {vector[1]:7.4f}, {vector[2]:7.4f}] Å")
    
    if atom_env['bond_angles']:
        output.append("  Bond Angles:")
        for angle_info in sorted(atom_env['bond_angles'], key=lambda x: x['angle']):
            elements = angle_info['elements']
            angle = angle_info['angle']
            vectors = angle_info['vectors']
            output.append(f"    {elements[0]}-{elements[1]}-{elements[2]}: {angle:.2f}°  |  Vector 1: [{vectors[0][0]:7.4f}, {vectors[0][1]:7.4f}, {vectors[0][2]:7.4f}] Å  |  Vector 2: [{vectors[1][0]:7.4f}, {vectors[1][1]:7.4f}, {vectors[1][2]:7.4f}] Å")
    
    return "\n".join(output)

def format_structure_bonds(atom_environments):
    """
    Format bond analysis results for all atoms.
    
    Args:
        atom_environments (list): List of atom environment dictionaries
    
    Returns:
        str: Formatted analysis results
    """
    output = []
    output.append("\nBond Analysis by Atom:")
    output.append("=" * 40)
    
    for i, env in enumerate(atom_environments, 1):
        output.append(format_atom_bonds(i, env))
    
    return "\n".join(output)
