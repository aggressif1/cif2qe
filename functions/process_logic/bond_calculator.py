import numpy as np
from .periodicity_calculator import PeriodicityCalculator

class BondCalculator:
    """Class for calculating bond information between atoms"""
    
    def __init__(self):
        # Default bond distance thresholds for each element (Å) - approximately twice the covalent radius of each element
        self.bond_thresholds = {
            'H': 1.5,   # Hydrogen
            'C': 1.8,   # Carbon
            'N': 1.8,   # Nitrogen
            'O': 1.8,   # Oxygen
            'Si': 2.5,  # Silicon
            'Ge': 2.6,  # Germanium
            'Sn': 3.0,  # Tin
            'Pb': 3.2,  # Lead
            'default': 3.0  # Default value
        }
    
    def calculate_distance(self, atom1, atom2):
        """
        Calculate distance between two atoms.
        
        Args:
            atom1, atom2 (dict): Atomic information {'cart_x': float, 'cart_y': float, 'cart_z': float}
        
        Returns:
            float: Distance between two atoms (Å)
        """
        dx = atom1['cart_x'] - atom2['cart_x']
        dy = atom1['cart_y'] - atom2['cart_y']
        dz = atom1['cart_z'] - atom2['cart_z']
        
        return np.sqrt(dx*dx + dy*dy + dz*dz)
    
    def calculate_bond_vector(self, atom1, atom2):
        """
        Calculate bond vector from atom1 to atom2.
        
        Args:
            atom1, atom2 (dict): Atomic information
        
        Returns:
            tuple: (dx, dy, dz) bond vector
        """
        dx = atom2['cart_x'] - atom1['cart_x']
        dy = atom2['cart_y'] - atom1['cart_y']
        dz = atom2['cart_z'] - atom1['cart_z']
        
        return (dx, dy, dz)
    
    def get_bond_threshold(self, element1, element2):
        """
        Return bond distance threshold between two elements.
        
        Args:
            element1, element2 (str): Element symbols
        
        Returns:
            float: Bond distance threshold (Å)
        """
        threshold1 = self.bond_thresholds.get(element1, self.bond_thresholds['default'])
        threshold2 = self.bond_thresholds.get(element2, self.bond_thresholds['default'])
        
        # Use average threshold of two elements
        return (threshold1 + threshold2) / 2
    
    def find_bonds(self, atoms, cell_vectors=None, max_bonds_per_atom=12, min_bond_distance=None, max_bond_distance=None, target_indices=None):
        """
        Find all bonds in atom list. Considers periodic boundary conditions.
        
        Args:
            atoms (list): List of atomic information
            cell_vectors (dict): Cell vector information {'a': [x,y,z], 'b': [x,y,z], 'c': [x,y,z]}
            max_bonds_per_atom (int): Maximum number of bonds per atom
            min_bond_distance (float): Minimum bond distance (Å)
            max_bond_distance (float): Maximum bond distance (Å)
            target_indices (list): List of specific atom indices to analyze (if None, analyze all atoms)
        
        Returns:
            dict: Bond information for each atom {atom_index: [bond_info, ...]}
        """
        bonds = {}
        
        # Initialize periodicity calculator
        periodicity_calc = PeriodicityCalculator(cell_vectors)
        
        # Generate atoms from neighboring unit cells to consider periodicity
        extended_atoms = periodicity_calc.create_periodic_images(atoms, 1, 1, 1) if cell_vectors else atoms
        
        # Determine atom indices to analyze
        if target_indices is not None:
            target_atom_indices = target_indices
        else:
            target_atom_indices = range(len(atoms))
        
        for i in target_atom_indices:
            if i >= len(atoms):  # Check index range
                continue
            
            atom1 = atoms[i]
            bonds[i] = []
            distances = []
            
            # Calculate distances with extended atom list
            for atom2 in extended_atoms:
                if not (atom2.get('is_original', False) and atom2.get('original_index', -1) == i):  # Don't bond with itself
                    distance = self.calculate_distance(atom1, atom2)
                    
                    # Check bond distance range
                    is_valid_bond = False
                    threshold = None  # Set default value
                    
                    if min_bond_distance is not None and max_bond_distance is not None:
                        # Use user-specified range
                        is_valid_bond = min_bond_distance <= distance <= max_bond_distance
                        threshold = max_bond_distance  # Record maximum value as threshold when range is set
                    else:
                        # Use default threshold
                        threshold = self.get_bond_threshold(atom1['element'], atom2['element'])
                        is_valid_bond = distance <= threshold
                    
                    if is_valid_bond:
                        bond_vector = self.calculate_bond_vector(atom1, atom2)
                        
                        # Original unit cell atom index
                        original_index = atom2.get('original_index', 0)
                        
                        distances.append({
                            'partner_index': original_index,
                            'partner_element': atom2['element'],
                            'distance': distance,
                            'bond_vector': bond_vector,
                            'threshold': threshold,
                            'is_periodic': not atom2.get('is_original', False)  # Indicate if it's a periodic image
                        })
            
            # Remove duplicates ensuring uniqueness of bond vectors
            unique_bonds = self._remove_duplicate_bonds(distances)
            
            # Sort by distance and select only maximum number of bonds
            unique_bonds.sort(key=lambda x: x['distance'])
            bonds[i] = unique_bonds[:max_bonds_per_atom]
        
        return bonds
    
    def _remove_duplicate_bonds(self, bond_list):
        """
        Remove duplicate bonds based on bond vector similarity.
        
        Args:
            bond_list (list): List of bond information
            
        Returns:
            list: List of bond information with duplicates removed
        """
        if not bond_list:
            return []
        
        unique_bonds = []
        tolerance = 0.01  # Tolerance for bond vector comparison
        
        for bond in bond_list:
            vector = bond['bond_vector']
            is_duplicate = False
            
            # Compare with existing unique bonds
            for unique_bond in unique_bonds:
                unique_vector = unique_bond['bond_vector']
                
                # Check if each component of the vector matches within tolerance
                if (abs(vector[0] - unique_vector[0]) < tolerance and
                    abs(vector[1] - unique_vector[1]) < tolerance and
                    abs(vector[2] - unique_vector[2]) < tolerance):
                    is_duplicate = True
                    break
            
            # Add if not duplicate
            if not is_duplicate:
                unique_bonds.append(bond)
        
        return unique_bonds
    
    def format_bond_info(self, bonds, atom_index):
        """
        Return formatted string of bond information for specific atom.
        
        Args:
            bonds (dict): Bond information
            atom_index (int): Atom index
        
        Returns:
            list: List of formatted bond information strings
        """
        if atom_index not in bonds or not bonds[atom_index]:
            return ["  No bonds"]
        
        bond_info = []
        for bond in bonds[atom_index]:
            partner_idx = bond['partner_index']
            partner_elem = bond['partner_element']
            distance = bond['distance']
            bv = bond['bond_vector']
            
            info_line = f"  → {partner_elem}{partner_idx+1}: {distance:.3f}Å, vector({bv[0]:6.3f}, {bv[1]:6.3f}, {bv[2]:6.3f})"
            bond_info.append(info_line)
        
        return bond_info 