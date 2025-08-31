import numpy as np
import random
from .bond_calculator import BondCalculator
from .periodicity_calculator import PeriodicityCalculator

class SupercellVerifier:
    """Class that performs supercell structure verification"""
    
    def __init__(self):
        """Initialize SupercellVerifier"""
        self.bond_calc = BondCalculator()
    
    def verify_supercell_structure(self, original_atoms, supercell_atoms, cell_vectors, supercell_vectors):
        """
        Verifies the validity of the supercell structure.
        
        Args:
            original_atoms (list): Original atom information
            supercell_atoms (list): Supercell atom information  
            cell_vectors (dict): Original cell vectors
            supercell_vectors (dict): Supercell vectors
            
        Returns:
            bool: Verification pass status
        """
        print("\n" + "="*80)
        print("Supercell structure verification")
        print("="*80)
        
        # Step 1: Analyze bonding information of original structure
        print("\n[ Step 1: Analyze bonding information of original structure ]")
        original_bond_info = self._analyze_original_structure(original_atoms, cell_vectors)
        
        # Step 2: Select random atoms from supercell
        print("\n[ Step 2: Select random atoms from supercell ]")
        selected_atoms = self._select_random_atoms(supercell_atoms, 100)
        
        # Step 3: Strict analysis of bonding information of selected atoms
        print("\n[ Step 3: Strict analysis of bonding information of selected atoms ]")
        try:
            supercell_bond_info = self._analyze_selected_atoms_strictly(
                selected_atoms, supercell_atoms, supercell_vectors
            )
        except Exception as e:
            print(f"❌ Error found in strict bonding analysis: {str(e)}")
            return False
        
        # Step 4: Compare and verify bonding information
        print("\n[ Step 4: Compare and verify bonding information ]")
        verification_result = self._compare_bond_patterns(original_bond_info, supercell_bond_info)
        
        return verification_result
    
    def _analyze_original_structure(self, original_atoms, cell_vectors):
        """
        Analyzes bonding information of the original structure.
        
        Args:
            original_atoms (list): Original atom information
            cell_vectors (dict): Cell vector information
            
        Returns:
            dict: Bonding information of original structure
        """
        print(f"Original structure: {len(original_atoms)} atoms")
        
        # Calculate bonding information for all atoms
        original_bonds = self.bond_calc.find_bonds(original_atoms, cell_vectors, max_bonds_per_atom=12)
        
        # Store bonding pattern of each atom
        bond_patterns = {}
        for i, atom in enumerate(original_atoms):
            element = atom['element']
            bonds = original_bonds.get(i, [])
            
            # Organize bonding pattern (count and distance by partner element)
            pattern = self._extract_bond_pattern(bonds)
            
            if element not in bond_patterns:
                bond_patterns[element] = []
            bond_patterns[element].append(pattern)
            
            print(f"  Atom {i+1} ({element}): {len(bonds)} bonds")
            for j, bond in enumerate(bonds[:5]):  # Print only first 5
                bond_vector = bond['bond_vector']
                print(f"    Bond {j+1}: {bond['partner_element']}, distance {bond['distance']:.3f}Å, "
                      f"vector ({bond_vector[0]:.3f}, {bond_vector[1]:.3f}, {bond_vector[2]:.3f})")
            if len(bonds) > 5:
                print(f"    ... and {len(bonds)-5} more")
        
        return bond_patterns
    
    def _select_random_atoms(self, supercell_atoms, num_atoms):
        """
        Selects atoms randomly from the supercell.
        
        Args:
            supercell_atoms (list): Supercell atom information
            num_atoms (int): Number of atoms to select
            
        Returns:
            list: Selected atoms and indices
        """
        # Random selection
        random.seed(42)  # Reproducible results
        selected_count = min(num_atoms, len(supercell_atoms))
        selected_indices = random.sample(range(len(supercell_atoms)), selected_count)
        
        selected_atoms_with_idx = []
        for i, idx in enumerate(selected_indices):
            atom = supercell_atoms[idx]
            selected_atoms_with_idx.append((idx, atom))
        
        print(f"Randomly selected {len(selected_atoms_with_idx)} atoms out of {len(supercell_atoms)} total atoms")
        print("Selected atoms:")
        for i, (idx, atom) in enumerate(selected_atoms_with_idx[:10]):  # Print only first 10
            print(f"  {i+1}. {atom['element']}{idx+1}: "
                  f"coordinates ({atom['cart_x']:.3f}, {atom['cart_y']:.3f}, {atom['cart_z']:.3f})")
        if len(selected_atoms_with_idx) > 10:
            print(f"  ... and {len(selected_atoms_with_idx)-10} more")
        
        return selected_atoms_with_idx
    
    def _calculate_supercell_center(self, supercell_atoms):
        """
        Calculates the center point of the supercell.
        
        Args:
            supercell_atoms (list): Supercell atom information
            
        Returns:
            np.array: Center point coordinates
        """
        coordinates = []
        for atom in supercell_atoms:
            coord = np.array([atom['cart_x'], atom['cart_y'], atom['cart_z']])
            coordinates.append(coord)
        
        center = np.mean(coordinates, axis=0)
        return center
    
    def _analyze_selected_atoms_strictly(self, selected_atoms_with_idx, all_supercell_atoms, supercell_vectors):
        """
        Strictly analyzes bonding information of selected atoms.
        
        Args:
            selected_atoms_with_idx (list): Selected atoms and indices
            all_supercell_atoms (list): All supercell atoms
            supercell_vectors (dict): Supercell vectors
            
        Returns:
            dict: Bonding information analysis results
        """
        print(f"Strictly analyzing bonding information of {len(selected_atoms_with_idx)} selected atoms...")
        
        # Extract tolerance range from Step1 bonding information
        bond_constraints = self._extract_bond_constraints()
        
        # Calculate bonding information for each selected atom
        bond_patterns = {}
        unusual_bonds_found = []
        
        for i, (atom_idx, atom) in enumerate(selected_atoms_with_idx[:10]):  # Analyze only first 10
            element = atom['element']
            
            # Calculate bonding information for this atom (considering periodicity)
            atom_bonds = self._find_bonds_for_atom_with_periodicity(
                atom, all_supercell_atoms, supercell_vectors
            )
            
            # Strict validation for each bond
            for bond in atom_bonds:
                is_valid, error_msg = self._validate_bond_strictly(bond, bond_constraints)
                if not is_valid:
                    unusual_bonds_found.append({
                        'atom_idx': atom_idx,
                        'atom_element': element,
                        'bond': bond,
                        'error': error_msg
                    })
            
            # Extract bonding pattern
            pattern = self._extract_bond_pattern(atom_bonds)
            
            if element not in bond_patterns:
                bond_patterns[element] = []
            bond_patterns[element].append(pattern)
            
            print(f"  Selected atom {i+1} ({element}{atom_idx+1}): {len(atom_bonds)} bonds")
            for j, bond in enumerate(atom_bonds[:5]):  # Print only first 5
                bond_vector = bond['bond_vector']
                print(f"    Bond {j+1}: {bond['partner_element']}, distance {bond['distance']:.3f}Å, "
                      f"vector ({bond_vector[0]:.3f}, {bond_vector[1]:.3f}, {bond_vector[2]:.3f})")
            if len(atom_bonds) > 5:
                print(f"    ... and {len(atom_bonds)-5} more")
        
        # Raise error if abnormal bonds are found
        if unusual_bonds_found:
            self._report_unusual_bonds(unusual_bonds_found)
            raise Exception(f"Total {len(unusual_bonds_found)} abnormal bonds were found.")
        
        print("✓ All bonds are within Step1 standard range.")
        return bond_patterns
    
    def _find_bonds_for_atom(self, target_atom, all_atoms):
        """
        Finds bonding information for a specific atom.
        
        Args:
            target_atom (dict): Target atom
            all_atoms (list): All atom list
            
        Returns:
            list: Bonding information list
        """
        target_coord = np.array([target_atom['cart_x'], target_atom['cart_y'], target_atom['cart_z']])
        bonds = []
        
        for other_atom in all_atoms:
            other_coord = np.array([other_atom['cart_x'], other_atom['cart_y'], other_atom['cart_z']])
            
            # Exclude self
            if np.allclose(target_coord, other_coord, atol=1e-6):
                continue
            
            # Calculate distance
            distance = np.linalg.norm(target_coord - other_coord)
            
            # Check bonding threshold
            threshold = self.bond_calc.get_bond_threshold(target_atom['element'], other_atom['element'])
            
            if distance <= threshold:
                bond_vector = other_coord - target_coord
                
                bonds.append({
                    'partner_element': other_atom['element'],
                    'distance': distance,
                    'bond_vector': bond_vector.tolist()
                })
        
        # Sort by distance and return maximum 12
        bonds.sort(key=lambda x: x['distance'])
        return bonds[:12]
    
    def _extract_bond_pattern(self, bonds):
        """
        Extracts pattern from bonding list.
        
        Args:
            bonds (list): Bonding information list
            
        Returns:
            dict: Bonding pattern (count and average distance by element)
        """
        pattern = {}
        
        for bond in bonds:
            partner = bond['partner_element']
            distance = bond['distance']
            
            if partner not in pattern:
                pattern[partner] = {'count': 0, 'distances': []}
            
            pattern[partner]['count'] += 1
            pattern[partner]['distances'].append(distance)
        
        # Calculate average distance
        for partner in pattern:
            distances = pattern[partner]['distances']
            pattern[partner]['avg_distance'] = np.mean(distances)
            pattern[partner]['distance_range'] = (min(distances), max(distances))
        
        return pattern
    
    def _compare_bond_patterns(self, original_patterns, supercell_patterns):
        """
        Compares bonding patterns between original and supercell.
        
        Args:
            original_patterns (dict): Bonding pattern of original structure
            supercell_patterns (dict): Bonding pattern of supercell
            
        Returns:
            bool: Verification success status
        """
        print("Bonding pattern comparison analysis:")
        
        success_count = 0
        total_comparisons = 0
        tolerance = 0.1  # Distance tolerance (Å)
        
        for element in original_patterns:
            if element not in supercell_patterns:
                print(f"  Warning: {element} element not found in supercell center")
                continue
            
            original_list = original_patterns[element]
            supercell_list = supercell_patterns[element]
            
            print(f"\n  {element} element comparison:")
            print(f"    Original: {len(original_list)} atoms")
            print(f"    Supercell: {len(supercell_list)} atoms")
            
            # Compare bonding pattern of each atom
            for i, original_pattern in enumerate(original_list):
                total_comparisons += 1
                best_match = self._find_best_pattern_match(original_pattern, supercell_list, tolerance)
                
                if best_match:
                    success_count += 1
                    print(f"    Original atom {i+1}: Matching successful")
                else:
                    print(f"    Original atom {i+1}: Matching failed")
                    self._print_pattern_details(original_pattern, "Original")
        
        # Verification results
        success_rate = (success_count / total_comparisons * 100) if total_comparisons > 0 else 0
        
        print(f"\n[ Verification Results ]")
        print(f"Total comparisons: {total_comparisons}")
        print(f"Matching successful: {success_count}")
        print(f"Success rate: {success_rate:.1f}%")
        
        # Consider success if 80% or more
        is_success = success_rate >= 80.0
        
        if is_success:
            print("✓ Supercell structure verification successful: Original structure was correctly replicated.")
        else:
            print("✗ Supercell structure verification failed: There may be issues with structure replication.")
        
        return is_success
    
    def _find_best_pattern_match(self, target_pattern, candidate_patterns, tolerance):
        """
        Finds the best matching pattern.
        
        Args:
            target_pattern (dict): Target pattern
            candidate_patterns (list): Candidate patterns
            tolerance (float): Tolerance
            
        Returns:
            dict or None: Matched pattern or None
        """
        for candidate in candidate_patterns:
            if self._patterns_match(target_pattern, candidate, tolerance):
                return candidate
        return None
    
    def _patterns_match(self, pattern1, pattern2, tolerance):
        """
        Checks if two patterns match.
        
        Args:
            pattern1, pattern2 (dict): Patterns to compare
            tolerance (float): Tolerance
            
        Returns:
            bool: Matching status
        """
        # Check if bonded element types are the same
        if set(pattern1.keys()) != set(pattern2.keys()):
            return False
        
        # Check if bonding count and distance are similar for each element
        for partner in pattern1:
            p1 = pattern1[partner]
            p2 = pattern2[partner]
            
            # Check bonding count
            if p1['count'] != p2['count']:
                return False
            
            # Check average distance
            distance_diff = abs(p1['avg_distance'] - p2['avg_distance'])
            if distance_diff > tolerance:
                return False
        
        return True
    
    def _print_pattern_details(self, pattern, label):
        """
        Prints pattern details.
        
        Args:
            pattern (dict): Pattern to print
            label (str): Pattern label
        """
        print(f"      {label} pattern:")
        for partner, info in pattern.items():
            print(f"        {partner}: {info['count']} bonds, average distance {info['avg_distance']:.3f}Å") 
    
    def _extract_bond_constraints(self):
        """
        Extracts bonding constraints from Step1 bonding information.
        
        Returns:
            dict: Bonding constraints
        """
        # Set default values
        return {
            'default_min_distance': 0.5,
            'default_max_distance': 15.0
        }
    
    def _find_bonds_for_atom_with_periodicity(self, target_atom, all_atoms, supercell_vectors):
        """
        Finds bonding information for a specific atom considering periodicity.
        
        Args:
            target_atom (dict): Target atom
            all_atoms (list): All atom list
            supercell_vectors (dict): Supercell vectors
            
        Returns:
            list: Bonding information list
        """
        target_coord = np.array([target_atom['cart_x'], target_atom['cart_y'], target_atom['cart_z']])
        bonds = []
        
        # Supercell vectors
        a_vec = np.array(supercell_vectors['a'])
        b_vec = np.array(supercell_vectors['b'])
        c_vec = np.array(supercell_vectors['c'])
        
        for other_atom in all_atoms:
            other_coord = np.array([other_atom['cart_x'], other_atom['cart_y'], other_atom['cart_z']])
            
            # Exclude self
            if np.allclose(target_coord, other_coord, atol=1e-6):
                continue
            
            # Calculate distance considering periodicity (check 27 images)
            min_distance = float('inf')
            best_vector = None
            
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        # Coordinates of periodic image
                        periodic_coord = other_coord + i*a_vec + j*b_vec + k*c_vec
                        bond_vector = periodic_coord - target_coord
                        distance = np.linalg.norm(bond_vector)
                        
                        if distance < min_distance:
                            min_distance = distance
                            best_vector = bond_vector
            
            # Check bonding threshold
            threshold = self.bond_calc.get_bond_threshold(target_atom['element'], other_atom['element'])
            
            if min_distance <= threshold:
                bonds.append({
                    'partner_element': other_atom['element'],
                    'distance': min_distance,
                    'bond_vector': best_vector.tolist()
                })
        
        # Sort by distance and return maximum 12
        bonds.sort(key=lambda x: x['distance'])
        return bonds[:12]
    
    def _validate_bond_strictly(self, bond, bond_constraints):
        """
        Strictly validates bonding.
        
        Args:
            bond (dict): Bond to validate
            bond_constraints (dict): Bonding constraints
            
        Returns:
            tuple: (Validation success status, error message)
        """
        distance = bond['distance']
        partner_element = bond['partner_element']
        
        # Find constraints for this bond type
        constraints = None
        for bond_type, constraint in bond_constraints.items():
            if bond_type == 'default':
                continue
            if partner_element in bond_type or bond_type in partner_element:
                constraints = constraint
                break
        
        # Use default constraints
        if constraints is None:
            constraints = bond_constraints.get('default', {
                'min_distance': 0.5,
                'max_distance': 15.0
            })
        
        min_dist = constraints['min_distance']
        max_dist = constraints['max_distance']
        
        # Distance validation
        if distance < min_dist:
            return False, f"Bond distance {distance:.3f}Å is less than minimum value {min_dist:.3f}Å"
        
        if distance > max_dist:
            return False, f"Bond distance {distance:.3f}Å is greater than maximum value {max_dist:.3f}Å"
        
        return True, ""
    
    def _report_unusual_bonds(self, unusual_bonds):
        """
        Reports abnormal bonds.
        
        Args:
            unusual_bonds (list): Abnormal bond list
        """
        print(f"\n❌ {len(unusual_bonds)} abnormal bonds were found:")
        
        for i, unusual in enumerate(unusual_bonds[:5]):  # Print only first 5
            atom_idx = unusual['atom_idx']
            atom_element = unusual['atom_element']
            bond = unusual['bond']
            error = unusual['error']
            
            print(f"  {i+1}. Atom {atom_element}{atom_idx+1} → {bond['partner_element']}")
            print(f"     Distance: {bond['distance']:.3f}Å")
            print(f"     Error: {error}")
        
        if len(unusual_bonds) > 5:
            print(f"     ... and {len(unusual_bonds)-5} more")
        
        print("\n⚠️ The supercell may not have been generated correctly.") 