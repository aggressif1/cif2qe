import numpy as np
from .supercell_verifier import SupercellVerifier
from .bond_calculator import BondCalculator
from ..utils.element_counter import ElementCounter

class MillerUnitCellVerifier:
    """Miller index-based unit cell dedicated verification class"""
    
    def __init__(self):
        """Initialize MillerUnitCellVerifier"""
        self.supercell_verifier = SupercellVerifier()
        self.bond_calc = BondCalculator()
    
    def verify_miller_unit_cell(self, original_unit_atoms, miller_unit_atoms, original_vectors, miller_vectors):
        """
        Compare and verify Miller index-based unit cell with original unit cell.
        
        Args:
            original_unit_atoms (list): Original unit cell atoms
            miller_unit_atoms (list): Miller index-based unit cell atoms
            original_vectors (dict): Original unit cell vectors
            miller_vectors (dict): Miller index-based unit cell vectors
            
        Returns:
            dict: Verification results
        """
        print("Starting Miller index-based unit cell dedicated verification")
        
        verification_results = {
            'element_consistency': False,
            'bond_pattern_consistency': False,
            'volume_reasonability': False,
            'geometric_consistency': False,
            'overall_status': False
        }
        
        try:
            # 1. Element consistency verification
            verification_results['element_consistency'] = self._verify_element_consistency(
                original_unit_atoms, miller_unit_atoms
            )
            
            # 2. Bond pattern consistency verification
            verification_results['bond_pattern_consistency'] = self._verify_bond_pattern_consistency(
                original_unit_atoms, miller_unit_atoms, original_vectors, miller_vectors
            )
            
            # 3. Volume reasonability verification
            verification_results['volume_reasonability'] = self._verify_volume_reasonability(
                original_vectors, miller_vectors
            )
            
            # 4. Geometric consistency verification
            verification_results['geometric_consistency'] = self._verify_geometric_consistency(
                original_unit_atoms, miller_unit_atoms, original_vectors, miller_vectors
            )
            
            # Overall verification result
            verification_results['overall_status'] = all([
                verification_results['element_consistency'],
                verification_results['bond_pattern_consistency'],
                verification_results['volume_reasonability'],
                verification_results['geometric_consistency']
            ])
            
            self._display_verification_details(verification_results)
            
            return verification_results
            
        except Exception as e:
            print(f"ERROR: Error during Miller unit cell verification: {str(e)}")
            return verification_results
    
    def _verify_element_consistency(self, original_atoms, miller_atoms):
        """
        Verify element composition consistency.
        
        Args:
            original_atoms (list): Original atoms
            miller_atoms (list): Miller atoms
            
        Returns:
            bool: Element consistency verification result
        """
        print("\nElement consistency verification:")
        
        # Calculate element counts
        original_elements = ElementCounter.count_elements(original_atoms)
        miller_elements = ElementCounter.count_elements(miller_atoms)
        
        print(f"   Original unit cell: {original_elements}")
        print(f"   Miller unit cell: {miller_elements}")
        
        # Compare element types and ratios
        element_match = self._compare_element_ratios(original_elements, miller_elements)
        
        print(f"   Element consistency: {'SUCCESS' if element_match else 'ERROR'}")
        
        return element_match
    
    def _verify_bond_pattern_consistency(self, original_atoms, miller_atoms, original_vectors, miller_vectors):
        """
        Verify bond pattern consistency.
        
        Args:
            original_atoms (list): Original atoms
            miller_atoms (list): Miller atoms
            original_vectors (dict): Original vectors
            miller_vectors (dict): Miller vectors
            
        Returns:
            bool: Bond pattern consistency verification result
        """
        print("\nBond pattern consistency verification:")
        
        try:
            # Analyze bond patterns of original unit cell
            original_bonds = self.bond_calc.find_bonds(original_atoms, original_vectors, max_bonds_per_atom=12)
            original_patterns = self._extract_bond_patterns_from_bonds(original_atoms, original_bonds)
            
            # Analyze bond patterns of Miller unit cell
            miller_bonds = self.bond_calc.find_bonds(miller_atoms, miller_vectors, max_bonds_per_atom=12)
            miller_patterns = self._extract_bond_patterns_from_bonds(miller_atoms, miller_bonds)
            
            # Compare patterns
            pattern_match = self._compare_bond_patterns(original_patterns, miller_patterns)
            
            print(f"   Bond pattern consistency: {'SUCCESS' if pattern_match else 'ERROR'}")
            
            return pattern_match
            
        except Exception as e:
            print(f"   ERROR: Bond pattern verification error: {str(e)}")
            return False
    
    def _verify_volume_reasonability(self, original_vectors, miller_vectors):
        """
        Verify unit cell volume reasonability.
        
        Args:
            original_vectors (dict): Original vectors
            miller_vectors (dict): Miller vectors
            
        Returns:
            bool: Volume reasonability verification result
        """
        print("\nVolume reasonability verification:")
        
        try:
            # Original unit cell volume
            orig_a = np.array(original_vectors['a'])
            orig_b = np.array(original_vectors['b'])
            orig_c = np.array(original_vectors['c'])
            original_volume = abs(np.dot(orig_a, np.cross(orig_b, orig_c)))
            
            # Miller unit cell volume
            mill_a = np.array(miller_vectors['a_user'])
            mill_b = np.array(miller_vectors['b_user'])
            mill_c = np.array(miller_vectors['c_user'])
            miller_volume = abs(np.dot(mill_a, np.cross(mill_b, mill_c)))
            
            print(f"   Original volume: {original_volume:.4f} Å³")
            print(f"   Miller volume: {miller_volume:.4f} Å³")
            
            if original_volume > 0:
                volume_ratio = miller_volume / original_volume
                print(f"   Volume ratio: {volume_ratio:.4f}")
                
                # Reasonable range: 0.25 ~ 16.0 (large variations possible depending on Miller indices)
                reasonable = 0.25 <= volume_ratio <= 16.0
                print(f"   Reasonable range (0.25~16.0): {'SUCCESS' if reasonable else 'ERROR'}")
                
                return reasonable
            else:
                print("   ERROR: Original volume calculation error")
                return False
                
        except Exception as e:
            print(f"   ERROR: Volume verification error: {str(e)}")
            return False
    
    def _verify_geometric_consistency(self, original_atoms, miller_atoms, original_vectors, miller_vectors):
        """
        Verify geometric consistency.
        
        Args:
            original_atoms (list): Original atoms
            miller_atoms (list): Miller atoms
            original_vectors (dict): Original vectors
            miller_vectors (dict): Miller vectors
            
        Returns:
            bool: Geometric consistency verification result
        """
        print("\nGeometric consistency verification:")
        
        try:
            # Calculate atomic density
            orig_volume = self._calculate_volume(original_vectors)
            miller_volume = self._calculate_volume_miller(miller_vectors)
            
            if orig_volume > 0 and miller_volume > 0:
                orig_density = len(original_atoms) / orig_volume
                miller_density = len(miller_atoms) / miller_volume
                
                print(f"   Original density: {orig_density:.4f} atoms/Å³")
                print(f"   Miller density: {miller_density:.4f} atoms/Å³")
                
                # Density should be similar (within ±50%)
                if orig_density > 0:
                    density_ratio = miller_density / orig_density
                    print(f"   Density ratio: {density_ratio:.4f}")
                    
                    reasonable_density = 0.5 <= density_ratio <= 2.0
                    print(f"   Reasonable density range (0.5~2.0): {'SUCCESS' if reasonable_density else 'ERROR'}")
                    
                    return reasonable_density
                else:
                    print("   ERROR: Original density calculation error")
                    return False
            else:
                print("   ERROR: Volume calculation error")
                return False
                
        except Exception as e:
            print(f"   ERROR: Geometric consistency verification error: {str(e)}")
            return False
    
    def _compare_element_ratios(self, original_elements, miller_elements):
        """
        Compare element ratios between two structures.
        
        Args:
            original_elements (dict): Original element counts
            miller_elements (dict): Miller element counts
            
        Returns:
            bool: Whether element ratios are similar
        """
        try:
            # Check if both structures have same element types
            orig_types = set(original_elements.keys())
            miller_types = set(miller_elements.keys())
            
            if orig_types != miller_types:
                print(f"   Element types differ: Original {orig_types} vs Miller {miller_types}")
                return False
            
            # Compare ratios
            orig_total = sum(original_elements.values())
            miller_total = sum(miller_elements.values())
            
            for element in orig_types:
                orig_ratio = original_elements[element] / orig_total
                miller_ratio = miller_elements[element] / miller_total
                
                # Allow ±20% difference
                if abs(orig_ratio - miller_ratio) > 0.2:
                    print(f"   Ratio differs for {element}: {orig_ratio:.3f} vs {miller_ratio:.3f}")
                    return False
            
            print(f"   Element ratios match within tolerance")
            return True
            
        except Exception as e:
            print(f"   ERROR: Element ratio comparison error: {str(e)}")
            return False
    
    def _extract_bond_patterns_from_bonds(self, atoms, bonds_dict):
        """
        Extract bond patterns from bond information.
        
        Args:
            atoms (list): Atom list
            bonds_dict (dict): Bond information dictionary
            
        Returns:
            dict: Bond patterns per element
        """
        try:
            patterns = {}
            
            for i, atom in enumerate(atoms):
                element = atom['element']
                if element not in patterns:
                    patterns[element] = []
                
                # Get bonds for this atom
                atom_bonds = bonds_dict.get(i, [])
                
                # Calculate bond pattern (coordination number and bond distances)
                bond_distances = []
                for bond in atom_bonds:
                    bond_distances.append(bond['distance'])
                
                bond_distances.sort()  # Sort for consistent comparison
                
                patterns[element].append({
                    'coordination': len(bond_distances),
                    'distances': bond_distances[:6]  # Keep up to 6 shortest bonds for comparison
                })
            
            return patterns
            
        except Exception as e:
            print(f"   ERROR: Bond pattern extraction error: {str(e)}")
            return {}
    
    def _compare_bond_patterns(self, original_patterns, miller_patterns):
        """
        Compare bond patterns between two structures.
        
        Args:
            original_patterns (dict): Original bond patterns
            miller_patterns (dict): Miller bond patterns
            
        Returns:
            bool: Whether bond patterns are similar
        """
        try:
            # Check if both have same elements
            if set(original_patterns.keys()) != set(miller_patterns.keys()):
                print(f"   Pattern elements differ")
                return False
            
            for element in original_patterns.keys():
                orig_element_patterns = original_patterns[element]
                miller_element_patterns = miller_patterns[element]
                
                # Get average patterns for comparison
                orig_avg = self._get_average_pattern(orig_element_patterns)
                miller_avg = self._get_average_pattern(miller_element_patterns)
                
                # Compare patterns
                if not self._patterns_similar(orig_avg, miller_avg):
                    print(f"   Bond patterns differ for {element}")
                    return False
            
            print(f"   Bond patterns match within tolerance")
            return True
            
        except Exception as e:
            print(f"   ERROR: Bond pattern comparison error: {str(e)}")
            return False
    
    def _get_average_pattern(self, patterns):
        """
        Calculate average pattern from multiple patterns.
        
        Args:
            patterns (list): List of patterns
            
        Returns:
            dict: Average pattern
        """
        if not patterns:
            return {'coordination': 0, 'distances': []}
        
        avg_coordination = sum(p['coordination'] for p in patterns) / len(patterns)
        
        # Average distances (only for patterns with same coordination)
        max_coord = max(p['coordination'] for p in patterns)
        same_coord_patterns = [p for p in patterns if p['coordination'] == max_coord]
        
        if same_coord_patterns:
            avg_distances = []
            max_dist_len = max(len(p['distances']) for p in same_coord_patterns)
            
            for i in range(max_dist_len):
                distances_at_i = [p['distances'][i] for p in same_coord_patterns if i < len(p['distances'])]
                if distances_at_i:
                    avg_distances.append(sum(distances_at_i) / len(distances_at_i))
        else:
            avg_distances = []
        
        return {
            'coordination': avg_coordination,
            'distances': avg_distances
        }
    
    def _patterns_similar(self, pattern1, pattern2, tolerance=0.5):
        """
        Check if two patterns are similar.
        
        Args:
            pattern1 (dict): First pattern
            pattern2 (dict): Second pattern
            tolerance (float): Tolerance for comparison
            
        Returns:
            bool: Whether patterns are similar
        """
        # Compare coordination numbers (allow ±1 difference)
        coord_diff = abs(pattern1['coordination'] - pattern2['coordination'])
        if coord_diff > 1:
            return False
        
        # Compare bond distances (allow tolerance)
        dist1 = pattern1['distances']
        dist2 = pattern2['distances']
        
        min_len = min(len(dist1), len(dist2))
        if min_len == 0:
            return True
        
        for i in range(min_len):
            if abs(dist1[i] - dist2[i]) > tolerance:
                return False
        
        return True
    
    def _calculate_volume(self, vectors):
        """Calculate volume from original vectors"""
        a, b, c = vectors['a'], vectors['b'], vectors['c']
        return abs(np.dot(a, np.cross(b, c)))
    
    def _calculate_volume_miller(self, vectors):
        """Calculate volume from Miller vectors"""
        a, b, c = vectors['a_user'], vectors['b_user'], vectors['c_user']
        return abs(np.dot(a, np.cross(b, c)))
    
    def _display_verification_details(self, results):
        """
        Display detailed verification results.
        
        Args:
            results (dict): Verification results
        """
        print("\nMiller unit cell verification results:")
        print("=" * 50)
        
        status_map = {True: "SUCCESS", False: "ERROR"}
        
        print(f"Element consistency: {status_map[results['element_consistency']]}")
        print(f"Bond pattern consistency: {status_map[results['bond_pattern_consistency']]}")
        print(f"Volume reasonability: {status_map[results['volume_reasonability']]}")
        print(f"Geometric consistency: {status_map[results['geometric_consistency']]}")
        print("-" * 50)
        print(f"Overall verification: {status_map[results['overall_status']]}")
        print("=" * 50) 