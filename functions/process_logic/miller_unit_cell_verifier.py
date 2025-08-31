import numpy as np
from .supercell_verifier import SupercellVerifier
from .bond_calculator import BondCalculator
from ..utils.element_counter import ElementCounter

class MillerUnitCellVerifier:
    """Dedicated verification class for Miller index-based unit cells"""
    
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
        print("🔍 Starting Miller index-based unit cell dedicated verification")
        
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
            print(f"❌ Error during Miller unit cell verification: {str(e)}")
            return verification_results
    
    def _verify_element_consistency(self, original_atoms, miller_atoms):
        """
        Verify consistency of element composition.
        
        Args:
            original_atoms (list): Original atoms
            miller_atoms (list): Miller atoms
            
        Returns:
            bool: Element consistency verification result
        """
        print("\n📊 Element Consistency Verification:")
        
        # Calculate element counts
        original_elements = ElementCounter.count_elements(original_atoms)
        miller_elements = ElementCounter.count_elements(miller_atoms)
        
        print(f"   Original unit cell: {original_elements}")
        print(f"   Miller unit cell: {miller_elements}")
        
        # Compare element types and ratios
        element_match = self._compare_element_ratios(original_elements, miller_elements)
        
        print(f"   Element consistency: {'✅ Pass' if element_match else '❌ Fail'}")
        
        return element_match
    
    def _verify_bond_pattern_consistency(self, original_atoms, miller_atoms, original_vectors, miller_vectors):
        """
        Verify consistency of bond patterns.
        
        Args:
            original_atoms (list): Original atoms
            miller_atoms (list): Miller atoms
            original_vectors (dict): Original vectors
            miller_vectors (dict): Miller vectors
            
        Returns:
            bool: Bond pattern consistency verification result
        """
        print("\n🔗 Bond Pattern Consistency Verification:")
        
        try:
            # Analyze bond patterns of original unit cell
            original_bonds = self.bond_calc.find_bonds(original_atoms, original_vectors, max_bonds_per_atom=12)
            original_patterns = self._extract_bond_patterns_from_bonds(original_atoms, original_bonds)
            
            # Analyze bond patterns of Miller unit cell
            miller_bonds = self.bond_calc.find_bonds(miller_atoms, miller_vectors, max_bonds_per_atom=12)
            miller_patterns = self._extract_bond_patterns_from_bonds(miller_atoms, miller_bonds)
            
            # Compare patterns
            pattern_match = self._compare_bond_patterns(original_patterns, miller_patterns)
            
            print(f"   Bond pattern consistency: {'✅ Pass' if pattern_match else '❌ Fail'}")
            
            return pattern_match
            
        except Exception as e:
            print(f"   ❌ Bond pattern verification error: {str(e)}")
            return False
    
    def _verify_volume_reasonability(self, original_vectors, miller_vectors):
        """
        Verify reasonability of unit cell volume.
        
        Args:
            original_vectors (dict): Original vectors
            miller_vectors (dict): Miller vectors
            
        Returns:
            bool: Volume reasonability verification result
        """
        print("\n📦 Volume Reasonability Verification:")
        
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
                
                # Reasonable range: 0.25 ~ 16.0 (large variation possible depending on Miller indices)
                reasonable = 0.25 <= volume_ratio <= 16.0
                print(f"   Reasonable range (0.25~16.0): {'✅ Pass' if reasonable else '❌ Fail'}")
                
                return reasonable
            else:
                print("   ❌ Original volume calculation error")
                return False
                
        except Exception as e:
            print(f"   ❌ Volume verification error: {str(e)}")
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
        print("\n📐 Geometric Consistency Verification:")
        
        try:
            # Calculate atomic density
            orig_volume = self._calculate_volume(original_vectors)
            miller_volume = self._calculate_volume_miller(miller_vectors)
            
            orig_density = len(original_atoms) / orig_volume if orig_volume > 0 else 0
            miller_density = len(miller_atoms) / miller_volume if miller_volume > 0 else 0
            
            print(f"   Original density: {orig_density:.6f} atoms/Å³")
            print(f"   Miller density: {miller_density:.6f} atoms/Å³")
            
            # Compare density (±50% tolerance)
            if orig_density > 0:
                density_ratio = miller_density / orig_density
                density_reasonable = 0.5 <= density_ratio <= 2.0
                print(f"   Density ratio: {density_ratio:.4f}")
                print(f"   Density reasonability (0.5~2.0): {'✅ Pass' if density_reasonable else '❌ Fail'}")
                
                return density_reasonable
            else:
                print("   ⚠️  Density calculation not possible")
                return True  # Rely on other verifications
                
        except Exception as e:
            print(f"   ❌ Geometric verification error: {str(e)}")
            return False
    

    
    def _compare_element_ratios(self, original_elements, miller_elements):
        """Compare element ratios."""
        # Check if element types are the same
        if set(original_elements.keys()) != set(miller_elements.keys()):
            return False
        
        # Calculate and compare ratios
        for element in original_elements:
            orig_count = original_elements[element]
            miller_count = miller_elements[element]
            
            # Simplify ratios using greatest common divisor
            from math import gcd
            element_gcd = gcd(orig_count, miller_count)
            orig_ratio = orig_count // element_gcd
            miller_ratio = miller_count // element_gcd
            
            # Reasonable ratio range (1:1 ~ 1:8)
            max_ratio = max(orig_ratio, miller_ratio)
            min_ratio = min(orig_ratio, miller_ratio)
            
            if min_ratio == 0 or max_ratio / min_ratio > 8:
                return False
        
        return True
    
    def _extract_bond_patterns_from_bonds(self, atoms, bonds_dict):
        """Extract patterns from bond information."""
        patterns = {}
        
        for i, atom in enumerate(atoms):
            element = atom.get('element', atom.get('atom_type', 'Unknown'))
            atom_bonds = bonds_dict.get(i, [])
            
            # Calculate count and average distance by bonding partner
            partner_counts = {}
            partner_distances = {}
            
            for bond in atom_bonds:
                partner = bond['partner_element']
                distance = bond['distance']
                
                if partner not in partner_counts:
                    partner_counts[partner] = 0
                    partner_distances[partner] = []
                
                partner_counts[partner] += 1
                partner_distances[partner].append(distance)
            
            # Calculate average distances
            avg_distances = {}
            for partner, distances in partner_distances.items():
                avg_distances[partner] = sum(distances) / len(distances)
            
            pattern = {
                'coordination': len(atom_bonds),
                'partners': partner_counts,
                'avg_distances': avg_distances
            }
            
            if element not in patterns:
                patterns[element] = []
            patterns[element].append(pattern)
        
        return patterns
    
    def _compare_bond_patterns(self, original_patterns, miller_patterns):
        """Compare bond patterns."""
        # Compare patterns by element
        for element in original_patterns:
            if element not in miller_patterns:
                return False
            
            orig_element_patterns = original_patterns[element]
            miller_element_patterns = miller_patterns[element]
            
            # Compare pattern counts (considering ratios)
            orig_count = len(orig_element_patterns)
            miller_count = len(miller_element_patterns)
            
            if orig_count == 0 or miller_count == 0:
                continue
            
            # Extract and compare representative patterns
            orig_avg_pattern = self._get_average_pattern(orig_element_patterns)
            miller_avg_pattern = self._get_average_pattern(miller_element_patterns)
            
            if not self._patterns_similar(orig_avg_pattern, miller_avg_pattern):
                return False
        
        return True
    
    def _get_average_pattern(self, patterns):
        """Calculate average of patterns."""
        if not patterns:
            return {}
        
        total_coordination = sum(p['coordination'] for p in patterns)
        avg_coordination = total_coordination / len(patterns)
        
        # Find most common partner composition
        all_partners = {}
        for pattern in patterns:
            for partner, count in pattern['partners'].items():
                if partner not in all_partners:
                    all_partners[partner] = []
                all_partners[partner].append(count)
        
        avg_partners = {}
        for partner, counts in all_partners.items():
            avg_partners[partner] = sum(counts) / len(counts)
        
        return {
            'coordination': avg_coordination,
            'partners': avg_partners
        }
    
    def _patterns_similar(self, pattern1, pattern2, tolerance=0.5):
        """Compare whether two patterns are similar."""
        # Compare coordination numbers
        coord_diff = abs(pattern1.get('coordination', 0) - pattern2.get('coordination', 0))
        if coord_diff > tolerance * 2:
            return False
        
        # Compare partner composition
        partners1 = pattern1.get('partners', {})
        partners2 = pattern2.get('partners', {})
        
        for partner in set(partners1.keys()) | set(partners2.keys()):
            count1 = partners1.get(partner, 0)
            count2 = partners2.get(partner, 0)
            
            if abs(count1 - count2) > tolerance:
                return False
        
        return True
    
    def _calculate_volume(self, vectors):
        """Calculate volume from original vectors."""
        a = np.array(vectors['a'])
        b = np.array(vectors['b'])
        c = np.array(vectors['c'])
        return abs(np.dot(a, np.cross(b, c)))
    
    def _calculate_volume_miller(self, vectors):
        """Calculate volume from Miller vectors."""
        a = np.array(vectors['a_user'])
        b = np.array(vectors['b_user'])
        c = np.array(vectors['c_user'])
        return abs(np.dot(a, np.cross(b, c)))
    
    def _display_verification_details(self, results):
        """Display detailed verification results."""
        print("\n📋 Miller Unit Cell Verification Details:")
        print(f"   📊 Element consistency: {'✅' if results['element_consistency'] else '❌'}")
        print(f"   🔗 Bond pattern consistency: {'✅' if results['bond_pattern_consistency'] else '❌'}")
        print(f"   📦 Volume reasonability: {'✅' if results['volume_reasonability'] else '❌'}")
        print(f"   📐 Geometric consistency: {'✅' if results['geometric_consistency'] else '❌'}")
        print(f"   🎯 Overall verification: {'✅ Success' if results['overall_status'] else '❌ Failed'}") 