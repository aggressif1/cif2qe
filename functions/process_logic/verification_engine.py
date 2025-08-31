"""
Phase4 verification engine related functions
"""
import numpy as np
from ..utils.element_counter import ElementCounter

class VerificationEngine:
    """Phase4 verification engine class"""
    
    def __init__(self):
        """Initialize Phase4VerificationEngine"""
        pass
    
    def verify_miller_unit_cell_vs_supercell(self, unit_cell_group, expanded_supercell, miller_indices):
        """
        Verify consistency between Miller index-based unit cell and supercell.
        
        Args:
            unit_cell_group (dict): Miller index-based unit cell
            expanded_supercell (dict): Miller index-based supercell
            miller_indices (dict): Miller index information
            
        Returns:
            dict: Verification results
        """
        print(f"Miller index-based unit cell vs supercell verification:")
        
        verification_results = {
            'atom_count_consistency': False,
            'element_composition_consistency': False,
            'geometric_consistency': False,
            'supercell_size_validation': False,
            'overall_status': False
        }
        
        try:
            unit_atoms = unit_cell_group.get('atoms', [])
            super_atoms = expanded_supercell.get('atoms', [])
            supercell_size = expanded_supercell.get('supercell_size', {})
            
            print(f"   - Unit cell atom count: {len(unit_atoms)} atoms")
            print(f"   - Supercell atom count: {len(super_atoms)} atoms")
            print(f"   - Supercell size: {supercell_size.get('nx', 'N/A')}×{supercell_size.get('ny', 'N/A')}×{supercell_size.get('nz', 'N/A')}")
            
            # 1. Verify atom count consistency
            verification_results['atom_count_consistency'] = self.verify_atom_count_consistency(
                unit_atoms, super_atoms, supercell_size
            )
            
            # 2. Verify element composition consistency  
            verification_results['element_composition_consistency'] = self.verify_element_composition_consistency(
                unit_atoms, super_atoms, supercell_size
            )
            
            # 3. Verify geometric consistency
            verification_results['geometric_consistency'] = self.verify_geometric_consistency(
                unit_cell_group, expanded_supercell
            )
            
            # 4. Verify supercell size validation
            verification_results['supercell_size_validation'] = self.verify_supercell_size_validation(
                supercell_size
            )
            
            # Overall verification result
            verification_results['overall_status'] = all([
                verification_results['atom_count_consistency'],
                verification_results['element_composition_consistency'],
                verification_results['geometric_consistency'],
                verification_results['supercell_size_validation']
            ])
            
            return verification_results
            
        except Exception as e:
            print(f"ERROR: Miller unit cell vs supercell verification error: {str(e)}")
            return verification_results
    
    def verify_atom_count_consistency(self, unit_atoms, super_atoms, supercell_size):
        """
        Verify atom count consistency.
        
        Args:
            unit_atoms (list): Unit cell atoms
            super_atoms (list): Supercell atoms  
            supercell_size (dict): Supercell size
            
        Returns:
            bool: Atom count consistency verification result
        """
        print(f"\nAtom count consistency verification:")
        
        na = supercell_size.get('nx', 1)
        nb = supercell_size.get('ny', 1)
        nc = supercell_size.get('nz', 1)
        
        expected_super_atoms = len(unit_atoms) * na * nb * nc
        actual_super_atoms = len(super_atoms)
        
        print(f"   - Unit cell atom count: {len(unit_atoms)} atoms")
        print(f"   - Supercell size: {na}×{nb}×{nc}")
        print(f"   - Expected supercell atom count: {expected_super_atoms} atoms")
        print(f"   - Actual supercell atom count: {actual_super_atoms} atoms")
        
        atom_count_match = expected_super_atoms == actual_super_atoms
        print(f"   - Atom count match: {'SUCCESS: Passed' if atom_count_match else 'ERROR: Failed'}")
        
        return atom_count_match
    
    def verify_element_composition_consistency(self, unit_atoms, super_atoms, supercell_size):
        """
        Verify element composition consistency.
        
        Args:
            unit_atoms (list): Unit cell atoms
            super_atoms (list): Supercell atoms
            supercell_size (dict): Supercell size
            
        Returns:
            bool: Element composition consistency verification result
        """
        print(f"\nElement composition consistency verification:")
        
        # Unit cell element count
        unit_elements = ElementCounter.count_elements(unit_atoms)
        super_elements = ElementCounter.count_elements(super_atoms)
        
        print(f"   - Unit cell element composition: {unit_elements}")
        print(f"   - Supercell element composition: {super_elements}")
        
        # Expected element composition considering supercell size
        na = supercell_size.get('nx', 1)
        nb = supercell_size.get('ny', 1)
        nc = supercell_size.get('nz', 1)
        expansion_factor = na * nb * nc
        
        expected_super_elements = {}
        for element, count in unit_elements.items():
            expected_super_elements[element] = count * expansion_factor
        
        print(f"   - Expected supercell element composition: {expected_super_elements}")
        
        element_match = expected_super_elements == super_elements
        print(f"   - Element composition match: {'SUCCESS: Passed' if element_match else 'ERROR: Failed'}")
        
        return element_match
    
    def verify_geometric_consistency(self, unit_cell_group, expanded_supercell):
        """
        Verify geometric consistency.
        
        Args:
            unit_cell_group (dict): Unit cell group
            expanded_supercell (dict): Expanded supercell
            
        Returns:
            bool: Geometric consistency verification result
        """
        print(f"\nGeometric consistency verification:")
        
        try:
            # Unit cell vectors
            unit_a = np.array(unit_cell_group['a_user'])
            unit_b = np.array(unit_cell_group['b_user'])  
            unit_c = np.array(unit_cell_group['c_user'])
            unit_volume = abs(np.dot(unit_a, np.cross(unit_b, unit_c)))
            
            # Supercell vectors (expected)
            supercell_size = expanded_supercell.get('supercell_size', {})
            na = supercell_size.get('nx', 1)
            nb = supercell_size.get('ny', 1)
            nc = supercell_size.get('nz', 1)
            
            expected_super_volume = unit_volume * na * nb * nc
            
            print(f"   - Unit cell volume: {unit_volume:.4f} Å³")
            print(f"   - Supercell expansion factor: {na}×{nb}×{nc} = {na*nb*nc}")
            print(f"   - Expected supercell volume: {expected_super_volume:.4f} Å³")
            
            # Check volume ratio (±5% tolerance)
            volume_ratio = expected_super_volume / unit_volume if unit_volume > 0 else 0
            expected_ratio = na * nb * nc
            ratio_diff = abs(volume_ratio - expected_ratio) / expected_ratio if expected_ratio > 0 else 1
            
            geometric_consistent = ratio_diff <= 0.05  # 5% tolerance
            print(f"   - Volume ratio: {volume_ratio:.4f} (expected: {expected_ratio})")
            print(f"   - Geometric consistency: {'SUCCESS: Passed' if geometric_consistent else 'ERROR: Failed'}")
            
            return geometric_consistent
            
        except Exception as e:
            print(f"   - ERROR: Geometric verification error: {str(e)}")
            return False
    
    def verify_supercell_size_validation(self, supercell_size):
        """
        Verify validity of supercell size.
        
        Args:
            supercell_size (dict): Supercell size
            
        Returns:
            bool: Supercell size validity verification result
        """
        print(f"\nSupercell size validity verification:")
        
        na = supercell_size.get('nx', 0)
        nb = supercell_size.get('ny', 0)
        nc = supercell_size.get('nz', 0)
        
        print(f"   - Supercell size: {na}×{nb}×{nc}")
        
        # Check if each direction is positive
        size_valid = na > 0 and nb > 0 and nc > 0
        
        # Check if within reasonable range (1~20)
        size_reasonable = 1 <= na <= 20 and 1 <= nb <= 20 and 1 <= nc <= 20
        
        overall_valid = size_valid and size_reasonable
        
        print(f"   - Size validity (positive): {'SUCCESS' if size_valid else 'ERROR'}")
        print(f"   - Size reasonableness (1~20): {'SUCCESS' if size_reasonable else 'ERROR'}")
        print(f"   - Overall validity: {'SUCCESS: Passed' if overall_valid else 'ERROR: Failed'}")
        
        return overall_valid
    
 