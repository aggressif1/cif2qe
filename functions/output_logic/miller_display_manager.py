"""
Phase3 Miller index-based output management related functions
"""
from ..utils.miller_display_helper import MillerDisplayHelper
from .display_manager_base import BaseDisplayManager

class Phase3MillerDisplayManager(BaseDisplayManager):
    """Class responsible for Phase3 Miller index-based output management"""
    
    def __init__(self):
        """Initialize Phase3MillerDisplayManager"""
        pass
    
    def display_input_info(self, step3_data):
        """
        Display Step4 input information
        
        Args:
            step3_data (dict): Step3 data
        """
        unit_cell_group = step3_data.get('unit_cell_group')
        expanded_supercell = step3_data.get('expanded_supercell')
        miller_indices = step3_data.get('miller_indices')
        
        print("Input Data:")
        
        MillerDisplayHelper.display_miller_indices(step3_data, default_value='N/A')
        
        if unit_cell_group:
            unit_atoms = unit_cell_group.get('atoms', [])
            print(f"   - Miller index-based unit cell atoms: {len(unit_atoms)}")
        else:
            print("   - Miller index-based unit cell: None")
            
        if expanded_supercell:
            super_atoms = expanded_supercell.get('atoms', [])
            size = expanded_supercell.get('supercell_size', {})
            nx = size.get('nx', 'N/A')
            ny = size.get('ny', 'N/A') 
            nz = size.get('nz', 'N/A')
            print(f"   - Miller index-based supercell atoms: {len(super_atoms)}")
            print(f"   - Supercell size: {nx}×{ny}×{nz}")
        else:
            print("   - Miller index-based supercell: None")
    
    def display_verification_summary(self, result_data):
        """
        Display verification result summary.
        
        Args:
            result_data (dict): verification result data
        """
        verification_results = result_data['verification_results']
        stats = result_data['comparison_statistics']
        
        print("\n" + "=" * 80)
        print("PHASE4 VERIFICATION RESULT SUMMARY")
        print("=" * 80)
        
        print("Verification Results:")
        print(f"   - Atom count consistency: {'PASS' if verification_results['atom_count_consistency'] else 'FAIL'}")
        print(f"   - Element composition consistency: {'PASS' if verification_results['element_composition_consistency'] else 'FAIL'}")
        print(f"   - Geometric consistency: {'PASS' if verification_results['geometric_consistency'] else 'FAIL'}")
        print(f"   - Supercell size validation: {'PASS' if verification_results['supercell_size_validation'] else 'FAIL'}")
        print(f"   - Overall verification: {'SUCCESS' if verification_results['overall_status'] else 'FAIL'}")
        
        print("\nComparison Statistics:")
        print(f"   - Unit cell atoms: {stats['unit_cell_total']}")
        print(f"   - Supercell atoms: {stats['supercell_total']}")
        print(f"   - Expansion factor: {stats['expansion_factor']}x")
        
        size = stats['supercell_size']
        print(f"   - Supercell size: {size.get('na')}×{size.get('nb')}×{size.get('nc')}")
        
        MillerDisplayHelper.display_miller_indices(result_data, prefix="   - Verified ", default_value='N/A')
        
        print("=" * 80)
    
    def display_skip_message(self, reason):
        """
        Display verification skip message.

        Args:
            reason (str): reason for skipping
        """
        BaseDisplayManager.show_warning(reason) 