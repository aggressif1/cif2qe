"""
Phase3 output management related functions
"""
from ..utils.miller_display_helper import MillerDisplayHelper
from typing import Dict, Optional, Tuple
from .display_manager_base import BaseDisplayManager

class Phase3DisplayManager(BaseDisplayManager):
    """Class responsible for Phase3 output management"""

    def __init__(self):
        """Initialize Phase3DisplayManager"""
        pass

    @staticmethod
    def show_phase3_start():
        """Display Phase3 start message."""
        BaseDisplayManager.show_phase_start(3, "Miller Index-based Supercell Generation")

    @staticmethod
    def show_input_data_info(filename: str, atom_count: int):
        """Display input data information."""
        BaseDisplayManager.show_input_data_info(filename, atom_count)

    @staticmethod
    def show_element_counts(element_counts: Dict[str, int]):
        """Display element counts."""
        BaseDisplayManager.show_element_counts(element_counts)

    @staticmethod
    def show_miller_input_start():
        """Display Miller index input start message."""
        print("\nMiller Index Input:")

    @staticmethod
    def show_plane_analysis_start():
        """Display plane analysis start message."""
        print("\nPlane Analysis:")

    @staticmethod
    def show_plane_info(plane_info: Dict):
        """Display plane information."""
        if plane_info:
            print(f"   - Reference plane found: {plane_info.get('plane_equation', 'N/A')}")
            print(f"   - Miller indices: ({plane_info.get('h', 'N/A')}, {plane_info.get('k', 'N/A')}, {plane_info.get('l', 'N/A')})")

    @staticmethod
    def show_unit_cell_generation_start():
        """Display unit cell generation start message."""
        print("\nMiller Index-based Unit Cell Generation:")

    @staticmethod
    def show_unit_cell_info(unit_cell_atoms: int, unit_cell_volume: float):
        """Display unit cell information."""
        print(f"   - Unit cell atoms: {unit_cell_atoms}")
        print(f"   - Unit cell volume: {unit_cell_volume:.4f} Å³")

    @staticmethod
    def show_supercell_expansion_start():
        """Display supercell expansion start message."""
        print("\nMiller Index-based Supercell Expansion:")

    @staticmethod
    def show_supercell_info(supercell_atoms: int, nx: int, ny: int, nz: int):
        """Display supercell information."""
        print(f"   - Supercell atoms: {supercell_atoms}")
        print(f"   - Expansion size: {nx}×{ny}×{nz}")

    @staticmethod
    def show_verification_start():
        """Display verification start message."""
        print("\nStructure Verification:")

    @staticmethod
    def show_verification_result(success: bool):
        """Display verification result."""
        if success:
            print("   - Verification: PASSED")
        else:
            print("   - Verification: FAILED")

    @staticmethod
    def show_qe_generation_start():
        """Display QE generation start message."""
        print("\nQuantum ESPRESSO File Generation:")

    @staticmethod
    def show_qe_file_info(filename: str):
        """Display QE file information."""
        print(f"   - Generated file: {filename}")

    @staticmethod
    def show_phase3_completion():
        """Display Phase3 completion message."""
        print("\n" + "=" * 80)
        print("PHASE3: Miller Index-based Supercell Generation Complete")
        print("=" * 80)

    @staticmethod
    def show_error_during_phase3(error: Exception):
        """Display error message during Phase3 execution."""
        BaseDisplayManager.show_error(error, "Error occurred during Phase3 execution")
    
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
        print("WARNING: Miller index-based unit cell or supercell was not generated.")
        print("Skipping verification.") 