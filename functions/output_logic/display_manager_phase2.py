"""
Display Manager for Phase2: Pattern Detection
Handles all UI/output logic for Phase2 operations.
"""

from typing import Dict, Tuple, Optional
from .display_manager_base import BaseDisplayManager


class Phase2DisplayManager(BaseDisplayManager):
    """
    Display Manager for Phase2 Pattern Detection operations.
    Handles Phase2-specific UI and inherits common functionality from BaseDisplayManager.
    """

    @staticmethod
    def show_phase2_start():
        """
        Display Phase2 start message.
        """
        BaseDisplayManager.show_phase_start(2, "Pattern Detection Supercell Generation")

    @staticmethod
    def show_phase1_data_error():
        """
        Display Phase1 data validation error message.
        """
        BaseDisplayManager.show_error(Exception("Phase1 data is not valid"))

    @staticmethod
    def show_input_data_info(filename: str, atom_count: int, volume: str):
        """
        Display input data information.

        Args:
            filename (str): CIF filename
            atom_count (int): Number of atoms in unit cell
            volume (str): Unit cell volume
        """
        BaseDisplayManager.show_input_data_info(filename, atom_count, volume)

    @staticmethod
    def show_supercell_generation_info(nx: int, ny: int, nz: int):
        """
        Display supercell generation information.

        Args:
            nx (int): X-direction multiplier
            ny (int): Y-direction multiplier
            nz (int): Z-direction multiplier
        """
        print(f"\nGenerating {nx}x{ny}x{nz} supercell...")

    @staticmethod
    def show_verification_start():
        """
        Display supercell verification start message.
        """
        print(f"\nVerifying supercell structure...")

    @staticmethod
    def show_verification_warning():
        """
        Display verification warning message.
        """
        print("WARNING: Issues found in supercell structure verification, but continuing.")

    @staticmethod
    def show_csv_save_start():
        """
        Display CSV save start message.
        """
        print(f"\nSaving supercell coordinates...")

    @staticmethod
    def show_phase2_completion(success: bool):
        """
        Display Phase2 completion message.

        Args:
            success (bool): Whether verification passed
        """
        print("\n" + "=" * 80)
        if success:
            print("PHASE2: Pattern Detection Supercell Generation Complete")
        else:
            print("PHASE2: Pattern Detection Supercell Generation Complete (Some issues found in verification)")
        print("=" * 80)

    @staticmethod
    def show_result_summary(filename: str, nx: int, ny: int, nz: int,
                          supercell_atom_count: int, verification_passed: bool,
                          csv_path: Optional[str] = None):
        """
        Display result summary.

        Args:
            filename (str): Original filename
            nx (int): X-direction multiplier
            ny (int): Y-direction multiplier
            nz (int): Z-direction multiplier
            supercell_atom_count (int): Number of atoms in supercell
            verification_passed (bool): Verification result
            csv_path (str, optional): Path to saved CSV file
        """
        print(f"Processed file: {filename}")
        print(f"Supercell size: {nx}×{ny}×{nz}")
        print(f"Supercell atom count: {supercell_atom_count}")
        print(f"Structure verification: {'Passed' if verification_passed else 'Some issues found'}")
        if verification_passed and csv_path:
            print(f"CSV saved: {csv_path}")
        print("=" * 80)

    @staticmethod
    def show_error_during_supercell_generation(error: Exception):
        """
        Display error message during supercell generation.

        Args:
            error (Exception): The error that occurred
        """
        BaseDisplayManager.show_error(error, "Error occurred during supercell generation")
