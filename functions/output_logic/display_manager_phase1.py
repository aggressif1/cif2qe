"""
Display Manager for Phase1: CIF Analysis
Handles all UI/output logic for Phase1 operations.
"""

from typing import Dict, Optional, Tuple
from .display_manager_base import BaseDisplayManager


class Phase1DisplayManager(BaseDisplayManager):
    """
    Display Manager for Phase1 CIF Analysis operations.
    Handles Phase1-specific UI and inherits common functionality from BaseDisplayManager.
    """

    @staticmethod
    def show_action_choice_menu():
        """
        Display the action choice menu for Phase1.
        """
        print("\n" + "=" * 80)
        print("Choose next action:")
        print("=" * 80)
        print("1. Generate QE file with the as-is unit cell")
        print("2. Generate QE file with the supercell of user-defined size")
        print("3. Generate QE file with the supercell of user-defined crystal directions")
        print("4. Exit")
        print("=" * 80)

    @staticmethod
    def show_invalid_choice():
        """
        Display invalid choice message.
        """
        BaseDisplayManager.show_invalid_input("Invalid choice. Please enter 1, 2, 3, or 4.")

    @staticmethod
    def show_program_interruption():
        """
        Display program interruption message.
        """
        BaseDisplayManager.show_program_interruption()

    @staticmethod
    def show_input_error(error: Exception):
        """
        Display input error message.

        Args:
            error (Exception): The error that occurred
        """
        BaseDisplayManager.show_error(error, "Error reading input")

    @staticmethod
    def show_direct_qe_generation_start():
        """
        Display start message for direct QE file generation.
        """
        print("\n" + "=" * 50)
        print("Generating QE file from unit cell...")
        print("=" * 50)

    @staticmethod
    def show_unit_cell_info(filename: str, atom_count: int, volume: str, source_folder: str):
        """
        Display unit cell information.

        Args:
            filename (str): CIF filename
            atom_count (int): Number of atoms
            volume (str): Cell volume
            source_folder (str): Source folder name
        """
        print(f"Unit cell information:")
        print(f"  - Filename: {filename}")
        print(f"  - Number of atoms: {atom_count}")
        print(f"  - Cell volume: {volume} Å³")
        print(f"  - Source folder: {source_folder}")

    @staticmethod
    def show_qe_file_creation_success(filepath: str):
        """
        Display successful QE file creation message.

        Args:
            filepath (str): Path to created QE file
        """
        print(f"QE file successfully created: {filepath}")

    @staticmethod
    def show_qe_file_creation_error(message: str):
        """
        Display QE file creation error message.

        Args:
            message (str): Error message
        """
        print(f"ERROR: {message}")

    @staticmethod
    def show_supercell_qe_generation_start():
        """
        Display start message for supercell QE file generation.
        """
        print("\n" + "=" * 50)
        print("Generating supercell and QE file...")
        print("=" * 50)

    @staticmethod
    def show_supercell_generation_info(nx: int, ny: int, nz: int, source_folder: str):
        """
        Display supercell generation information.

        Args:
            nx (int): X-direction multiplier
            ny (int): Y-direction multiplier
            nz (int): Z-direction multiplier
            source_folder (str): Source folder name
        """
        print(f"Generating {nx}x{ny}x{nz} supercell...")
        print(f"  - Source folder: {source_folder}")

    @staticmethod
    def show_supercell_size_info(nx: int, ny: int, nz: int):
        """
        Display selected supercell size information.

        Args:
            nx (int): X-direction multiplier
            ny (int): Y-direction multiplier
            nz (int): Z-direction multiplier
        """
        print(f"  - Supercell size: {nx}×{ny}×{nz}")

    @staticmethod
    def show_supercell_input_prompt():
        """
        Display supercell size input prompt.
        """
        print("\nEnter supercell size:")

    @staticmethod
    def show_supercell_size_selected(nx: int, ny: int, nz: int):
        """
        Display selected supercell size confirmation.

        Args:
            nx (int): X-direction multiplier
            ny (int): Y-direction multiplier
            nz (int): Z-direction multiplier
        """
        print(f"Selected supercell size: {nx}×{ny}×{nz}")

    @staticmethod
    def show_invalid_supercell_size(axis: str):
        """
        Display invalid supercell size message.

        Args:
            axis (str): Axis name ('nx', 'ny', or 'nz')
        """
        print(f"{axis} must be positive integer")

    @staticmethod
    def show_cancellation_message():
        """
        Display cancellation message.
        """
        BaseDisplayManager.show_cancellation_message()

    @staticmethod
    def show_error_during_generation(error: Exception):
        """
        Display error message during file generation.

        Args:
            error (Exception): The error that occurred
        """
        BaseDisplayManager.show_error(error, "Error during file generation")

    @staticmethod
    def show_error_during_supercell_generation(error: Exception):
        """
        Display error message during supercell generation.

        Args:
            error (Exception): The error that occurred
        """
        BaseDisplayManager.show_error(error, "Error during supercell/QE generation")
