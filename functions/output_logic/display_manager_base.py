"""
Base Display Manager for CIF2QE
Contains common UI/output logic shared across all phases.
"""

import os
from typing import Dict, Optional, Tuple


class BaseDisplayManager:
    """
    Base Display Manager with common UI patterns used across all phases.

    This class provides standardized methods for:
    - Error message display
    - Phase start/completion messages
    - Input data information display
    - Element count display
    - Common output formatting
    """

    @staticmethod
    def show_error(error: Exception, context: str = ""):
        """
        Display standardized error message.

        Args:
            error (Exception): The error that occurred
            context (str, optional): Additional context for the error
        """
        if context:
            print(f"ERROR: {context}: {str(error)}")
        else:
            print(f"ERROR: {str(error)}")

    @staticmethod
    def show_phase_start(phase_num: int, description: str):
        """
        Display standardized phase start message.

        Args:
            phase_num (int): Phase number (1-5)
            description (str): Phase description
        """
        print("=" * 80)
        print(f"PHASE{phase_num}: {description} Start")
        print("=" * 80)

    @staticmethod
    def show_phase_completion(phase_num: int, description: str):
        """
        Display standardized phase completion message.

        Args:
            phase_num (int): Phase number (1-5)
            description (str): Phase description
        """
        print("\n" + "=" * 80)
        print(f"PHASE{phase_num}: {description} Complete")
        print("=" * 80)

    @staticmethod
    def show_input_data_info(filename: str, atom_count: int, volume: str = ""):
        """
        Display standardized input data information.

        Args:
            filename (str): Input filename
            atom_count (int): Number of atoms
            volume (str, optional): Cell volume information
        """
        print("Input data:")
        print(f"   - Filename: {filename}")
        print(f"   - Atom count: {atom_count}")
        if volume:
            print(f"   - Cell volume: {volume} Å³")

    @staticmethod
    def show_element_counts(element_counts: Dict[str, int], indent: str = "   - "):
        """
        Display element count information.

        Args:
            element_counts (Dict[str, int]): Dictionary of element counts
            indent (str): Indentation for display
        """
        if element_counts:
            print(f"{indent}Atoms by element:")
            for element, count in sorted(element_counts.items()):
                print(f"{indent}  {element}: {count}")

    @staticmethod
    def show_separator_line(length: int = 80):
        """
        Display separator line.

        Args:
            length (int): Length of separator line
        """
        print("=" * length)

    @staticmethod
    def show_section_header(title: str, length: int = 80):
        """
        Display section header with borders.

        Args:
            title (str): Section title
            length (int): Length of header line
        """
        print("\n" + "=" * length)
        print(title)
        print("=" * length)

    @staticmethod
    def show_success_message(message: str):
        """
        Display success message.

        Args:
            message (str): Success message
        """
        print(f"SUCCESS: {message}")

    @staticmethod
    def show_warning_message(message: str):
        """
        Display warning message.

        Args:
            message (str): Warning message
        """
        print(f"WARNING: {message}")

    @staticmethod
    def show_info_message(message: str, prefix: str = "   - "):
        """
        Display information message.

        Args:
            message (str): Information message
            prefix (str): Prefix for the message
        """
        print(f"{prefix}{message}")

    @staticmethod
    def show_file_info(filename: str, description: str = "Generated file"):
        """
        Display file information.

        Args:
            filename (str): Filename
            description (str): Description of the file
        """
        print(f"   - {description}: {filename}")

    @staticmethod
    def show_processing_info(message: str):
        """
        Display processing information.

        Args:
            message (str): Processing message
        """
        print(message)

    @staticmethod
    def show_verification_result(success: bool, context: str = "Verification"):
        """
        Display verification result.

        Args:
            success (bool): Verification success status
            context (str): Context for verification
        """
        status = "PASSED" if success else "FAILED"
        print(f"   - {context}: {status}")

    @staticmethod
    def show_cancellation_message():
        """
        Display operation cancellation message.
        """
        print("\nOperation cancelled by user.")

    @staticmethod
    def show_program_interruption():
        """
        Display program interruption message.
        """
        print("\nProgram interrupted by user.")

    @staticmethod
    def show_invalid_input(message: str = "Invalid input"):
        """
        Display invalid input message.

        Args:
            message (str): Invalid input message
        """
        print(f"ERROR: {message}")

    @staticmethod
    def show_file_saved(filename: str, path: Optional[str] = None):
        """
        Display file saved message.

        Args:
            filename (str): Saved filename
            path (str, optional): File path
        """
        if path:
            print(f"SUCCESS: File saved - {path}")
        else:
            print(f"SUCCESS: File saved - {filename}")

    @staticmethod
    def show_miller_indices(data, prefix: str = "   - ", default_value: str = 'N/A'):
        """
        Display Miller indices information.

        Args:
            data (dict): Data containing miller_indices
            prefix (str): Prefix for display
            default_value (str): Default value for missing indices
        """
        miller_indices = data.get('miller_indices')
        if miller_indices:
            h = miller_indices.get('h', default_value)
            k = miller_indices.get('k', default_value)
            l = miller_indices.get('l', default_value)
            print(f"{prefix}Miller indices: ({h}, {k}, {l})")

    @staticmethod
    def display_file_list(files, title: str = "Found Files"):
        """
        Display list of files in hierarchical tree structure.

        Args:
            files (dict): Dictionary of file paths by folder
            title (str): Title to display
        """
        print("\n" + "=" * 80)
        print(" " * (40 - len(title)//2) + title)
        print("=" * 80)

        file_counter = 1

        for folder, paths in files.items():
            print(f"{folder}/")

            if len(paths) == 1:
                filename = os.path.basename(paths[0])
                print(f"└── {file_counter}. {filename}")
                file_counter += 1
            else:
                for i, path in enumerate(paths):
                    filename = os.path.basename(path)
                    if i == len(paths) - 1:
                        print(f"└── {file_counter}. {filename}")
                    else:
                        print(f"├── {file_counter}. {filename}")
                    file_counter += 1
