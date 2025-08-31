#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase1: CIF Analysis
CIF file analysis utilizing the existing CIFModule from the functions folder
"""

import os
import sys
from typing import Dict, Optional, Tuple

# Add project root directory to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

# Import from functions package (centralized imports)
from functions import (
    CIFProcessor,
    TemporaryFileHandler,
    SupercellProcessor,
    QEWriter,
    ComponentFactory,
    Phase1DisplayManager
)


class Phase1CifAnalysis:
    """Phase1: CIF Analysis - Utilizing existing CIFModule from functions"""
    
    def __init__(self, factory=None):
        """
        Initialize Phase1CifAnalysis with factory pattern.

        Args:
            factory: ComponentFactory instance (optional, uses default if None)
        """
        # Use factory pattern for component creation
        self.factory = factory or ComponentFactory()
        components = self.factory.create_all_phase1_components()

        # Initialize components from factory
        self.cif_module = components['cif_processor']
        self.temp_handler = components['temp_handler']
        self.supercell_processor = components['supercell_processor']
        self.qe_writer = components['qe_writer']
        self.display_manager = components['display_manager']
    
    def get_user_choice_after_analysis(self) -> str:
        """
        Get user choice for what to do after CIF analysis.

        Returns:
            str: User choice ('direct_qe', 'supercell_qe', 'continue_phase2', 'exit')
        """
        while True:
            try:
                self.display_manager.show_action_choice_menu()
                choice = input("Enter your choice (1-4): ").strip()
                if choice == '1':
                    return 'direct_qe'
                elif choice == '2':
                    return 'supercell_qe'
                elif choice == '3':
                    return 'continue_phase2'
                elif choice == '4':
                    return 'exit'
                else:
                    self.display_manager.show_invalid_choice()
            except KeyboardInterrupt:
                self.display_manager.show_program_interruption()
                return 'exit'
            except Exception as e:
                self.display_manager.show_input_error(e)
                return 'exit'
    
    def generate_direct_qe_output(self, cif_data: Dict) -> bool:
        """
        Generate QE file directly from CIF unit cell data.
        
        Args:
            cif_data (Dict): CIF analysis data
            
        Returns:
            bool: True if successful, False otherwise
        """
        self.display_manager.show_direct_qe_generation_start()

        try:
            # Extract necessary data
            atoms = cif_data['atoms']
            cell_vectors = cif_data['cell_vectors']
            filename = cif_data.get('filename', 'unknown.cif')
            source_folder = cif_data.get('source_folder', 'INPUT_cif')

            self.display_manager.show_unit_cell_info(
                filename=filename,
                atom_count=len(atoms),
                volume=cif_data.get('volume', 'N/A'),
                source_folder=source_folder
            )
            
            # Prepare structure data for QEWriter
            structure = {
                'atoms': atoms,
                'cell_vectors': cell_vectors,
                'structure_type': 'unitcell'
            }
            
            # Default Miller indices for Phase 1 (no Miller analysis yet)
            miller_indices = {'h': 0, 'k': 0, 'l': 0}
            
            # Use QEWriter with existing folder structure functionality
            result = self.qe_writer.write_qe_input_files(
                structure=structure,
                miller_indices=miller_indices,
                source_folder=source_folder,
                cif_filename=filename,
                supercell_size=None,
                suffix="unitcell"
            )
            
            if result and result.get('success'):
                files_created = result.get('files', {})
                if 'full_input' in files_created:
                    qe_filepath = files_created['full_input']
                    self.display_manager.show_qe_file_creation_success(qe_filepath)
                    return True
                else:
                    self.display_manager.show_qe_file_creation_error("No QE file was created")
                    return False
            else:
                self.display_manager.show_qe_file_creation_error("Failed to create QE file")
                return False

        except Exception as e:
            self.display_manager.show_error_during_generation(e)
            return False
    
    def generate_supercell_qe_output(self, cif_data: Dict) -> bool:
        """
        Generate supercell first, then create QE file.
        
        Args:
            cif_data (Dict): CIF analysis data
            
        Returns:
            bool: True if successful, False otherwise
        """
        self.display_manager.show_supercell_qe_generation_start()

        try:
            # Get supercell size from user
            supercell_size = self.get_supercell_size_input()
            if not supercell_size:
                return False

            nx, ny, nz = supercell_size

            # Extract necessary data
            atoms = cif_data['atoms']
            cell_vectors = cif_data['cell_vectors']
            filename = cif_data.get('filename', 'unknown.cif')
            source_folder = cif_data.get('source_folder', 'INPUT_cif')

            self.display_manager.show_supercell_generation_info(nx, ny, nz, source_folder)
            
            # Create temporary output directory for supercell generation
            temp_output_dir = "output"
            os.makedirs(temp_output_dir, exist_ok=True)
            
            # Generate supercell
            supercell_atoms, supercell_vectors = self.supercell_processor.generate_supercell(
                atoms, cell_vectors, nx, ny, nz, temp_output_dir
            )
            
            # Note: Supercell generation info could be added to display_manager if needed
            self.display_manager.show_supercell_size_info(nx, ny, nz)
            
            # Prepare structure data for QEWriter
            structure = {
                'atoms': supercell_atoms,
                'cell_vectors': supercell_vectors,
                'structure_type': 'supercell'
            }
            
            # Default Miller indices for Phase 1 (no Miller analysis yet)
            miller_indices = {'h': 0, 'k': 0, 'l': 0}
            
            # Prepare supercell size info for filename generation
            supercell_size_dict = {'nx': nx, 'ny': ny, 'nz': nz}
            
            # Use QEWriter with existing folder structure functionality
            result = self.qe_writer.write_qe_input_files(
                structure=structure,
                miller_indices=miller_indices,
                source_folder=source_folder,
                cif_filename=filename,
                supercell_size=supercell_size_dict,
                suffix="supercell"
            )
            
            if result and result.get('success'):
                files_created = result.get('files', {})
                if 'full_input' in files_created:
                    qe_filepath = files_created['full_input']
                    self.display_manager.show_qe_file_creation_success(qe_filepath)
                    return True
                else:
                    self.display_manager.show_qe_file_creation_error("No QE file was created")
                    return False
            else:
                self.display_manager.show_qe_file_creation_error("Failed to create QE file")
                return False

        except Exception as e:
            self.display_manager.show_error_during_supercell_generation(e)
            return False
    

    def get_supercell_size_input(self) -> Optional[Tuple[int, int, int]]:
        """
        Get supercell size input from user.
        
        Returns:
            Optional[Tuple[int, int, int]]: (nx, ny, nz) or None if cancelled
        """
        self.display_manager.show_supercell_input_prompt()
        
        try:
            while True:
                nx_input = input("Enter nx (x-direction multiplier, default=2): ").strip()
                if not nx_input:
                    nx = 2
                else:
                    nx = int(nx_input)
                    if nx < 1:
                        self.display_manager.show_invalid_supercell_size("nx")
                        continue
                break
            
            while True:
                ny_input = input("Enter ny (y-direction multiplier, default=2): ").strip()
                if not ny_input:
                    ny = 2
                else:
                    ny = int(ny_input)
                    if ny < 1:
                        self.display_manager.show_invalid_supercell_size("ny")
                        continue
                break
            
            while True:
                nz_input = input("Enter nz (z-direction multiplier, default=2): ").strip()
                if not nz_input:
                    nz = 2
                else:
                    nz = int(nz_input)
                    if nz < 1:
                        self.display_manager.show_invalid_supercell_size("nz")
                        continue
                break
            
            self.display_manager.show_supercell_size_selected(nx, ny, nz)
            return (nx, ny, nz)

        except KeyboardInterrupt:
            self.display_manager.show_cancellation_message()
            return None
        except ValueError:
            self.display_manager.show_qe_file_creation_error("Invalid input. Please enter integer values.")
            return None
        except Exception as e:
            self.display_manager.show_error_during_generation(e)
            return None
    
    def execute(self) -> Optional[Dict]:
        """
        Execute Phase1 CIF analysis.
        
        Returns:
            Optional[Dict]: Analyzed CIF data with processing choice, None if failed
        """
        print("=" * 80)
        print("PHASE1: CIF Analysis Start")
        print("=" * 80)
        
        try:
            # 1. Discover CIF files
            print("\nSearching for CIF files...")
            cif_files = self.cif_module.find_cif_files()
            
            if not cif_files:
                print("WARNING: No readable CIF files found.")
                return None
            
            # Display number of discovered files
            total_files = sum(len(files) for files in cif_files.values())
            print(f"\nTotal {total_files} CIF files discovered.")
            
            # 2. File selection
            selected_folder, selected_path = self.cif_module.display_file_menu(cif_files)
            
            if selected_folder == "exit":
                print("Program terminated by user.")
                return None
            elif selected_folder == "all":
                print("WARNING: Phase1 can only process single files. Selecting the first file.")
                # Select first file
                first_folder = list(cif_files.keys())[0]
                selected_path = cif_files[first_folder][0]
            
            print(f"Selected file: {os.path.basename(selected_path)}")
            
            # 3. CIF analysis
            print("\nAnalyzing CIF file...")
            cif_data = self.cif_module.analyze_cif(selected_path)
            
            if not cif_data:
                print("ERROR: CIF analysis failed.")
                return None
            
            # Add filename and folder information
            cif_data['filename'] = os.path.basename(selected_path)
            cif_data['source_folder'] = selected_folder
            cif_data['full_path'] = selected_path
            
            # 4. Save temporary files (using TemporaryFileHandler from functions)
            self.temp_handler.save_phase_results("phase1", cif_data)
            
            print("=" * 80)
            print("PHASE1: CIF Analysis Complete")
            print("=" * 80)
            
            # 5. Get user choice for next action
            user_choice = self.get_user_choice_after_analysis()
            
            if user_choice == 'exit':
                print("Program terminated by user.")
                return None
            elif user_choice == 'direct_qe':
                success = self.generate_direct_qe_output(cif_data)
                if success:
                    print("\nDirect QE file generation completed successfully!")
                    cif_data['processing_choice'] = 'direct_qe'
                    cif_data['final_output'] = True
                    return cif_data
                else:
                    print("ERROR: Direct QE file generation failed.")
                    return None
            elif user_choice == 'supercell_qe':
                success = self.generate_supercell_qe_output(cif_data)
                if success:
                    print("\nSupercell QE file generation completed successfully!")
                    cif_data['processing_choice'] = 'supercell_qe'
                    cif_data['final_output'] = True
                    return cif_data
                else:
                    print("ERROR: Supercell QE file generation failed.")
                    return None
            elif user_choice == 'continue_phase2':
                print("\nContinuing to Phase2...")
                cif_data['processing_choice'] = 'continue_phase2'
                cif_data['final_output'] = False
                return cif_data
            else:
                print("ERROR: Unknown user choice.")
                return None
            
        except Exception as e:
            print(f"\nError occurred during PHASE1 execution: {str(e)}")
            print("=" * 80)
            return None

def main():
    """Main function used when running Phase1 independently"""
    phase1 = Phase1CifAnalysis()
    result = phase1.execute()
    
    if result:
        print("\nPhase1 CIF analysis completed successfully!")
        return result
    else:
        print("\nPhase1 CIF analysis failed.")
        return None


if __name__ == "__main__":
    main() 