#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase4 structure processor
Responsible for handling bulk and slab structure processing.
"""

from ..utils.temporary_file_handler import TemporaryFileHandler
from ..output_logic.display_manager_phase4 import Phase4DisplayManager

class StructureProcessor:
    """Structure processor"""

    def __init__(self):
        """Initialize StructureProcessor"""
        from .writers.qe_writer import QEWriter
        
        self.qe_writer = QEWriter()
        self.display_manager = Phase4DisplayManager()
        self.temp_handler = TemporaryFileHandler()
    
    def validate_input_data(self, phase3_data):
        """Validate input data"""
        if not phase3_data:
            print("ERROR: No Phase3 data available.")
            return False
        
        if not phase3_data.get('expanded_supercell'):
            print("ERROR: No expanded supercell data available.")
            return False
        
        return True
    
    def ask_structure_type(self):
        """Ask user to select structure type"""
        print("\nWhich structure would you like to generate?")
        print("   1) Bulk Structure")
        print("      - Uses bulk_group atoms and user-defined supercell vectors")
        print("      - 3D structure with periodic boundary conditions")
        print("   2) Slab Structure")
        print("      - 2D slab structure including surfaces")
        print("      - Structure with vacuum layer (to be implemented)")
        
        while True:
            user_choice = input("\nPlease choose (1: Bulk structure, 2: Slab structure, q: Cancel): ").strip().lower()
            
            if user_choice in ['1', 'bulk']:
                return 'bulk'
            elif user_choice in ['2', 'slab']:
                return 'slab'
            elif user_choice in ['q', 'quit', 'cancel']:
                return None
            else:
                print("ERROR: Invalid input. Please enter 1, 2, or q.")
    
    def process_bulk_structure(self, phase3_data, expanded_supercell, miller_indices):
        """Process bulk structure"""
        print("\nGenerating bulk structure...")
        
        try:
            # Check atom group data
            atom_groups = expanded_supercell.get('atom_groups')
            if not atom_groups:
                print("ERROR: No atom group data available.")
                return None
            
            # Extract bulk_group atoms
            bulk_atoms = atom_groups.get('bulk_group', [])
            if not bulk_atoms:
                print("ERROR: No bulk_group atoms available.")
                return None
            
            # Extract user-defined supercell vectors
            user_supercell_vectors = expanded_supercell.get('user_defined_supercell_vectors')
            if not user_supercell_vectors:
                print("ERROR: No user-defined supercell vectors available.")
                return None
            
            a_user_extended = user_supercell_vectors.get('a_user_extended')
            b_user_extended = user_supercell_vectors.get('b_user_extended')
            c_user_extended = user_supercell_vectors.get('c_user_extended')
            
            print(f"Bulk structure information:")
            print(f"   - bulk_group atom count: {len(bulk_atoms)}")
            print(f"   - a_user_extended: {a_user_extended}")
            print(f"   - b_user_extended: {b_user_extended}")
            print(f"   - c_user_extended: {c_user_extended}")
            
            # Construct bulk structure data
            bulk_structure = {
                'atoms': bulk_atoms,
                'cell_vectors': {
                    'a': a_user_extended,
                    'b': b_user_extended,
                    'c': c_user_extended
                },
                'structure_type': 'bulk',
                'atom_groups': atom_groups,
                'user_defined_supercell_vectors': user_supercell_vectors
            }
            
            # Extract file information from Phase1 data
            source_folder = None
            cif_filename = None
            if phase3_data.get('file_info'):
                source_folder = phase3_data['file_info'].get('source_folder')
                cif_filename = phase3_data['file_info'].get('cif_filename')
            
            # Load file information directly from Phase 1 temporary data
            phase1_data = self.temp_handler.load_phase_results("phase1")
            
            if phase1_data:
                source_folder = phase1_data.get('source_folder')
                cif_filename = phase1_data.get('filename')
                print(f"   Restored file information from Phase 1: {cif_filename} (folder: {source_folder})")
            else:
                print("   WARNING: Phase 1 data not found.")
            
            # Extract supercell size information
            supercell_size = expanded_supercell.get('supercell_size', {})
            
            # Generate Quantum ESPRESSO bulk structure file
            qe_output_result = self.generate_bulk_qe_output(
                bulk_structure, miller_indices, source_folder, cif_filename, supercell_size
            )
            
            # Display QE output result
            if qe_output_result:
                self.display_manager.display_qe_output_summary(qe_output_result)
            
            # Construct result data
            result_data = {
                'structure_type': 'bulk',
                'bulk_structure': bulk_structure,
                'qe_output': qe_output_result,
                'miller_indices': miller_indices,
                'phase3_data': phase3_data
            }
            
            print("SUCCESS: Bulk structure generation completed")
            return result_data
            
        except Exception as e:
            print(f"ERROR: Error during bulk structure processing: {str(e)}")
            return None
    
    def generate_bulk_qe_output(self, bulk_structure, miller_indices, source_folder=None, cif_filename=None, supercell_size=None):
        """Generate Quantum ESPRESSO output file"""
        # Get user permission
        user_consent = input("Would you like to generate output files in Quantum ESPRESSO format? (y/n): ").strip().lower()
        
        if user_consent in ['y', 'yes']:
            qe_result = self.qe_writer.write_qe_input_files(bulk_structure, miller_indices, source_folder, cif_filename, supercell_size, suffix="bulk")
            if qe_result:
                return qe_result
            else:
                print("ERROR: Quantum ESPRESSO output file generation failed")
                return None
        else:
            print("Skipping Quantum ESPRESSO output file generation.")
            return {'skipped': True}
    
    def ask_for_hydrogenation(self):
        """Check whether to proceed with hydrogenation process"""
        print("\nWould you like to proceed with the hydrogenation process?")
        print("   - Adds hydrogen atoms to the surface to saturate dangling bonds.")
        print("   - Phase5 will automatically calculate and add hydrogen positions.")
        
        user_choice = input("Would you like to proceed with the hydrogenation process? (y/n): ").strip().lower()
        
        if user_choice in ['y', 'yes']:
            print("SUCCESS: Proceeding to Phase5 hydrogenation process.")
            return True
        else:
            print("Skipping hydrogenation process. Terminating program.")
            return False 