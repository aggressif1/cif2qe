#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quantum ESPRESSO output generator
Responsible for generating QE files for expanded supercells and final hydrogenated structures.
"""

class QEOutputGenerator:
    """Quantum ESPRESSO output generator"""
    
    def __init__(self):
        """Initialize QEOutputGenerator"""
        from .writers.qe_writer import QEWriter
        from ..output_logic.display_manager_phase5 import Phase5DisplayManager

        self.qe_writer = QEWriter()
        self.display_manager = Phase5DisplayManager()
    
    def generate_supercell_qe_output(self, expanded_supercell, miller_indices, file_info):
        """
        Generate Quantum ESPRESSO files for expanded supercell.
        
        Args:
            expanded_supercell (dict): Expanded supercell information
            miller_indices (dict): Miller index information
            file_info (dict): File information
            
        Returns:
            dict: QE output result
        """
        try:
            print("\nGenerating supercell Quantum ESPRESSO file...")
            
            # Import QE Writer
            from .writers.qe_writer import QEWriter
            qe_writer = QEWriter()
            
            # Convert supercell structure data to QE Writer format
            qe_structure = self._convert_supercell_to_qe_structure(expanded_supercell)
            
            # Load file information directly from Phase 1 temporary data
            from ..utils.temporary_file_handler import TemporaryFileHandler
            temp_handler = TemporaryFileHandler()
            phase1_data = temp_handler.load_phase_results("phase1")
            
            # Extract file information
            source_folder = None
            cif_filename = None
            
            if phase1_data:
                source_folder = phase1_data.get('source_folder')
                cif_filename = phase1_data.get('filename')  # Use filename saved from Phase 1
                print(f"   Restored file information from Phase 1: {cif_filename} (folder: {source_folder})")
            else:
                print("   WARNING: Phase 1 data not found.")
            
            supercell_size = expanded_supercell.get('supercell_size', {})
            
            # Generate QE file
            qe_result = qe_writer.write_qe_input_files(
                qe_structure,
                miller_indices,
                source_folder,
                cif_filename,
                supercell_size
            )
            
            if qe_result and qe_result.get('success'):
                print("SUCCESS: Supercell Quantum ESPRESSO file generation completed")
                print(f"   File location: {qe_result.get('base_filename', 'unknown')}")
                return qe_result
            else:
                print("ERROR: Supercell Quantum ESPRESSO file generation failed")
                return None
                
        except Exception as e:
            print(f"ERROR: Error during supercell QE file generation: {str(e)}")
            return None
    
    def _convert_supercell_to_qe_structure(self, expanded_supercell):
        """
        Convert supercell data to structure format used by QE Writer.
        
        Args:
            expanded_supercell (dict): Expanded supercell information
            
        Returns:
            dict: Structure data in QE Writer format
        """
        try:
            atoms = expanded_supercell.get('atoms', [])
            
            # Convert atom data
            qe_atoms = []
            for atom in atoms:
                qe_atom = {
                    'element': atom.get('element', 'X'),
                    'x': atom.get('x', 0.0),
                    'y': atom.get('y', 0.0),
                    'z': atom.get('z', 0.0)
                }
                qe_atoms.append(qe_atom)
            
            # Convert cell vectors (modified to fit new structure)
            cell_vectors = {
                'a': expanded_supercell.get('a_user', [0, 0, 0]),
                'b': expanded_supercell.get('b_user', [0, 0, 0]),
                'c': expanded_supercell.get('c_user', [0, 0, 0])
            }
            
            qe_structure = {
                'atoms': qe_atoms,
                'cell_vectors': cell_vectors,
                'supercell_info': {
                    'total_atoms': len(qe_atoms),
                    'supercell_size': expanded_supercell.get('supercell_size', {}),
                    'expansion_ratio': expanded_supercell.get('total_atoms', 0) / max(expanded_supercell.get('unit_cell_atoms', 1), 1)
                }
            }
            
            return qe_structure
            
        except Exception as e:
            print(f"ERROR: Error during supercell structure conversion: {str(e)}")
            return None
    
    def display_qe_output_summary(self, result_data):
        """
        Display summary of QE file output results.
        
        Args:
            result_data (dict): Result data
        """
        try:
            supercell_qe_output = result_data.get('supercell_qe_output')
            expanded_supercell = result_data.get('expanded_supercell')
            
            if supercell_qe_output and supercell_qe_output.get('success'):
                print(f"\nSupercell Quantum ESPRESSO File Generation Summary:")
                print(f"   Generated file: {supercell_qe_output.get('base_filename', 'unknown')}.in")
                
                if expanded_supercell:
                    supercell_size = expanded_supercell.get('supercell_size', {})
                    total_atoms = expanded_supercell.get('total_atoms', 0)
                    unit_cell_atoms = expanded_supercell.get('unit_cell_atoms', 0)
                    expansion_ratio = total_atoms / max(unit_cell_atoms, 1) if unit_cell_atoms > 0 else 0
                    
                    print(f"   Supercell size: {supercell_size.get('nx', 1)} × {supercell_size.get('ny', 1)} × {supercell_size.get('nz', 1)}")
                    print(f"   Total atoms: {total_atoms} atoms (unit cell {unit_cell_atoms} → {expansion_ratio:.1f}x expansion)")
                    
                files = supercell_qe_output.get('files', {})
                if files:
                    print(f"   File list:")
                    for file_type, filepath in files.items():
                        filename = filepath.split('/')[-1] if '/' in filepath else filepath.split('\\')[-1]
                        print(f"     - {file_type}: {filename}")
                
                print(f"   Ready for Quantum ESPRESSO calculations!")
                
            elif expanded_supercell:
                print(f"\nWARNING: Supercell was generated but QE file generation failed.")
                supercell_size = expanded_supercell.get('supercell_size', {})
                total_atoms = expanded_supercell.get('total_atoms', 0)
                print(f"   Supercell size: {supercell_size.get('nx', 1)} × {supercell_size.get('ny', 1)} × {supercell_size.get('nz', 1)}")
                print(f"   Total atoms: {total_atoms} atoms")
            else:
                print(f"\nWARNING: Supercell expansion and QE file generation were not performed.")
                
        except Exception as e:
            print(f"WARNING: Error displaying QE output summary: {str(e)}")

    def generate_final_qe_output(self, optimized_structure, miller_indices, source_folder=None, cif_filename=None, supercell_size=None):
        """Generate final QE output file"""
        try:
            print("\nGenerating final Quantum ESPRESSO output file...")

            qe_result = self.qe_writer.write_qe_input_files(optimized_structure, miller_indices, source_folder, cif_filename, supercell_size, suffix="hydrogenated")

            if qe_result:
                print("SUCCESS: Final QE output file generation completed")
                self.display_manager.display_final_qe_summary(qe_result)
                return qe_result
            else:
                print("ERROR: Final QE output file generation failed")
                return None

        except Exception as e:
            print(f"ERROR: Error during final QE output file generation: {str(e)}")
            return None 