#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase4 slab structure processor
Responsible for slab structure generation and vacuum layer addition.
"""

import numpy as np
from ..utils.temporary_file_handler import TemporaryFileHandler
from ..utils.coordinate_transformer import CoordinateTransformer
from ..output_logic.display_manager_phase4 import Phase4DisplayManager

class SlabProcessor:
    """Slab structure processor"""

    def __init__(self):
        """Initialize SlabProcessor"""
        from .writers.qe_writer import QEWriter
        
        self.qe_writer = QEWriter()
        self.display_manager = Phase4DisplayManager()
        self.temp_handler = TemporaryFileHandler()
        self.coord_transformer = CoordinateTransformer()
    
    def process_slab_structure(self, phase3_data, expanded_supercell, miller_indices):
        """Process slab structure"""
        print("\n🔧 Creating slab structure...")
        
        try:
            # Check atom group data
            atom_groups = expanded_supercell.get('atom_groups')
            if not atom_groups:
                print("❌ No atom group data.")
                return None
            
            # Extract bulk_group and boundary_group atoms
            bulk_atoms = atom_groups.get('bulk_group', [])
            boundary_atoms = atom_groups.get('boundary_group', [])
            
            # Slab atoms = bulk_group + boundary_group
            slab_atoms = bulk_atoms + boundary_atoms
            
            if not slab_atoms:
                print("❌ No atoms for slab.")
                return None
            
            # Extract user-defined supercell cell vectors
            user_supercell_vectors = expanded_supercell.get('user_defined_supercell_vectors')
            if not user_supercell_vectors:
                print("❌ No user-defined supercell cell vectors.")
                return None
            
            a_user_extended = user_supercell_vectors.get('a_user_extended')
            b_user_extended = user_supercell_vectors.get('b_user_extended')
            c_user_extended = user_supercell_vectors.get('c_user_extended')
            
            print(f"📊 Slab structure information:")
            print(f"   - Slab atom count: {len(slab_atoms)} (bulk: {len(bulk_atoms)}, boundary: {len(boundary_atoms)})")
            print(f"   - a_user_extended: {a_user_extended}")
            print(f"   - b_user_extended: {b_user_extended}")
            print(f"   - c_user_extended: {c_user_extended}")
            
            # Input vacuum layer size
            vacuum_sizes = self.get_vacuum_sizes_for_slab()
            if not vacuum_sizes:
                print("❌ Vacuum layer size input was cancelled.")
                return None
            
            # Convert fractional coordinates to absolute coordinates
            print("\n🔧 Converting fractional coordinates → absolute coordinates...")
            absolute_atoms = self.coord_transformer.convert_fractional_to_absolute_coordinates(
                slab_atoms, a_user_extended, b_user_extended, c_user_extended
            )
            
            if not absolute_atoms:
                print("❌ Failed to convert fractional coordinates → absolute coordinates")
                return None
            
            # Apply vacuum layers (atom movement + cell vector expansion)
            print("\n🔧 Applying vacuum layers...")
            final_atoms, new_cell_vectors = self.apply_vacuum_layers_method1(
                absolute_atoms, vacuum_sizes, a_user_extended, b_user_extended, c_user_extended
            )
            
            if not final_atoms or not new_cell_vectors:
                print("❌ Failed to apply vacuum layers")
                return None
            
            # Compose slab structure data
            slab_structure = {
                'atoms': final_atoms,
                'cell_vectors': new_cell_vectors,
                'structure_type': 'slab',
                'vacuum_applied': True,
                'vacuum_sizes': vacuum_sizes
            }
            
            print(f"✅ Slab structure creation completed: {len(final_atoms)} atoms")
            
            # Generate QE file
            phase1_data = self.temp_handler.load_phase_results("phase1")
            source_folder = None
            cif_filename = None
            
            if phase1_data:
                source_folder = phase1_data.get('source_folder')
                cif_filename = phase1_data.get('filename')
                print(f"   📝 Restored file information from Phase 1: {cif_filename} (folder: {source_folder})")
            
            supercell_size = expanded_supercell.get('supercell_size', {})
            
            qe_output_result = self.generate_slab_qe_output(
                slab_structure, miller_indices, source_folder, cif_filename, supercell_size
            )
            
            if qe_output_result:
                self.display_manager.display_qe_output_summary(qe_output_result)
                slab_structure['qe_output'] = qe_output_result
            
            return {
                'structure_type': 'slab',
                'slab_structure': slab_structure,
                'vacuum_sizes': vacuum_sizes,
                'absolute_atoms': absolute_atoms,
                'final_atoms': final_atoms,
                'new_cell_vectors': new_cell_vectors,
                'atom_groups': atom_groups,
                'user_supercell_vectors': user_supercell_vectors,
                'phase3_data': phase3_data
            }
            
        except Exception as e:
            print(f"❌ Error during slab structure processing: {str(e)}")
            return None
    
    def get_vacuum_sizes_for_slab(self):
        """Input vacuum layer sizes for slab structure"""
        print("\n🔧 Slab structure vacuum layer size configuration")
        print("=" * 50)
        print("Enter vacuum layer sizes in positive and negative directions along user-defined cell vectors.")
        print("Enter how much to expand in each direction in Å (angstrom) units.")
        print("Enter 0 to not add vacuum layers in that direction.")
        print("\n💡 Processing method:")
        print("   - Positive direction: Expand cell vector in that direction")
        print("   - Negative direction: Move all atoms in positive direction + expand cell vector")
        
        vacuum_sizes = {}
        
        try:
            # Vacuum layer in a_user_extended direction
            print(f"\n📐 Vacuum layer in a_user_extended direction:")
            while True:
                try:
                    a_pos_input = input("  Positive direction vacuum layer size (Å): ").strip()
                    if a_pos_input.lower() in ['q', 'quit', 'cancel']:
                        return None
                    a_pos_vacuum = float(a_pos_input)
                    if a_pos_vacuum < 0:
                        print("  ❌ Cannot enter negative numbers. Enter positive number or 0.")
                        continue
                    break
                except ValueError:
                    print("  ❌ Please enter a valid number.")
            
            while True:
                try:
                    a_neg_input = input("  Negative direction vacuum layer size (Å): ").strip()
                    if a_neg_input.lower() in ['q', 'quit', 'cancel']:
                        return None
                    a_neg_vacuum = float(a_neg_input)
                    if a_neg_vacuum < 0:
                        print("  ❌ Cannot enter negative numbers. Enter positive number or 0.")
                        continue
                    break
                except ValueError:
                    print("  ❌ Please enter a valid number.")
            
            vacuum_sizes['a_positive'] = a_pos_vacuum
            vacuum_sizes['a_negative'] = a_neg_vacuum
            
            # Vacuum layer in b_user_extended direction
            print(f"\n📐 Vacuum layer in b_user_extended direction:")
            while True:
                try:
                    b_pos_input = input("  Positive direction vacuum layer size (Å): ").strip()
                    if b_pos_input.lower() in ['q', 'quit', 'cancel']:
                        return None
                    b_pos_vacuum = float(b_pos_input)
                    if b_pos_vacuum < 0:
                        print("  ❌ Cannot enter negative numbers. Enter positive number or 0.")
                        continue
                    break
                except ValueError:
                    print("  ❌ Please enter a valid number.")
            
            while True:
                try:
                    b_neg_input = input("  Negative direction vacuum layer size (Å): ").strip()
                    if b_neg_input.lower() in ['q', 'quit', 'cancel']:
                        return None
                    b_neg_vacuum = float(b_neg_input)
                    if b_neg_vacuum < 0:
                        print("  ❌ Cannot enter negative numbers. Enter positive number or 0.")
                        continue
                    break
                except ValueError:
                    print("  ❌ Please enter a valid number.")
            
            vacuum_sizes['b_positive'] = b_pos_vacuum
            vacuum_sizes['b_negative'] = b_neg_vacuum
            
            # Vacuum layer in c_user_extended direction
            print(f"\n📐 Vacuum layer in c_user_extended direction:")
            while True:
                try:
                    c_pos_input = input("  Positive direction vacuum layer size (Å): ").strip()
                    if c_pos_input.lower() in ['q', 'quit', 'cancel']:
                        return None
                    c_pos_vacuum = float(c_pos_input)
                    if c_pos_vacuum < 0:
                        print("  ❌ Cannot enter negative numbers. Enter positive number or 0.")
                        continue
                    break
                except ValueError:
                    print("  ❌ Please enter a valid number.")
            
            while True:
                try:
                    c_neg_input = input("  Negative direction vacuum layer size (Å): ").strip()
                    if c_neg_input.lower() in ['q', 'quit', 'cancel']:
                        return None
                    c_neg_vacuum = float(c_neg_input)
                    if c_neg_vacuum < 0:
                        print("  ❌ Cannot enter negative numbers. Enter positive number or 0.")
                        continue
                    break
                except ValueError:
                    print("  ❌ Please enter a valid number.")
            
            vacuum_sizes['c_positive'] = c_pos_vacuum
            vacuum_sizes['c_negative'] = c_neg_vacuum
            
            # Confirm input
            print(f"\n📋 Entered vacuum layer sizes:")
            print(f"   a_user_extended direction: +{vacuum_sizes['a_positive']:.2f} Å (positive), +{vacuum_sizes['a_negative']:.2f} Å (negative)")
            print(f"   b_user_extended direction: +{vacuum_sizes['b_positive']:.2f} Å (positive), +{vacuum_sizes['b_negative']:.2f} Å (negative)")
            print(f"   c_user_extended direction: +{vacuum_sizes['c_positive']:.2f} Å (positive), +{vacuum_sizes['c_negative']:.2f} Å (negative)")
            
            while True:
                confirm = input("\nProceed with this? (y/n): ").strip().lower()
                if confirm in ['y', 'yes', 'yes']:
                    return vacuum_sizes
                elif confirm in ['n', 'no', 'no']:
                    print("Restarting vacuum layer size input.")
                    return self.get_vacuum_sizes_for_slab()
                else:
                    print("❌ Please enter y or n.")
            
        except KeyboardInterrupt:
            print("\n❌ User cancelled input.")
            return None
        except Exception as e:
            print(f"❌ Error during vacuum layer size input: {str(e)}")
            return None



    def apply_vacuum_layers_method1(self, atoms, vacuum_sizes, a_user_extended, b_user_extended, c_user_extended):
        """Atom movement + cell vector expansion"""
        try:
            print(f"   Number of atoms to process: {len(atoms)}")
            
            # Calculate original cell vector lengths
            a_length = np.linalg.norm(np.array(a_user_extended))
            b_length = np.linalg.norm(np.array(b_user_extended))
            c_length = np.linalg.norm(np.array(c_user_extended))
            
            # Calculate unit vectors
            a_unit = np.array(a_user_extended) / a_length if a_length > 0 else np.array([1, 0, 0])
            b_unit = np.array(b_user_extended) / b_length if b_length > 0 else np.array([0, 1, 0])
            c_unit = np.array(c_user_extended) / c_length if c_length > 0 else np.array([0, 0, 1])
            
            print(f"   Original cell vector lengths: |a|={a_length:.2f}Å, |b|={b_length:.2f}Å, |c|={c_length:.2f}Å")
            
            # Extract vacuum layer sizes
            a_pos = vacuum_sizes.get('a_positive', 0.0)
            a_neg = vacuum_sizes.get('a_negative', 0.0)
            b_pos = vacuum_sizes.get('b_positive', 0.0)
            b_neg = vacuum_sizes.get('b_negative', 0.0)
            c_pos = vacuum_sizes.get('c_positive', 0.0)
            c_neg = vacuum_sizes.get('c_negative', 0.0)
            
            print(f"   Vacuum layer sizes: a(+{a_pos:.0f}/-{a_neg:.0f}Å), b(+{b_pos:.0f}/-{b_neg:.0f}Å), c(+{c_pos:.0f}/-{c_neg:.0f}Å)")
            
            # Calculate new cell vectors
            new_a_length = a_length + a_pos + a_neg
            new_b_length = b_length + b_pos + b_neg
            new_c_length = c_length + c_pos + c_neg
            
            # Calculate atom movement vector
            shift_vector = (a_neg * a_unit + b_neg * b_unit + c_neg * c_unit)
            
            # Move atoms
            final_atoms = []
            for atom in atoms:
                new_atom = atom.copy()
                abs_pos = np.array([atom['abs_x'], atom['abs_y'], atom['abs_z']])
                new_pos = abs_pos + shift_vector
                
                new_atom['x'] = new_pos[0]
                new_atom['y'] = new_pos[1]
                new_atom['z'] = new_pos[2]
                
                final_atoms.append(new_atom)
            
            # New cell vectors (expanded cell vectors with vacuum layers added)
            new_cell_vectors = {
                'a': (new_a_length * a_unit).tolist(),
                'b': (new_b_length * b_unit).tolist(),
                'c': (new_c_length * c_unit).tolist(),
                # Also add key names needed for Phase5
                'a_user_extended_withvac': (new_a_length * a_unit).tolist(),
                'b_user_extended_withvac': (new_b_length * b_unit).tolist(),
                'c_user_extended_withvac': (new_c_length * c_unit).tolist()
            }
            
            print(f"   ✅ Vacuum layers applied successfully: {len(final_atoms)} atoms moved")
            
            return final_atoms, new_cell_vectors
            
        except Exception as e:
            print(f"   ❌ Error during vacuum layer application: {str(e)}")
            return None, None

    def generate_slab_qe_output(self, slab_structure, miller_indices, source_folder=None, cif_filename=None, supercell_size=None):
        """Generate QE output files for slab structure"""
        try:
            qe_result = self.qe_writer.write_qe_input_files(slab_structure, miller_indices, source_folder, cif_filename, supercell_size, suffix="slab")
            
            if qe_result:
                return qe_result
            else:
                print("❌ Failed to create slab structure QE file")
                return None
                
        except Exception as e:
            print(f"❌ Error during slab structure QE file creation: {str(e)}")
            return None 