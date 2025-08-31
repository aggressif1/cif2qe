"""
Quantum ESPRESSO output-related functions
"""
import os
from datetime import datetime
import numpy as np

class QEWriter:
    """Class responsible for Quantum ESPRESSO output"""
    
    def __init__(self):
        """Initialize QEWriter"""
        self.base_output_dir = "output"
        self.output_dir = self.base_output_dir
        self._ensure_output_directory()
    
    def _ensure_output_directory(self):
        """Ensure output directory exists and create if necessary"""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
    
    def _generate_base_filename(self, cif_filename, miller_indices, supercell_size, structure_type=None, suffix=None):
        """
        Generate filename (unified logic)
        
        Args:
            cif_filename (str, optional): Original CIF filename
            miller_indices (dict): Miller indices information
            supercell_size (dict, optional): Supercell size information
            structure_type (str, optional): Structure type ('bulk', 'slab', etc.)
            suffix (str, optional): Filename suffix (e.g., "slab", "hydrogenated")
            
        Returns:
            str: Generated base filename
        """
        # Generate filename (based on CIF filename)
        if cif_filename:
            # Remove extension from CIF filename
            base_name = cif_filename.replace('.cif', '') if cif_filename.endswith('.cif') else cif_filename
        else:
            # Load filename directly from Phase 1 temporary data (fallback)
            from functions.output_logic.temporary_file_handler import TemporaryFileHandler
            temp_handler = TemporaryFileHandler()
            phase1_data = temp_handler.load_phase_results("phase1")
            
            if phase1_data and phase1_data.get('filename'):
                filename = phase1_data['filename']
                base_name = filename.replace('.cif', '') if filename.endswith('.cif') else filename
                print(f"   📝 Filename restored from Phase 1 (fallback): {filename}")
            else:
                base_name = "unknown"
        
        # Miller indices string
        miller_str = f"{miller_indices.get('h', 0)}{miller_indices.get('k', 0)}{miller_indices.get('l', 0)}"
        
        # Supercell size string
        if supercell_size:
            nx = supercell_size.get('nx', 1)
            ny = supercell_size.get('ny', 1) 
            nz = supercell_size.get('nz', 1)
            supercell_str = f"{nx}x{ny}x{nz}"
        else:
            supercell_str = "1x1x1"
        
        # Determine suffix (priority: suffix parameter > structure_type)
        final_suffix = ""
        if suffix:
            final_suffix = f"_{suffix}"
            print(f"   📝 User-specified suffix: {suffix}")
        elif structure_type == 'bulk':
            final_suffix = "_bulk"
            print(f"   📝 Bulk structure suffix: bulk")
        elif structure_type == 'slab':
            final_suffix = "_slab"
            print(f"   📝 Slab structure suffix: slab")
        
        # Generate final filename: CIFfilename_MillerIndices_SupercellSize_suffix
        base_filename = f"{base_name}_{miller_str}_{supercell_str}{final_suffix}"
        
        return base_filename
    
    def write_qe_input_files(self, structure, miller_indices, source_folder=None, cif_filename=None, supercell_size=None, suffix=None):
        """
        Generate Quantum ESPRESSO input files.
        
        Args:
            structure (dict): Structure data with vacuum layers added or bulk structure data
            miller_indices (dict): Miller indices information
            source_folder (str, optional): Source folder name of original CIF file
            cif_filename (str, optional): Original CIF filename (without extension)
            supercell_size (dict, optional): Supercell size information
            suffix (str, optional): Filename suffix (e.g., "slab", "hydrogenated")
            
        Returns:
            dict: Information about generated files
        """
        try:
            # Check structure type (bulk structure vs others)
            structure_type = structure.get('structure_type', 'unknown')
            
            # Set output directory (INPUT_cif → output)
            if source_folder:
                # Replace INPUT_cif with output
                output_folder = source_folder.replace('INPUT_cif', 'output')
                self.output_dir = os.path.join(self.base_output_dir, output_folder)
                self._ensure_output_directory()
            
            # Generate filename (unified logic)
            if suffix is None and structure_type == 'slab':
                suffix = "slab"
            base_filename = self._generate_base_filename(cif_filename, miller_indices, supercell_size, structure_type, suffix)
            
            # Generate QE input files (main .in file only)
            files_created = {}
            
            # Generate unified QE input file
            full_input_file = self._write_full_input_file(structure, miller_indices, base_filename)
            if full_input_file:
                files_created['full_input'] = full_input_file
            
            if files_created:
                return {
                    'success': True,
                    'files': files_created,
                    'base_filename': base_filename,
                    'structure_type': structure_type
                }
            else:
                print("❌ Failed to generate Quantum ESPRESSO input files")
                return None
                
        except Exception as e:
            print(f"❌ Error during Quantum ESPRESSO input file generation: {str(e)}")
            return None
    

    
    def _write_full_input_file(self, structure, miller_indices, base_filename):
        """Write complete Quantum ESPRESSO input file"""
        try:
            filename = f"{base_filename}.in"
            filepath = os.path.join(self.output_dir, filename)
            
            atoms = structure.get('atoms', [])
            cell_vectors = structure.get('cell_vectors', {})
            vacuum_info = structure.get('vacuum_info', {})
            
            # Calculate element counts
            element_counts = {}
            for atom in atoms:
                element = atom.get('element', 'X')
                element_counts[element] = element_counts.get(element, 0) + 1
            
            nat = len(atoms)
            ntyp = len(element_counts)
            
            with open(filepath, 'w', encoding='utf-8') as f:
                # Header comments
                f.write(f"! Quantum ESPRESSO input file\n")
                f.write(f"! Generated by CIF2QE\n")
                f.write(f"! Miller indices: ({miller_indices.get('h', 0)} {miller_indices.get('k', 0)} {miller_indices.get('l', 0)})\n")
                f.write(f"! Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                
                if vacuum_info:
                    f.write(f"! Vacuum layers added: {vacuum_info.get('coordinate_system', 'unknown')} system\n")
                
                f.write("\n")
                
                # &SYSTEM section
                f.write("&SYSTEM\n")
                f.write("  ibrav = 0\n")
                f.write(f"  nat = {nat}\n")
                f.write(f"  ntyp = {ntyp}\n")
                f.write("/\n\n")
                
                # ATOMIC_POSITIONS card
                f.write("ATOMIC_POSITIONS angstrom\n")
                for atom in atoms:
                    element = atom.get('element', 'X')
                    
                    # Use final coordinates moved due to vacuum layer insertion in slab structure
                    if 'final_x' in atom and 'final_y' in atom and 'final_z' in atom:
                        # Use final coordinates moved due to vacuum layer insertion (slab structure)
                        x = atom.get('final_x')
                        y = atom.get('final_y') 
                        z = atom.get('final_z')
                        
                        # Output debugging information (first 3 atoms only)
                        debug_count = len([a for a in atoms[:atoms.index(atom)+1] if a.get('element') == element])
                        if debug_count <= 3:
                            orig_x = atom.get('abs_x', atom.get('x', 0.0))
                            orig_y = atom.get('abs_y', atom.get('y', 0.0))
                            orig_z = atom.get('abs_z', atom.get('z', 0.0))
                            print(f"   Atom {element}: Original ({orig_x:.6f}, {orig_y:.6f}, {orig_z:.6f}) → Final ({x:.6f}, {y:.6f}, {z:.6f})")
                    
                    # If fractional coordinates exist (bulk structure)
                    elif 'frac_coords_extended' in atom and 'structure_type' in structure and structure['structure_type'] == 'bulk':
                        # Convert fractional coordinates to Cartesian coordinates
                        frac_coords = atom['frac_coords_extended']
                        cell_vectors = structure.get('cell_vectors', {})
                        
                        if cell_vectors and len(frac_coords) == 3:
                            a_vec = np.array(cell_vectors.get('a', [0, 0, 0]))
                            b_vec = np.array(cell_vectors.get('b', [0, 0, 0]))
                            c_vec = np.array(cell_vectors.get('c', [0, 0, 0]))
                            
                            # Convert fractional coordinates → Cartesian coordinates
                            cart_pos = (frac_coords[0] * a_vec + 
                                       frac_coords[1] * b_vec + 
                                       frac_coords[2] * c_vec)
                            
                            x, y, z = cart_pos[0], cart_pos[1], cart_pos[2]
                            

                        else:
                            # Use existing Cartesian coordinates if fractional coordinate conversion fails
                            x = atom.get('x', 0.0)
                            y = atom.get('y', 0.0)
                            z = atom.get('z', 0.0)
                            print(f"   ⚠️ Atom {element}: Fractional coordinate conversion failed, using existing coordinates")
                    else:
                        # Use existing Cartesian coordinates
                        x = atom.get('x', 0.0)
                        y = atom.get('y', 0.0)
                        z = atom.get('z', 0.0)
                    
                    f.write(f"{element:2s} {x:12.8f} {y:12.8f} {z:12.8f}\n")
                f.write("\n")
                
                # CELL_PARAMETERS card
                f.write("CELL_PARAMETERS angstrom\n")
                
                # Handle cell vector key names according to structure type
                # slab structure: a_user_extended_withvac, b_user_extended_withvac, c_user_extended_withvac
                # bulk structure: a, b, c (existing method)
                structure_type = structure.get('structure_type', 'unknown')
                
                if structure_type == 'slab':
                    # For slab structure, use extended cell vectors with vacuum layers added
                    vector_keys = ['a_user_extended_withvac', 'b_user_extended_withvac', 'c_user_extended_withvac']
                else:
                    # For bulk structure, use existing cell vectors
                    vector_keys = ['a', 'b', 'c']
                
                for vector_key in vector_keys:
                    vector = cell_vectors.get(vector_key, [0, 0, 0])
                    f.write(f"{vector[0]:12.8f} {vector[1]:12.8f} {vector[2]:12.8f}\n")
            
            return filepath
            
        except Exception as e:
            print(f"❌ Error writing complete QE input file: {str(e)}")
            return None
    
    def _get_atomic_mass(self, element):
        """Return atomic mass of element (simple version)"""
        atomic_masses = {
            'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
            'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
            'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
            'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
            'Ti': 47.867, 'Fe': 55.845, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38
        }
        return atomic_masses.get(element, 1.0)  # Default value 1.0 