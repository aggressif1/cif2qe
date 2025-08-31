"""
Phase4 output management related functions
"""
from ..utils.element_counter import ElementCounter
from ..utils.miller_display_helper import MillerDisplayHelper
from .display_manager_base import BaseDisplayManager

class Phase4DisplayManager(BaseDisplayManager):
    """Class responsible for Phase4 output management"""
    
    def __init__(self):
        """Initialize Phase4DisplayManager"""
        pass
    
    def display_input_info(self, phase3_data):
        """Display input information received from Phase3"""
        try:
            print("\nPhase4 Input Data:")
            print("-" * 40)
            
            # Miller index information
            MillerDisplayHelper.display_miller_indices(phase3_data)
            
            # Supercell information
            expanded_supercell = phase3_data.get('expanded_supercell', {})
            if expanded_supercell:
                atoms = expanded_supercell.get('atoms', [])
                supercell_size = expanded_supercell.get('supercell_size', {})
                
                print(f"   - Supercell atoms: {len(atoms)}")
                
                if supercell_size:
                    nx = supercell_size.get('nx', 'N/A')
                    ny = supercell_size.get('ny', 'N/A')
                    nz = supercell_size.get('nz', 'N/A')
                    print(f"   - Supercell size: {nx}×{ny}×{nz}")
                
                # Element count
                element_counts = ElementCounter.count_elements(atoms)
                
                if element_counts:
                    print("   - Atoms by element:")
                    for element, count in sorted(element_counts.items()):
                        print(f"     {element}: {count}")
            
            # Unit cell information
            unit_cell_group = phase3_data.get('unit_cell_group', {})
            if unit_cell_group:
                unit_atoms = unit_cell_group.get('atoms', [])
                print(f"   - Unit cell atoms: {len(unit_atoms)}")
                
                # Check user-defined cell vectors
                a_user = unit_cell_group.get('a_user')
                b_user = unit_cell_group.get('b_user')
                c_user = unit_cell_group.get('c_user')
                
                # Check if numpy arrays and not None
                vectors_available = []
                for vector in [a_user, b_user, c_user]:
                    if vector is not None:
                        # If numpy array, check if size is greater than 0
                        if hasattr(vector, '__len__') and len(vector) > 0:
                            vectors_available.append(True)
                        else:
                            vectors_available.append(False)
                    else:
                        vectors_available.append(False)
                
                if all(vectors_available):
                    print("   - User-defined cell vectors available")
                else:
                    print("   - No user-defined cell vectors")
            
        except Exception as e:
            print(f"ERROR: Error displaying input information: {str(e)}")
    
    def display_vacuum_addition_results(self, final_structure, vacuum_config):
        """Display vacuum layer addition results"""
        try:
            print("\nVacuum Layer Addition Results:")
            print("-" * 40)
            
            # Vacuum configuration summary
            coord_system = vacuum_config.get('coordinate_system', 'unknown')
            vacuum_layers = vacuum_config.get('vacuum_layers', {})
            
            if coord_system == 'cartesian':
                print("   Coordinate system: Cartesian coordinates (+x, +y, +z)")
            else:
                print("   Coordinate system: User-defined cell vectors (+a_user, +b_user, +c_user)")
            
            print("   Vacuum layer thickness:")
            total_vacuum = 0
            for direction, thickness in vacuum_layers.items():
                print(f"     +{direction}: {thickness:.2f} Å")
                total_vacuum += thickness
            
            print(f"   Total vacuum thickness: {total_vacuum:.2f} Å")
            
            # Final structure information
            atoms = final_structure.get('atoms', [])
            cell_vectors = final_structure.get('cell_vectors', {})
            vacuum_info = final_structure.get('vacuum_info', {})
            
            print(f"\n   Final structure:")
            print(f"     Number of atoms: {len(atoms)}")
            
            if cell_vectors:
                print("     Cell vectors:")
                for vector_name in ['a', 'b', 'c']:
                    vector = cell_vectors.get(vector_name, [0, 0, 0])
                    # Convert numpy array to list if needed
                    if hasattr(vector, 'tolist'):
                        vector = vector.tolist()
                    elif not isinstance(vector, (list, tuple)):
                        vector = [0, 0, 0]
                    
                    # Check if vector has 3 elements
                    if len(vector) >= 3:
                        length = (vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5
                        print(f"       {vector_name}: [{vector[0]:.4f}, {vector[1]:.4f}, {vector[2]:.4f}] (length: {length:.4f} Å)")
                    else:
                        print(f"       {vector_name}: No vector information")
            
            # Vacuum information
            if vacuum_info:
                print("     Vacuum information saved")
                
                if coord_system == 'cartesian' and 'new_cell_size' in vacuum_info:
                    new_size = vacuum_info['new_cell_size']
                    print(f"       New cell size: [{new_size[0]:.4f}, {new_size[1]:.4f}, {new_size[2]:.4f}] Å")
                elif coord_system == 'user_defined' and 'new_vectors' in vacuum_info:
                    print("       User-defined vectors updated")
            
        except Exception as e:
            print(f"ERROR: Error displaying vacuum addition results: {str(e)}")
    
    def display_qe_output_summary(self, qe_result):
        """Display Quantum ESPRESSO output result summary"""
        try:
            if not qe_result:
                return
            
            if qe_result.get('skipped'):
                print("Quantum ESPRESSO output file generation skipped.")
                return
            
            files = qe_result.get('files', {})
            
            # Extract filename from file path
            for file_type, filepath in files.items():
                if file_type == 'full_input':
                    print(f"\nQuantum ESPRESSO Output File Summary: {filepath}")
                    break
            
        except Exception as e:
            print(f"ERROR: Error displaying QE output result: {str(e)}")
    
    def display_phase4_summary(self, result_data):
        """Display Phase4 overall result summary"""
        try:
            print("\nPhase4 Completion Summary:")
            print("=" * 50)
            
            # Vacuum information
            vacuum_config = result_data.get('vacuum_config', {})
            if vacuum_config:
                coord_system = vacuum_config.get('coordinate_system', 'unknown')
                vacuum_layers = vacuum_config.get('vacuum_layers', {})
                total_vacuum = sum(vacuum_layers.values())
                
                print(f"Vacuum layer added: {coord_system} coordinate system, total {total_vacuum:.2f} Å")
            
            # QE output information
            qe_output = result_data.get('qe_output', {})
            if qe_output and not qe_output.get('skipped'):
                files_count = len(qe_output.get('files', {}))
                print(f"Quantum ESPRESSO output: {files_count} files generated")
            elif qe_output and qe_output.get('skipped'):
                print("Quantum ESPRESSO output: Skipped")
            
            # Phase5 progression
            proceed_to_phase5 = result_data.get('proceed_to_phase5', False)
            if proceed_to_phase5:
                print("Next step: Proceeding to Phase5 hydrogenation process")
            else:
                print("Program end: Hydrogenation process skipped")
            
            # Final structure information
            final_structure = result_data.get('final_structure', {})
            if final_structure:
                atoms = final_structure.get('atoms', [])
                print(f"Final structure: {len(atoms)} atoms")
            
        except Exception as e:
            print(f"ERROR: Error displaying Phase4 summary: {str(e)}") 