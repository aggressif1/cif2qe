"""
Phase5 output management related functions
"""
from ..utils.element_counter import ElementCounter
from ..utils.miller_display_helper import MillerDisplayHelper
from .display_manager_base import BaseDisplayManager

class Phase5DisplayManager(BaseDisplayManager):
    """Class responsible for Phase5 output management"""
    
    def __init__(self):
        """Initialize Phase5DisplayManager"""
        pass
    
    def display_input_info(self, phase4_data):
        """Display input information received from Phase4"""
        try:
            print("\nPhase5 Input Data:")
            print("-" * 40)
            
            # Check structure type
            structure_type = phase4_data.get('structure_type', 'unknown')
            print(f"   - Structure type: {structure_type}")
            
            if structure_type == 'slab':
                # Vacuum layer information (slab structure)
                vacuum_sizes = phase4_data.get('vacuum_sizes', {})
                if vacuum_sizes:
                    # Display a, b, c direction vacuum layers in 2 lines
                    a_info = " | ".join([f"{k}: {v:.2f} Å" for k, v in vacuum_sizes.items() if k.startswith('a_')])
                    bc_info = " | ".join([f"{k}: {v:.2f} Å" for k, v in vacuum_sizes.items() if k.startswith(('b_', 'c_'))])
                    print("   - Vacuum layer configuration:")
                    if a_info:
                        print(f"     a-direction: {a_info}")
                    if bc_info:
                        print(f"     b,c-directions: {bc_info}")
                
                # Slab structure information
                slab_structure = phase4_data.get('slab_structure', {})
                if slab_structure:
                    atoms = slab_structure.get('atoms', [])
                    print(f"   - Atoms with vacuum layer added: {len(atoms)}")
                    
                    # Element count
                    element_counts = ElementCounter.count_elements(atoms)
                    
                    if element_counts:
                        print("   - Atoms by element:")
                        for element, count in sorted(element_counts.items()):
                            print(f"     {element}: {count}")
            
            else:
                # Legacy method (backward compatibility)
                vacuum_config = phase4_data.get('vacuum_config', {})
                if vacuum_config:
                    coord_system = vacuum_config.get('coordinate_system', 'unknown')
                    vacuum_layers = vacuum_config.get('vacuum_layers', {})
                    total_vacuum = sum(vacuum_layers.values())
                    
                    print(f"   - Vacuum layer: {coord_system} coordinate system, total {total_vacuum:.2f} Å")
                
                # Final structure information
                final_structure = phase4_data.get('final_structure', {})
                if final_structure:
                    atoms = final_structure.get('atoms', [])
                    print(f"   - Atoms with vacuum layer added: {len(atoms)}")
                    
                    # Element count
                    element_counts = ElementCounter.count_elements(atoms)
                    
                    if element_counts:
                        print("   - Atoms by element:")
                        for element, count in sorted(element_counts.items()):
                            print(f"     {element}: {count}")
            
            # Miller indices
            MillerDisplayHelper.display_miller_indices(phase4_data)
            
        except Exception as e:
            print(f"ERROR: Error displaying input information: {str(e)}")
    
    def display_hydrogenation_results(self, structure, hydrogen_config):
        """Display hydrogenation results"""
        try:
            print("\nHydrogenation Results:")
            print("-" * 40)
            
            # Hydrogenation configuration summary
            direction_config = hydrogen_config.get('direction_config', {})
            hydrogen_count = hydrogen_config.get('hydrogen_count', 0)
            
            # Check active directions
            active_directions = []
            for vector_name, directions in direction_config.items():
                if directions['positive']['enabled']:
                    active_directions.append(f"{vector_name}(+)")
                if directions['negative']['enabled']:
                    active_directions.append(f"{vector_name}(-)")
            
            if active_directions:
                print(f"   Hydrogenation directions: {', '.join(active_directions)}")
            else:
                print("   Hydrogenation directions: None")
            
            print(f"   Requested hydrogen count: {hydrogen_count}")
            
            # Actual hydrogenation results
            hydrogenation_info = structure.get('hydrogenation_info', {})
            if hydrogenation_info:
                actual_h_count = hydrogenation_info.get('hydrogen_count', 0)
                original_count = hydrogenation_info.get('original_atom_count', 0)
                total_count = hydrogenation_info.get('total_atom_count', 0)
                
                print(f"   Actually added hydrogen: {actual_h_count}")
                print(f"   Atom count change: {original_count} → {total_count}")
            
            # Final element count
            atoms = structure.get('atoms', [])
            element_counts = ElementCounter.count_elements(atoms)
            
            print("   Final atoms by element:")
            for element, count in sorted(element_counts.items()):
                print(f"     {element}: {count}")
            
        except Exception as e:
            print(f"ERROR: Error displaying hydrogenation results: {str(e)}")
    
    def display_final_qe_summary(self, qe_result):
        """Display final QE output result summary"""
        try:
            if not qe_result:
                print("ERROR: No final QE output result.")
                return
            
            print("\nFinal Quantum ESPRESSO Output Files:")
            print("-" * 50)
            
            files = qe_result.get('files', {})
            base_filename = qe_result.get('base_filename', 'unknown')
            
            print(f"   Base filename: {base_filename}")
            print("   Generated files:")
            
            for file_type, filepath in files.items():
                filename = filepath.split('/')[-1] if '/' in filepath else filepath
                if file_type == 'full_input':
                    print(f"     Final QE input file: {filename}")
                else:
                    print(f"     {file_type}: {filename}")
            
            print("\n   Final structure features:")
            print("     - Surface structure with vacuum layer added")
            print("     - Dangling bonds saturated with hydrogen atoms")
            print("     - Ready for Quantum ESPRESSO calculation")
            
        except Exception as e:
            print(f"ERROR: Error displaying final QE output result: {str(e)}")
    
    def display_phase5_summary(self, result_data):
        """Display Phase5 overall result summary"""
        try:
            print("\nPhase5 Completion Summary:")
            print("=" * 50)
            
            # Hydrogenation information
            hydrogen_config = result_data.get('hydrogen_config', {})
            if hydrogen_config:
                direction_config = hydrogen_config.get('direction_config', {})
                hydrogen_count = hydrogen_config.get('hydrogen_count', 0)
                
                # Check active directions
                active_directions = []
                for vector_name, directions in direction_config.items():
                    if directions['positive']['enabled']:
                        pos_config = directions['positive']
                        active_directions.append(f"{vector_name}(+): {pos_config['min_range']:.1f}~{pos_config['max_range']:.1f}Å")
                    if directions['negative']['enabled']:
                        neg_config = directions['negative']
                        active_directions.append(f"{vector_name}(-): {neg_config['min_range']:.1f}~{neg_config['max_range']:.1f}Å")
                
                if active_directions:
                    print(f"Hydrogenation: {len(active_directions)} directions, {hydrogen_count} hydrogen atoms added")
                    for direction in active_directions:
                        print(f"   - {direction}")
                else:
                    print(f"Hydrogenation: Disabled")
            
            # Final structure information
            hydrogenated_structure = result_data.get('hydrogenated_structure', {})
            if hydrogenated_structure:
                hydrogenation_info = hydrogenated_structure.get('hydrogenation_info', {})
                if hydrogenation_info:
                    actual_h = hydrogenation_info.get('hydrogen_count', 0)
                    total_atoms = hydrogenation_info.get('total_atom_count', 0)
                    print(f"Final structure: {total_atoms} atoms (including {actual_h} hydrogen atoms)")
            
            # QE output information
            final_qe_output = result_data.get('final_qe_output', {})
            if final_qe_output and final_qe_output.get('success'):
                files_count = len(final_qe_output.get('files', {}))
                print(f"Final QE output: {files_count} files generated")
            
            print("\nCIF2QE program completed successfully!")
            print("   Check the generated files in the output/ folder")
            print("   Ready to start Quantum ESPRESSO calculation")
            
        except Exception as e:
            print(f"ERROR: Error displaying Phase5 summary: {str(e)}") 