#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase5: Hydrogenation process manager
- Add hydrogen atoms to surface to saturate dangling bonds
- Optimize hydrogen positions and prevent collisions
- Final Quantum ESPRESSO output
"""

import os
import sys
from typing import Dict, Optional

# Add project root directory to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

# Import from functions package (centralized imports)
from functions import (
    ComponentFactory
)

class Phase5Hydrogenation:
    """Phase5: Hydrogenation process"""
    
    def __init__(self, factory=None):
        """
        Initialize Phase5Hydrogenation with factory pattern.

        Args:
            factory: ComponentFactory instance (optional, uses default if None)
        """
        # Use factory pattern for component creation
        self.factory = factory or ComponentFactory()
        components = self.factory.create_all_phase5_components()

        # Initialize components from factory
        self.hydrogenation_input = components['hydrogenation_input']
        self.hydrogen_placer = components['hydrogen_placer']
        self.hydrogen_optimizer = components['hydrogen_optimizer']
        self.qe_output_generator = components['qe_output_generator']
        self.display_manager = components['display_manager']
        self.temp_handler = components['temp_handler']
        
    def execute(self, phase4_data):
        """Execute Phase5: Hydrogenation process"""
        print("PHASE5: Hydrogenation Process Start")
        
        try:
            # Validate input data
            if not self.hydrogen_optimizer.validate_input_data(phase4_data):
                return None
            
            # Extract Phase4 data (based on slab structure)
            structure_type = phase4_data.get('structure_type')
            if structure_type == 'slab':
                final_structure = phase4_data.get('slab_structure')
                vacuum_config = phase4_data.get('vacuum_sizes')
                miller_indices = phase4_data.get('phase3_data', {}).get('miller_indices')
            else:
                print("ERROR: Phase5 can only be executed with slab structures.")
                return None
            
            # Display input information
            self.display_manager.display_input_info(phase4_data)
            
            # Input hydrogenation settings
            hydrogen_config = self.hydrogenation_input.get_hydrogenation_configuration()
            if not hydrogen_config:
                print("ERROR: Hydrogenation configuration input failed")
                return None
            
            # Define hydrogen placement region
            print("\nDefining hydrogen placement region...")
            placement_region = self.hydrogen_placer.define_placement_region(
                final_structure, hydrogen_config
            )
            
            if not placement_region:
                print("ERROR: Hydrogen placement region definition failed")
                return None
            
            # Place hydrogen atoms
            print("\nPlacing hydrogen atoms...")
            hydrogenated_structure = self.hydrogen_placer.place_hydrogen_atoms(
                final_structure, placement_region, hydrogen_config
            )
            
            if not hydrogenated_structure:
                print("ERROR: Hydrogen atom placement failed")
                return None
            
            # Check collisions and optimize positions
            print("\nOptimizing hydrogen positions...")
            optimized_structure = self.hydrogen_optimizer.optimize_hydrogen_positions(hydrogenated_structure)
            
            # Display results
            self.display_manager.display_hydrogenation_results(optimized_structure, hydrogen_config)
            
            # Extract required information from Phase1 data
            source_folder = None
            cif_filename = None
            if phase4_data.get('phase3_data', {}).get('phase1_data'):
                source_folder = phase4_data['phase3_data']['phase1_data'].get('source_folder')
                cif_filename = phase4_data['phase3_data']['phase1_data'].get('filename')
            
            # Load file information directly from Phase 1 temporary data
            phase1_data = self.temp_handler.load_phase_results("phase1")
            
            if phase1_data:
                source_folder = phase1_data.get('source_folder')
                cif_filename = phase1_data.get('filename')
                print(f"   File information restored from Phase 1: {cif_filename} (folder: {source_folder})")
            else:
                print("   WARNING: Phase 1 data not found.")
            
            # Extract supercell size information
            supercell_size = None
            if phase4_data.get('phase3_data', {}).get('expanded_supercell'):
                supercell_size = phase4_data['phase3_data']['expanded_supercell'].get('supercell_size', {})
            
            # Generate final QE output file
            final_qe_output = self.qe_output_generator.generate_final_qe_output(optimized_structure, miller_indices, source_folder, cif_filename, supercell_size)
            
            # Compose result data
            result_data = {
                'original_structure': final_structure,
                'hydrogen_config': hydrogen_config,
                'placement_region': placement_region,
                'hydrogenated_structure': optimized_structure,
                'final_qe_output': final_qe_output,
                'miller_indices': miller_indices,
                'phase4_data': phase4_data
            }
            
            # Display Phase5 completion summary
            self.display_manager.display_phase5_summary(result_data)
            
            print("PHASE5: Hydrogenation Process Complete")
            print("All Phases Complete! Terminating program.")
            return result_data
            
        except Exception as e:
            print(f"ERROR: Error during Phase5 execution: {str(e)}")
            return None


def main():
    """Main function used when running Phase5 independently"""
    print("WARNING: Phase5 requires results from Phase4.")
    print("Please run the full program (start.py).")
    return None


if __name__ == "__main__":
    main() 