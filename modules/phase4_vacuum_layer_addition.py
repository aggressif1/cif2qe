#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase4: Vacuum layer addition and Quantum ESPRESSO output manager
- Add vacuum layers around newly generated supercell
- Generate output files in Quantum ESPRESSO format
- Connect to Phase5 hydrogenation process
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

class Phase4VacuumLayerAddition:
    """Phase4: Vacuum layer addition and Quantum ESPRESSO output"""
    
    def __init__(self, factory=None):
        """
        Initialize Phase4VacuumLayerAddition with factory pattern.

        Args:
            factory: ComponentFactory instance (optional, uses default if None)
        """
        # Use factory pattern for component creation
        self.factory = factory or ComponentFactory()
        components = self.factory.create_all_phase4_components()

        # Initialize components from factory
        self.vacuum_input = components['vacuum_input']
        self.vacuum_processor = components['vacuum_processor']
        self.coordinate_transformer = components['coordinate_transformer']
        self.structure_processor = components['structure_processor']
        self.qe_writer = components['qe_writer']
        self.display_manager = components['display_manager']
        
    def execute(self, phase3_data):
        """Execute Phase4: Vacuum layer addition and QE output"""
        print("PHASE4: Vacuum Layer Addition and Quantum ESPRESSO Output Start")
        
        try:
            # Validate input data
            if not self.structure_processor.validate_input_data(phase3_data):
                return None
            
            # Extract Phase3 data
            unit_cell_group = phase3_data.get('unit_cell_group')
            expanded_supercell = phase3_data.get('expanded_supercell')
            miller_indices = phase3_data.get('miller_indices')
            
            # Display input information
            self.display_manager.display_input_info(phase3_data)
            
            # Select structure type (bulk vs slab)
            structure_type = self.structure_processor.ask_structure_type()
            if not structure_type:
                print("ERROR: Structure type selection was cancelled.")
                return None
            
            print(f"Selected structure type: {structure_type}")
            
            # Generate bulk structure
            if structure_type == 'bulk':
                return self.structure_processor.process_bulk_structure(phase3_data, expanded_supercell, miller_indices)
            
            # Generate slab structure
            elif structure_type == 'slab':
                # Slab structure processing is handled by separate module (future implementation)
                from functions.process_logic.slab_processor import SlabProcessor
                slab_processor = SlabProcessor()
                return slab_processor.process_slab_structure(phase3_data, expanded_supercell, miller_indices)
            
            else:
                print("ERROR: Unknown structure type.")
                return None
                
        except Exception as e:
            print(f"ERROR: Error during Phase4 execution: {str(e)}")
            return None


def main():
    """Main function used when running Phase4 independently"""
    print("WARNING: Phase4 requires results from Phase3.")
    print("Please run the full program (start.py).")
    return None


if __name__ == "__main__":
    main() 