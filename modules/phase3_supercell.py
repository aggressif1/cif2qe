#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase3: Miller index-based supercell generation and plane analysis manager
- Input Miller indices from user
- Generate plane equations
- Assign atoms to planes
- Utilize existing functions from functions folder
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


class Phase3MillerIndexSupercell:
    """Phase3: Miller index-based supercell generation and plane analysis"""
    
    def __init__(self, factory=None):
        """
        Initialize Phase3MillerIndexSupercell with factory pattern.

        Args:
            factory: ComponentFactory instance (optional, uses default if None)
        """
        # Use factory pattern for component creation
        self.factory = factory or ComponentFactory()
        components = self.factory.create_all_phase3_components()

        # Initialize components from factory
        self.miller_input = components['miller_input']
        self.plane_generator = components['plane_generator']
        self.atom_assigner = components['atom_assigner']
        self.reference_selector = components['reference_selector']
        self.comparison_analyzer = components['comparison_analyzer']
        self.plane_writer = components['plane_writer']
        self.data_validator = components['data_validator']
        self.reference_plane_manager = components['reference_plane_manager']
        self.supercell_expander = components['supercell_expander']
        self.miller_verification_manager = components['miller_verification_manager']
        self.sequence_and_cell_calculator = components['sequence_and_cell_calculator']
        self.qe_output_generator = components['qe_output_generator']
        self.display_manager = components['display_manager']
        self.temp_handler = components['temp_handler']

        # Initialize instance variables
        self.miller_indices = None
        self.reference_plane_info = None
        self.unit_cell_group = None
    
    def execute(self, step2_data):
        """Execute Phase3: Miller index-based plane analysis"""
        self.display_manager.show_phase3_start()
        
        try:
            # Validate and extract input data
            if not self.data_validator.validate_phase2_data(step2_data):
                return None
            
            supercell_atoms = step2_data['supercell_atoms']
            cell_vectors = step2_data['cell_vectors']
            supercell_vectors = step2_data['supercell_vectors']
            file_info = step2_data.get('file_info', {})
            
            self.cell_vectors = cell_vectors
            self.all_supercell_atoms = supercell_atoms
            
            # Display input information
            self.display_manager.show_input_data_info(
                filename=file_info.get('filename', 'unknown.cif'),
                atom_count=len(supercell_atoms)
            )
            
            # Calculate element counts
            element_counts = {}
            for atom in supercell_atoms:
                element = atom.get('element', atom.get('atom_type', 'Unknown'))
                element_counts[element] = element_counts.get(element, 0) + 1
            
            self.display_manager.show_element_counts(element_counts)
            
            # Input Miller indices (with crystal system info)
            # Get crystal system from Phase1 data
            phase1_data = self.temp_handler.load_phase_results("phase1")
            crystal_system = None
            if phase1_data:
                crystal_system = phase1_data.get('crystal_system', 'Unknown')
                print(f"   - Crystal system: {crystal_system}")
            
            miller_indices = self.miller_input.get_miller_indices(crystal_system)
            if not miller_indices:
                self.display_manager.show_error_during_phase3("Miller indices input failed")
                return None
            self.miller_indices = miller_indices
            
            # Save Miller indices to temporary file
            phase3_temp_data = {'miller_indices': miller_indices}
            self.temp_handler.save_phase_results("phase3", phase3_temp_data)

            
            # Generate plane equations
            plane_equation = self.plane_generator.generate_plane_equation(miller_indices, cell_vectors)
            if not plane_equation:
                self.display_manager.show_error_during_phase3("Plane equation generation failed")
                return None
            
            # Assign atoms to planes
            plane_assignments = self.atom_assigner.assign_atoms_to_planes(supercell_atoms, plane_equation)
            if not plane_assignments:
                print("ERROR: Atom plane assignment failed")
                return None
            
            # Display results and update CSV
            self.plane_writer.display_plane_analysis_results(miller_indices, plane_equation, plane_assignments)
            
            csv_file_path = "output/supercell_coordinates.csv"
            self.reference_selector.update_csv_with_all_plane_info(csv_file_path, plane_assignments)
            
            # Select reference plane
            reference_plane = self.reference_selector.select_reference_plane(plane_assignments)
            if not reference_plane:
                print("WARNING: Reference plane selection failed")
                return None
            
            self.reference_plane_info = reference_plane
            self.reference_selector.display_reference_plane_summary()
            self.reference_selector.update_csv_with_reference_plane(csv_file_path, plane_assignments)
            
            # Setup bond analysis
            phase1_data = self.temp_handler.load_phase_results("phase1")
            
            if phase1_data and 'bond_info' in phase1_data:
                bond_info = phase1_data['bond_info']
                bond_criteria = {
                    'bond_type': bond_info.get('bond_type', 'Unknown'),
                    'min_distance': bond_info.get('min_distance', 0.5),
                    'max_distance': bond_info.get('max_distance', 3.0)
                }
            else:
                from functions.output_logic.console_writer import ConsoleWriter
                from functions.process_logic.bond_calculator import BondCalculator
                
                bond_calc = BondCalculator()
                min_distance, max_distance, bond_type = ConsoleWriter.show_distance_analysis_and_get_bond_range(
                    bond_calc, supercell_atoms, cell_vectors
                )
                bond_criteria = {
                    'bond_type': bond_type if bond_type else 'Unknown',
                    'min_distance': min_distance,
                    'max_distance': max_distance
                }
            
            self.comparison_analyzer.set_bond_criteria(bond_criteria)
            self.comparison_analyzer.set_cell_vectors(supercell_vectors)
            
            # Iterative reference plane analysis (starting from 1st reference plane)
            self.reference_plane_manager.perform_iterative_reference_plane_analysis(
                plane_assignments, reference_plane, supercell_atoms, cell_vectors, miller_indices,
                self.comparison_analyzer
            )
            
            # Sequence analysis and cell vector calculation (unit cell generation)
            unit_cell_group = self.sequence_and_cell_calculator.perform_sequence_analysis_and_cell_vector_calculation()
            if unit_cell_group:
                self.unit_cell_group = unit_cell_group
            
            # Prepare result data
            result_data = {
                'miller_indices': miller_indices,
                'plane_equation': plane_equation,
                'plane_assignments': plane_assignments,
                'reference_plane': reference_plane,
                'supercell_atoms': supercell_atoms,
                'cell_vectors': cell_vectors,
                'file_info': file_info,
                'original_atoms': step2_data.get('original_atoms')
            }
            
            # Process unit cell group and Miller index-based verification
            if self.unit_cell_group:
                result_data['unit_cell_group'] = self.unit_cell_group
                
                # Generate supercell
                supercell_size = self.supercell_expander.get_supercell_size_input()
                if supercell_size:
                    expanded_supercell = self.supercell_expander.generate_expanded_supercell(
                        self.unit_cell_group, supercell_size
                    )
                    if expanded_supercell:
                        result_data['expanded_supercell'] = expanded_supercell
                        
                        # Generate QE files for supercell
                        qe_output_result = self.qe_output_generator.generate_supercell_qe_output(
                            expanded_supercell, miller_indices, file_info
                        )
                        if qe_output_result:
                            result_data['supercell_qe_output'] = qe_output_result
                        
                        # Perform Miller index-based unit cell vs supercell verification
                        self.miller_verification_manager.perform_miller_verification(result_data)
                        
                        # Display QE file generation result summary
                        self.qe_output_generator.display_qe_output_summary(result_data)
            
            print("PHASE3 Complete")
            return result_data
            
        except Exception as e:
            print(f"ERROR: Error during Phase3 execution: {str(e)}")
            return None


def main():
    """Main function used when running Phase3 independently"""
    print("WARNING: Phase3 requires results from Phase2.")
    print("Please run the full program (start.py).")
    return None


if __name__ == "__main__":
    main() 