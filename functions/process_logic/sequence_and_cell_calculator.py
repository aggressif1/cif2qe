#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase3 sequence analysis and cell vector calculator
Performs D-value based sequence analysis and user-defined cell vector calculations.
"""

class SequenceAndCellCalculator:
    """Phase3 sequence analysis and cell vector calculator"""
    
    def __init__(self):
        """Initialize Phase3SequenceAndCellCalculator"""
        from .sequence_analyzer import SequenceAnalyzer
        from .cell_vector_calculator import CellVectorCalculator
        from .unit_cell_creator import UnitCellCreator
        
        self.sequence_analyzer = SequenceAnalyzer()
        self.cell_vector_calculator = CellVectorCalculator()
        self.unit_cell_creator = UnitCellCreator()
    
    def perform_sequence_analysis_and_cell_vector_calculation(self):
        """Phase 3-4: D-value based sequence analysis and user-defined cell vector calculation"""
        print("Starting Phase 3-4: Sequence analysis and cell vector calculation")
        
        try:
            # Sequence analysis
            sequence_info = self.sequence_analyzer.analyze_d_value_sequence()
            if not sequence_info:
                print("ERROR: Sequence analysis failed")
                return None
            
            # Cell vector calculation (plane-based without transformation matrix)
            user_cell_vectors = self.cell_vector_calculator.initialize_user_cell_vectors()
            
            # In-plane vector calculation (without transformation matrix)
            a_user, b_user = self.cell_vector_calculator.calculate_in_plane_vectors(None)
            if a_user is None or b_user is None:
                print("ERROR: In-plane vector calculation failed")
                return None
            
            # Out-of-plane vector calculation (based on sequence information)
            c_user = self.cell_vector_calculator.calculate_out_of_plane_vector(sequence_info, None)
            if c_user is None:
                print("ERROR: Out-of-plane vector calculation failed")
                return None
            
            self.cell_vector_calculator.display_final_cell_vectors(a_user, b_user, c_user, sequence_info, None)
            
            # Create unit cell group
            unit_cell_group = self.unit_cell_creator.create_new_unit_cell_group(a_user, b_user, c_user)
            if unit_cell_group:
                print("SUCCESS: Phase 3-5: Unit cell creation completed")
                return unit_cell_group
            else:
                print("WARNING: Phase 3-5: Unit cell creation failed")
                return None
            
        except Exception as e:
            print(f"ERROR: Phase 3-4 execution failed: {str(e)}")
            return None 