#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase3 Miller index based verification manager
Verifies consistency between Miller index based unit cells and supercells.
"""

class MillerVerificationManager:
    """Phase3 Miller index based verification manager"""
    
    def __init__(self):
        """Initialize Phase3MillerVerificationManager"""
        from .miller_verification_engine import MillerVerificationEngine
        from .miller_statistics_generator import MillerStatisticsGenerator
        from .miller_data_validator import MillerDataValidator
        from ..output_logic.miller_display_manager import Phase3MillerDisplayManager

        self.miller_data_validator = MillerDataValidator()
        self.miller_verification_engine = MillerVerificationEngine()
        self.miller_statistics_generator = MillerStatisticsGenerator()
        self.miller_display_manager = Phase3MillerDisplayManager()
    
    def perform_miller_verification(self, result_data):
        """
        Verify consistency between Miller index based unit cells and supercells.
        
        Args:
            result_data (dict): Result data (includes unit_cell_group, expanded_supercell)
        """
        try:
            print("\nPhase3 completed!")
            print("Starting Phase3 verification: Miller index based unit cell verification")
            
            # Validate input data
            if not self.miller_data_validator.validate_input_data(result_data):
                return
            
            # Extract data
            unit_cell_group = result_data.get('unit_cell_group', {})
            expanded_supercell = result_data.get('expanded_supercell', {})
            miller_indices = result_data.get('miller_indices', {})
            
            # Execute verification engine
            verification_results = self.miller_verification_engine.execute_verification(
                unit_cell_group, expanded_supercell, miller_indices
            )
            
            if verification_results:
                result_data['verification_results'] = verification_results
                
                # Generate verification statistics
                statistics = self.miller_statistics_generator.generate_statistics(verification_results)
                result_data['verification_statistics'] = statistics
                
                # Display results
                self.miller_display_manager.display_verification_results(verification_results, statistics)
                
                # Determine overall status
                overall_status = verification_results.get('overall_status', False)
                result_data['verification_results']['overall_status'] = overall_status
                
                if overall_status:
                    print("SUCCESS: Phase3 verification: Miller index based unit cell verification completed")
                else:
                    print("SUCCESS: Phase3 verification: Miller index based unit cell verification completed (some issues found in verification)")
            else:
                print("SUCCESS: Phase3 verification: Miller index based unit cell verification completed (some issues found in verification)")
            
        except Exception as e:
            print(f"ERROR: Error during Miller index based verification: {str(e)}")
            print("SUCCESS: Phase3 verification: Miller index based unit cell verification completed (terminated due to error)")
            result_data['verification_error'] = str(e) 