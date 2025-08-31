import os
import sys
from modules.phase1_cif_analysis import Phase1CifAnalysis
from modules.phase2_pattern_detection import Phase2SupercellForPatternDetection
from modules.phase3_supercell import Phase3MillerIndexSupercell
from modules.phase4_vacuum_layer_addition import Phase4VacuumLayerAddition
from modules.phase5_hydrogenation import Phase5Hydrogenation


class TeeOutput:
    """Class for simultaneous output to screen and file"""
    def __init__(self, file_path):
        self.terminal = sys.stdout
        self.log_file = open(file_path, 'w', encoding='utf-8')
    
    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush()  # Write to file immediately
    
    def flush(self):
        self.terminal.flush()
        self.log_file.flush()
    
    def close(self):
        if hasattr(self, 'log_file'):
            self.log_file.close()

def setup_output_logging():
    """Redirect stdout so all output is saved to both file and console simultaneously"""
    tee = TeeOutput('process.log')
    sys.stdout = tee
    return tee

def execute_phase1_to_phase5():
    """
    Execute Phase1 CIF analysis, Phase2 supercell generation, Phase3 Miller index plane analysis and verification, Phase4 vacuum layer addition, and Phase5 hydrogenation in sequence.
    
    Returns:
        Dict or None: Final result data
    """
    print("Starting CIF2QE Program")
    print("Current Mode: Phase1+Phase2+Phase3+Phase4+Phase5 execution (Complete surface structure generation)")
    
    # Execute Phase1 CIF analysis
    phase1_analyzer = Phase1CifAnalysis()
    cif_data = phase1_analyzer.execute()
    
    if not cif_data:
        print("\nERROR: Phase1 execution failed.")
        return None
    
    print("\nPhase1 Complete!")
    input("\nPress Enter to proceed to the next Phase...")
    
    # Execute Phase2 supercell generation
    print("\n" + "="*80)
    print("PHASE2: Starting Supercell Generation")
    print("="*80)
    phase2_analyzer = Phase2SupercellForPatternDetection()
    supercell_data = phase2_analyzer.execute(cif_data)
    
    if not supercell_data:
        print("\nERROR: Phase2 execution failed.")
        return None
        
    print("\nPhase2 Complete!")
    input("\nPress Enter to proceed to the next Phase...")
    
    # Execute Phase3 Miller index plane analysis and verification
    print("\n" + "="*80)
    print("PHASE3: Starting Miller Index-based Supercell Generation")
    print("="*80)
    phase3_analyzer = Phase3MillerIndexSupercell()
    phase3_data = phase3_analyzer.execute(supercell_data)
    
    # Include Phase1 data in Phase3 results
    if phase3_data:
        phase3_data['phase1_data'] = cif_data
    
    if not phase3_data:
        print("\nERROR: Phase3 execution failed.")
        return None
    
    # Output messages based on verification results
    if phase3_data.get('verification_skipped'):
        print("WARNING: Verification skipped as Miller index-based unit cell was not generated.")
    elif phase3_data.get('verification_results', {}).get('overall_status'):
        print("SUCCESS: Miller index-based unit cell verification completed successfully!")
    elif 'verification_results' in phase3_data:
        print("WARNING: Some issues were found in Miller index-based unit cell verification.")
    
    print("\nPhase3 Complete!")
    input("\nPress Enter to proceed to the next Phase...")
    
    # Execute Phase4 vacuum layer addition and QE output
    print("\n" + "="*80)
    print("PHASE4: Starting Vacuum Layer Addition and Quantum ESPRESSO Output")
    print("="*80)
    phase4_analyzer = Phase4VacuumLayerAddition()
    phase4_data = phase4_analyzer.execute(phase3_data)
    
    if not phase4_data:
        print("\nERROR: Phase4 execution failed.")
        # Return Phase3 results at least
        return phase3_data
    
    print("\nPhase4 Complete!")
    
    # Execute Phase5 only for slab structures
    structure_type = phase4_data.get('structure_type')
    if structure_type == 'slab':
        input("\nPress Enter to proceed to the next Phase...")
        
        # Execute Phase5 hydrogenation (slab structures only)
        print("\n" + "="*80)
        print("PHASE5: Starting Hydrogenation")
        print("="*80)
        phase5_analyzer = Phase5Hydrogenation()
        phase5_data = phase5_analyzer.execute(phase4_data)
        
        if not phase5_data:
            print("\nERROR: Phase5 execution failed.")
            # Return Phase4 results at least
            return phase4_data
        
        print("\nPhase5 Complete!")
        # Completion message already output from Phase5
        return phase5_data
    
    elif structure_type == 'bulk':
        print("\nINFO: Bulk structures do not require surface hydrogenation, completing at Phase4.")
        print("SUCCESS: All Phases Complete! Terminating program.")
        return phase4_data
    
    else:
        print("\nWARNING: Unknown structure type. Completing at Phase4.")
        return phase4_data

def main():
    """Main function"""
    # Setup output logging
    tee_output = setup_output_logging()
    
    try:
        print("=== CIF2QE: Crystallographic Information File to Quantum ESPRESSO ===")
        print("Version: 1.0.0 (Phase1+Phase2+Phase3+Phase4+Phase5 Complete surface structure generation)")
        
        # Execute Phase1+Phase2+Phase3+Phase4+Phase5 (Complete surface structure generation)
        result = execute_phase1_to_phase5()
        
        if result:
            print("\n=== Program Completed Successfully ===")
        else:
            print("\n=== Program Execution Failed ===")
        
        print("\n=== All output has been saved to process.log file ===")
        
    except KeyboardInterrupt:
        print("\n\nWARNING: Program interrupted by user.")
    except Exception as e:
        print(f"\nERROR: Unexpected error occurred: {str(e)}")
    finally:
        # Restore stdout and close file
        sys.stdout = tee_output.terminal
        tee_output.close()

if __name__ == "__main__":
    main()
