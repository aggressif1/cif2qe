import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import from functions package (centralized imports)
from functions import (
    ComponentFactory,
    Phase2DisplayManager
)

class Phase2SupercellForPatternDetection:
    """
    Phase2: Supercell generation for pattern detection
    
    Receives CIF data analyzed in Step1 and generates 6x6x6 supercell.
    """
    
    def __init__(self, factory=None):
        """
        Initialize Phase2SupercellForPatternDetection with factory pattern.

        Args:
            factory: ComponentFactory instance (optional, uses default if None)
        """
        # Use factory pattern for component creation
        self.factory = factory or ComponentFactory()
        components = self.factory.create_all_phase2_components()

        # Initialize components from factory
        self.supercell_processor = components['supercell_processor']
        self.supercell_verifier = components['supercell_verifier']
        self.data_validator = components['data_validator']
        self.supercell_writer = components['supercell_writer']
        self.display_manager = components['display_manager']
        self.temp_handler = components['temp_handler']
        
    def execute(self, cif_data, output_dir='OUTPUT'):
        """
        Execute Step2 supercell generation.
        
        Args:
            cif_data (dict): CIF data analyzed in Step1
            output_dir (str): Output directory
            
        Returns:
            dict: Supercell generation result data
        """
        self.display_manager.show_phase2_start()
        
        # Validate Phase1 data (using DataValidator from functions)
        if not self.data_validator.validate_phase1_data(cif_data):
            self.display_manager.show_phase1_data_error()
            return None
            
        # Extract required information from CIF data
        atoms = cif_data['atoms']
        cell_vectors = cif_data['cell_vectors']
        filename = cif_data.get('filename', 'unknown.cif')
        
        self.display_manager.show_input_data_info(
            filename=filename,
            atom_count=len(atoms),
            volume=cif_data.get('volume', 'N/A')
        )
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # Generate 6x6x6 supercell (hard-coded)
            nx, ny, nz = 6, 6, 6
            self.display_manager.show_supercell_generation_info(nx, ny, nz)
            supercell_atoms, supercell_vectors = self.supercell_processor.generate_supercell(
                atoms, cell_vectors, nx, ny, nz, output_dir
            )
            
            # Display supercell statistics (using SupercellWriter from functions)
            self.supercell_writer.display_detailed_supercell_statistics(
                atoms, supercell_atoms, cif_data, (nx, ny, nz)
            )
            
            # Verify supercell structure
            self.display_manager.show_verification_start()
            verification_result = self.supercell_verifier.verify_supercell_structure(
                atoms, supercell_atoms, cell_vectors, supercell_vectors
            )

            if not verification_result:
                self.display_manager.show_verification_warning()
            
            # Save supercell coordinates as CSV file
            if verification_result:
                self.display_manager.show_csv_save_start()
                csv_filename = os.path.join(output_dir, "supercell_coordinates.csv")
                self.supercell_processor.writer.save_to_csv(supercell_atoms, csv_filename)
            
            # Compose result data (in format required by Step3)
            result_data = {
                'supercell_atoms': supercell_atoms,
                'supercell_vectors': supercell_vectors,
                'cell_vectors': cell_vectors,  # Required for plane equation calculation in Step3
                'file_info': {
                    'filename': filename,
                    'original_atom_count': len(atoms),
                    'supercell_atom_count': len(supercell_atoms),
                    'supercell_size': (nx, ny, nz)
                },
                'original_atoms': atoms,
                'supercell_size': (nx, ny, nz),
                'total_atoms': len(supercell_atoms),
                'output_directory': output_dir,
                'verification_passed': verification_result
            }
            
            # Display Phase2 completion message
            self.display_manager.show_phase2_completion(verification_result)
            
            # Save Phase2 results as temporary file
            self.temp_handler.save_phase_results("phase2", result_data)
            
            # Display result summary
            csv_path = os.path.join(output_dir, "supercell_coordinates.csv") if verification_result else None
            self.display_manager.show_result_summary(
                filename=filename,
                nx=nx, ny=ny, nz=nz,
                supercell_atom_count=len(supercell_atoms),
                verification_passed=verification_result,
                csv_path=csv_path
            )
            
            return result_data
            
        except Exception as e:
            self.display_manager.show_error_during_supercell_generation(e)
            return None