from ..input_logic.supercell_input import SupercellInput
from ..process_logic.supercell_generator import SupercellGenerator
from .writers.supercell_writer import SupercellWriter
import os

class SupercellProcessor:
    """Module responsible for supercell generation"""
    
    def __init__(self):
        """Initialize SupercellProcessor"""
        self.input_handler = SupercellInput()
        self.generator = SupercellGenerator()
        self.writer = SupercellWriter()
    
    def generate_supercell(self, atoms, cell_vectors, nx, ny, nz, output_dir='.'):
        """
        Generate supercell.
        
        Args:
            atoms (list): List of atomic information
            cell_vectors (dict): Cell vector information
            nx, ny, nz (int): Supercell size in each direction
            output_dir (str): Output directory path
            
        Returns:
            tuple: (supercell atomic information, supercell vectors)
        """
        print("\n" + "=" * 60)
        print(f"Supercell Information ({nx}x{ny}x{nz})")
        print("=" * 60)
        
        # Calculate supercell vectors
        supercell_vectors = self.generator.calculate_supercell_vectors(cell_vectors, nx, ny, nz)
        
        # Generate supercell atoms
        supercell_atoms = self.generator.generate_supercell_atoms(atoms, nx, ny, nz, cell_vectors)
        
        # DISABLED: Move atoms to center of cell - causes fractional coordinates > 1 in non-orthogonal lattices
        # supercell_atoms = self.generator.center_atoms(supercell_atoms, supercell_vectors)
        
        # Create QE output directory
        qe_output_dir = output_dir
        os.makedirs(qe_output_dir, exist_ok=True)
        
        # Display supercell information
        self.writer.display_supercell_info(
            atoms, supercell_atoms, cell_vectors, supercell_vectors
        )
        
        return supercell_atoms, supercell_vectors 