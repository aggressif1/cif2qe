from ..input_logic.file_reader import FileReader
from ..input_logic.cif_reader import CIFReader
from ..input_logic.data_parser import CIFParser
from ..process_logic.calculator import CellCalculator
from ..output_logic.console_writer import ConsoleWriter
from ..utils.coordinates import convert_atoms_to_cartesian
import os

class CIFProcessor:
    """Class that integrates all components for CIF file processing"""
    
    def __init__(self):
        """Initialize components"""
        self.file_reader = FileReader()
        self.cif_reader = CIFReader()
        self.parser = CIFParser()
        self.calculator = CellCalculator()
        self.writer = ConsoleWriter()
    
    def find_cif_files(self, base_dir='.'):
        """
        Find CIF files in directory.
        
        Args:
            base_dir (str): Base directory to start search
            
        Returns:
            dict: Dictionary in {folder_name: [CIF file paths]} format
        """
        return self.cif_reader.find_cif_files(base_dir)
    
    def display_file_menu(self, cif_files):
        """
        Display discovered CIF file list and get user selection.

        Args:
            cif_files (dict): Dictionary of CIF file paths by folder

        Returns:
            tuple: (selected folder, selected file path) or list of file paths for multiple selection
        """
        # Display CIF files using ConsoleWriter (no Selection Guide)
        self.writer.display_cif_files(cif_files)

        # Create a flat list of all files with their numbers
        all_files = []
        file_counter = 1

        for folder, paths in cif_files.items():
            for path in paths:
                all_files.append({
                    'number': file_counter,
                    'folder': folder,
                    'path': path,
                    'filename': os.path.basename(path)
                })
                file_counter += 1

        # Display Selection Guide
        print("\n" + "=" * 80)
        print("Selection Guide:")
        print("• \"3\" → Process single file")
        print("• \"3,4\" → Process multiple files")
        print("• \"0\" → Process all files")
        print("• \"x\" → Exit program")
        print("=" * 80)

        while True:
            choice = input("\nEnter your choice: ").strip().lower()

            if choice == 'x':
                return 'exit', None
            elif choice == '0':
                return 'all', None

            # Handle multiple file selection (e.g., "3,4,7")
            elif ',' in choice:
                try:
                    selected_numbers = [int(x.strip()) for x in choice.split(',')]
                    selected_files = []

                    for num in selected_numbers:
                        if 1 <= num <= len(all_files):
                            selected_files.append(all_files[num - 1])
                        else:
                            print(f"Invalid file number: {num}")
                            continue

                    if selected_files:
                        return 'multiple', selected_files
                    else:
                        print("No valid files selected. Please try again.")
                        continue

                except ValueError:
                    print("Invalid format. Use comma-separated numbers (e.g., 3,4,7)")
                    continue

            # Handle single file selection
            else:
                try:
                    file_num = int(choice)
                    if 1 <= file_num <= len(all_files):
                        selected_file = all_files[file_num - 1]
                        return selected_file['folder'], selected_file['path']
                    else:
                        print(f"Invalid file number. Please enter a number between 1 and {len(all_files)}")
                        continue
                except ValueError:
                    print("Invalid input. Please enter a number, comma-separated numbers, '0' for all, or 'x' to exit")
                    continue
    
    def analyze_cif(self, cif_path):
        """
        Analyze CIF file.
        
        Args:
            cif_path (str): CIF file path
            
        Returns:
            dict: Dictionary containing analyzed data
                 {'cell_params': {...}, 'cell_vectors': {...}, 'atoms': [...]}
        """
        # Read file
        content = self.file_reader.read_file(cif_path)
        if not content:
            print(f"Error: Could not read CIF file: {cif_path}")
            return None
        
        # Extract unit cell parameters
        cell_params = self.parser.parse_cell_parameters(content)
        if not cell_params:
            print("Error: Could not extract cell parameters from CIF file.")
            return None
        
        # Calculate cell vectors
        cell_vectors = self.calculator.calculate_cell_vectors(cell_params)
        if not cell_vectors:
            print("Error: Could not calculate cell vectors.")
            return None
        
        # Extract atomic positions
        base_atoms = self.parser.parse_atomic_positions(content)
        
        if base_atoms is None:
            print("\nError parsing atomic positions from CIF file.")
            return None
        elif len(base_atoms) == 0:
            print("\nNo atomic positions found in the CIF file.")
            return None
        
        print(f"Found {len(base_atoms)} base atomic positions")
        
        # Extract symmetry operations
        symmetry_ops = self.parser.parse_symmetry_operations(content)
        print(f"Found {len(symmetry_ops)} symmetry operations")
        if len(symmetry_ops) <= 5:  # Print only first few
            for i, op in enumerate(symmetry_ops):
                print(f"  {i+1}: {op}")
        
        # Generate all atoms by applying symmetry operations
        atoms = self.parser.apply_symmetry_operations(base_atoms, symmetry_ops)
        print(f"Generated {len(atoms)} atoms from {len(base_atoms)} unique positions")
        
        # Calculate Cartesian coordinates
        atoms = convert_atoms_to_cartesian(atoms, cell_vectors)
        
        # Output results
        cell_info = {
            'params': cell_params,
            'vectors': cell_vectors,
            'volume': self.calculator.calculate_volume(cell_params),
            'system': self.calculator.determine_crystal_system(cell_params)
        }
        self.writer.write_cell_info(cell_info)
        bond_info = self.writer.write_atomic_positions(atoms, cell_vectors)
        
        return {
            'cell_params': cell_params,
            'cell_vectors': cell_vectors,
            'atoms': atoms,
            'volume': cell_info['volume'],
            'crystal_system': cell_info['system'],
            'bond_info': bond_info
        } 