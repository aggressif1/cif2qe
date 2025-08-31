import os
import csv
import numpy as np
from ..reference_plane_calculator import ReferencePlaneCalculator

class SupercellWriter:
    """Class for outputting supercell information"""
    
    def __init__(self):
        self.coordinates_file = None
        self.csv_writer = None
        self.plane_count = 0
        self.atoms_data = []  # Store atomic data in memory
        self.output_dir = None
        self.filename = 'supercell_coordinates.csv'
        self.reference_plane_calc = ReferencePlaneCalculator()  # Primary reference plane calculator
    
    def initialize_csv(self, output_dir, filename='supercell_coordinates.csv'):
        """Initialize CSV file."""
        self.output_dir = output_dir
        self.filename = filename
        self.atoms_data = []  # Initialize data
        print(f"\nAtomic coordinates will be saved to {self.filename}.")
    

    
    def update_plane_info(self, d_value, atoms, plane_equation=None):
        """
        Update information for atoms belonging to a plane with specific D value.
        
        Args:
            d_value (float): D value of plane equation
            atoms (list): List of atoms belonging to the plane
            plane_equation (tuple): Plane equation coefficients (a, b, c, d)
        """
        if not atoms or not self.atoms_data:
            return
        
        self.plane_count += 1
        plane_id = f"Plane {self.plane_count}"
        
        # Generate complete plane equation
        plane_equation_str = ''
        if plane_equation:
            a, b, c, d = plane_equation
            plane_equation_str = f"{a:.4f}x + {b:.4f}y + {c:.4f}z + {d:.4f} = 0"
        
        # Create dictionary with coordinates of atoms belonging to the plane as keys
        plane_atoms_coords = {
            (float(atom['cart_x']), float(atom['cart_y']), float(atom['cart_z'])): 
            (plane_id, f"{d_value:.4f}", plane_equation_str) 
            for atom in atoms
        }
        
        # Update data in memory
        for atom_data in self.atoms_data:
            coord = (float(atom_data['Cart_X']), float(atom_data['Cart_Y']), float(atom_data['Cart_Z']))
            if coord in plane_atoms_coords:
                plane_id, projection_val, equation_str = plane_atoms_coords[coord]
                atom_data['Plane_ID'] = plane_id
                atom_data['Projection_Value'] = projection_val
                atom_data['Plane_Equation'] = equation_str
        
        # Save updated data to file
        self._write_to_file()
        
        # Display progress (every 10 planes)
        if self.plane_count % 10 == 0:
            print(f"Processing... {self.plane_count} planes processed")
    
    def _write_to_file(self):
        """Save data from memory to CSV file."""
        filepath = os.path.join(self.output_dir, self.filename)
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow([
                'Index', 'Element', 'Atom_ID', 'Cart_X', 'Cart_Y', 'Cart_Z',
                'Frac_X', 'Frac_Y', 'Frac_Z', 'Plane_ID', 'Projection_Value', 
                'Plane_Equation', 'Is_Reference_Plane'
            ])
            
            # Write data
            for atom_data in self.atoms_data:
                writer.writerow([
                    atom_data['Index'],
                    atom_data['Element'],
                    atom_data['Atom_ID'],
                    f"{atom_data['Cart_X']:.6f}",
                    f"{atom_data['Cart_Y']:.6f}",
                    f"{atom_data['Cart_Z']:.6f}",
                    f"{atom_data['Frac_X']:.6f}",
                    f"{atom_data['Frac_Y']:.6f}",
                    f"{atom_data['Frac_Z']:.6f}",
                    atom_data['Plane_ID'],
                    atom_data['Projection_Value'],
                    atom_data['Plane_Equation'],
                    atom_data['Is_Reference_Plane']
                ])
    
    def display_plane_statistics(self):
        """Display atom count statistics for each plane."""
        if not self.atoms_data:
            return
        
        # Calculate atom count per plane
        plane_counts = {}
        for atom_data in self.atoms_data:
            plane_id = atom_data.get('Plane_ID', '')
            if plane_id and plane_id != '':
                if plane_id not in plane_counts:
                    plane_counts[plane_id] = {
                        'count': 0,
                        'projection_value': atom_data.get('Projection_Value', ''),
                        'equation': atom_data.get('Plane_Equation', '')
                    }
                plane_counts[plane_id]['count'] += 1
        
        if plane_counts:
            print("\n" + "=" * 80)
            print("Atom Count Statistics by Plane")
            print("=" * 80)
            print(f"{'Plane ID':<12} {'Atom Count':<10} {'Projection':<12} {'Plane Equation'}")
            print("-" * 80)
            
            total_classified = 0
            for plane_id in sorted(plane_counts.keys(), key=lambda x: int(x.split()[1]) if x.split()[1].isdigit() else 0):
                info = plane_counts[plane_id]
                count = info['count']
                projection = info['projection_value']
                equation = info['equation']
                total_classified += count
                print(f"{plane_id:<12} {count:<10} {projection:<12} {equation}")
            
            print("-" * 80)
            print(f"{'Total':<12} {total_classified:<10}")
            
            # Unclassified atoms count
            unclassified = len(self.atoms_data) - total_classified
            if unclassified > 0:
                print(f"{'Unclassified':<12} {unclassified:<10}")
            
            print("=" * 80)
            
            # Set primary reference plane
            self.reference_plane_calc.set_reference_plane(plane_counts)
            self.reference_plane_calc.display_reference_plane_info()
            
            # Update reference plane flag display
            self._update_reference_plane_flag()
    
    def _update_reference_plane_flag(self):
        """Add markers to atoms belonging to the primary reference plane."""
        reference_plane = self.reference_plane_calc.get_reference_plane()
        if not reference_plane:
            return
        
        reference_plane_id = reference_plane['plane_id']
        
        # Set reference plane status for all atoms
        for atom_data in self.atoms_data:
            if atom_data['Plane_ID'] == reference_plane_id:
                atom_data['Is_Reference_Plane'] = 'Yes'
            else:
                atom_data['Is_Reference_Plane'] = 'No'
    
    def get_reference_plane_calculator(self):
        """Return the primary reference plane calculator."""
        return self.reference_plane_calc
    
    def close_csv(self):
        """Finalize CSV operations."""
        # Display plane statistics
        self.display_plane_statistics()
        
        self.atoms_data = []  # Memory cleanup
        print("\nPlane analysis completed.")
    
    @staticmethod
    def display_supercell_info(original_atoms, supercell_atoms, cell_vectors, supercell_vectors):
        """
        Display supercell information.
        
        Args:
            original_atoms (list): Original atom information
            supercell_atoms (list): Supercell atom information
            cell_vectors (dict): Original cell vectors
            supercell_vectors (dict): Supercell vectors
        """
        print("\n[ Supercell Generation (6x6x6 fixed) ]")
        print("[ Atom Count ]")
        print(f"Original structure: {len(original_atoms):4d} atoms")
        print(f"Supercell:          {len(supercell_atoms):4d} atoms")
        
        # Cell volume information
        original_volume = np.abs(np.dot(cell_vectors['a'], 
                                      np.cross(cell_vectors['b'], cell_vectors['c'])))
        print(f"[ Volume ] {original_volume:.4f} Å³")
        
        # Cell vector information
        print("[ Cell Vectors ]")
        for vector_name in ['a', 'b', 'c']:
            vector = cell_vectors[vector_name]
            print(f"{vector_name} = [{vector[0]:8.4f}, {vector[1]:8.4f}, {vector[2]:8.4f}] Å")
        
        # Lattice parameter information
        print("[ Lattice Parameters ]")
        print(f"a = {np.linalg.norm(cell_vectors['a']):8.4f} Å")
        print(f"b = {np.linalg.norm(cell_vectors['b']):8.4f} Å")
        print(f"c = {np.linalg.norm(cell_vectors['c']):8.4f} Å")
    
    @staticmethod
    def save_to_csv(atoms, filename):
        """
        Save atom information to CSV file.
        
        Args:
            atoms (list): List of dictionaries containing atom information
            filename (str): File name to save
        """
        # Check and add file extension
        if not filename.endswith('.csv'):
            filename += '.csv'
        
        # Create directory if it doesn't exist
        directory = os.path.dirname(filename)
        if directory and not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
        
        # Write CSV file
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Write header (including additional columns for Step3 support, index column added)
            writer.writerow(['index', 'element', 'atom_id', 'cart_x', 'cart_y', 'cart_z', 'frac_x', 'frac_y', 'frac_z', 'plane_id', 'plane_equation', 'reference_plane', 'plane_type'])
            
            # Write atom information
            for i, atom in enumerate(atoms, 1):
                atom_id = f"{atom['element']}{i}"
                writer.writerow([
                    i,  # index (starts from 1)
                    atom['element'],  # element (changed from atom_type)
                    atom_id,  # atom_id (element + index format)
                    f"{atom['cart_x']:.6f}",  # cart_x
                    f"{atom['cart_y']:.6f}",  # cart_y
                    f"{atom['cart_z']:.6f}",  # cart_z
                    f"{atom['x']:.6f}",  # frac_x
                    f"{atom['y']:.6f}",  # frac_y
                    f"{atom['z']:.6f}",   # frac_z
                    "",  # plane_id (filled in Step3)
                    "",  # plane_equation (filled in Step3)
                    "",  # reference_plane (filled in Step3)
                    ""   # plane_type (filled in Step3)
                ])
        
        print(f"\nSupercell coordinates saved to {filename}.")
    
    @staticmethod
    def display_composition_analysis(atoms):
        """
        Display atomic composition analysis results.
        
        Args:
            atoms (list): List of dictionaries containing atom information
        """
        print("\n============================================================")
        print("Atomic Composition Analysis")
        print("============================================================")
        
        # Count by element type
        composition = {}
        for atom in atoms:
            element = atom['element']
            composition[element] = composition.get(element, 0) + 1
        
        total_atoms = len(atoms)
        
        print("\nCount and percentage by element type:")
        for element, count in sorted(composition.items()):
            percentage = (count / total_atoms) * 100
            print(f"  {element}: {count} atoms ({percentage:.2f}%)")
        
        print(f"\nTotal atoms: {total_atoms} atoms")
    
    @staticmethod
    def display_detailed_supercell_statistics(original_atoms, supercell_atoms, cif_data, supercell_size):
        """
        Display detailed supercell generation statistics.
        (Functionality separated from Phase2)
        
        Args:
            original_atoms (list): Original atom information
            supercell_atoms (list): Supercell atom information
            cif_data (dict): Original CIF data
            supercell_size (tuple): Supercell size (nx, ny, nz)
        """
        nx, ny, nz = supercell_size
        multiplier = nx * ny * nz
        print(f"\nSupercell Statistics:")
        print(f"   - Expansion factor: {nx} × {ny} × {nz} = {multiplier}x")
        print(f"   - Unit cell atoms: {len(original_atoms)}")
        print(f"   - Supercell atoms: {len(supercell_atoms)}")
        print(f"   - Expected atoms: {len(original_atoms) * multiplier}")
        
        # Element-wise statistics
        original_composition = {}
        supercell_composition = {}
        
        for atom in original_atoms:
            element = atom['element']
            original_composition[element] = original_composition.get(element, 0) + 1
            
        for atom in supercell_atoms:
            element = atom['element']
            supercell_composition[element] = supercell_composition.get(element, 0) + 1
        
        print(f"\nAtom Count by Element:")
        print(f"   {'Element':<8} {'Unit Cell':<8} {'Supercell':<10} {'Factor':<8}")
        print(f"   {'-'*8} {'-'*8} {'-'*10} {'-'*8}")
        
        for element in sorted(original_composition.keys()):
            unit_count = original_composition[element]
            super_count = supercell_composition.get(element, 0)
            element_multiplier = super_count // unit_count if unit_count > 0 else 0
            print(f"   {element:<8} {unit_count:<8} {super_count:<10} {element_multiplier:<8}")
            
        # Volume information (if available)
        if 'volume' in cif_data:
            original_volume = cif_data['volume']
            supercell_volume = original_volume * multiplier
            print(f"\nVolume Information:")
            print(f"   - Unit cell volume: {original_volume:.4f} Å³")
            print(f"   - Supercell volume: {supercell_volume:.4f} Å³") 