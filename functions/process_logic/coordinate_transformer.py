"""
Phase4 coordinate transformation related functions
"""
import numpy as np

class CoordinateTransformer:
    """Phase4 coordinate transformation class"""
    
    def __init__(self):
        """Initialize Phase4CoordinateTransformer"""
        pass
    
    def transform_coordinates(self, structure, vacuum_config):
        """
        Transform coordinates if necessary.
        
        Args:
            structure (dict): Structure data with added vacuum layer
            vacuum_config (dict): Vacuum layer configuration
            
        Returns:
            dict: Structure data with transformed coordinates
        """
        try:
            # Currently no special coordinate transformation needed
            # Atom coordinates remain unchanged when adding vacuum layer, only cell is expanded
            
            print("Checking coordinate transformation...")
            
            # Check number of atoms
            atoms = structure.get('atoms', [])
            print(f"   - Number of atoms: {len(atoms)}")
            
            # Check cell vectors
            cell_vectors = structure.get('cell_vectors', {})
            if cell_vectors:
                print("   - Cell vector update completed")
                for key, vector in cell_vectors.items():
                    print(f"     {key}: [{vector[0]:.4f}, {vector[1]:.4f}, {vector[2]:.4f}]")
            
            # Check vacuum layer information
            vacuum_info = structure.get('vacuum_info', {})
            if vacuum_info:
                print("   - Vacuum layer information saved completed")
            
            print("SUCCESS: Coordinate transformation completed (transformation not necessary)")
            return structure
            
        except Exception as e:
            print(f"ERROR: Error during coordinate transformation: {str(e)}")
            return structure  # Return original on error
    
    def convert_to_fractional_coordinates(self, structure):
        """
        Convert Cartesian coordinates to fractional coordinates.
        
        Args:
            structure (dict): Structure data
            
        Returns:
            dict: Structure data with added fractional coordinates
        """
        try:
            atoms = structure.get('atoms', [])
            cell_vectors = structure.get('cell_vectors', {})
            
            if not atoms or not cell_vectors:
                print("ERROR: No atom data or cell vectors available.")
                return structure
            
            # Construct cell vector matrix
            a = np.array(cell_vectors.get('a', [1, 0, 0]))
            b = np.array(cell_vectors.get('b', [0, 1, 0]))
            c = np.array(cell_vectors.get('c', [0, 0, 1]))
            
            cell_matrix = np.column_stack([a, b, c])
            
            # Calculate inverse matrix
            if np.linalg.det(cell_matrix) == 0:
                print("ERROR: Cell vector matrix is singular.")
                return structure
            
            inv_cell_matrix = np.linalg.inv(cell_matrix)
            
            # Calculate fractional coordinates for each atom
            for atom in atoms:
                cart_coords = np.array([
                    atom.get('x', 0.0),
                    atom.get('y', 0.0),
                    atom.get('z', 0.0)
                ])
                
                frac_coords = inv_cell_matrix @ cart_coords
                
                atom['frac_x'] = frac_coords[0]
                atom['frac_y'] = frac_coords[1]
                atom['frac_z'] = frac_coords[2]
            
            print("SUCCESS: Fractional coordinate conversion completed")
            return structure
            
        except Exception as e:
            print(f"ERROR: Error during fractional coordinate conversion: {str(e)}")
            return structure 