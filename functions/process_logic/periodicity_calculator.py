import numpy as np

class PeriodicityCalculator:
    """Class for handling periodic boundary conditions"""
    
    def __init__(self, cell_vectors=None):
        """
        Initialize PeriodicityCalculator
        
        Args:
            cell_vectors (dict): Cell vector information {'a': [x,y,z], 'b': [x,y,z], 'c': [x,y,z]}
        """
        self.cell_vectors = cell_vectors
        if cell_vectors:
            self.a_vec = np.array(cell_vectors['a'])
            self.b_vec = np.array(cell_vectors['b'])
            self.c_vec = np.array(cell_vectors['c'])
        else:
            self.a_vec = None
            self.b_vec = None
            self.c_vec = None
    
    def set_cell_vectors(self, cell_vectors):
        """
        Set cell vector information.
        
        Args:
            cell_vectors (dict): Cell vector information
        """
        self.cell_vectors = cell_vectors
        if cell_vectors:
            self.a_vec = np.array(cell_vectors['a'])
            self.b_vec = np.array(cell_vectors['b'])
            self.c_vec = np.array(cell_vectors['c'])
    
    def create_periodic_images(self, atoms, range_a=1, range_b=1, range_c=1):
        """
        Generate periodic images within given range.
        
        Args:
            atoms (list): Original atom list
            range_a, range_b, range_c (int): Unit cell range to consider in each direction
        
        Returns:
            list: Extended atom list including original + periodic images
        """
        if not self.cell_vectors:
            return atoms
        
        extended_atoms = []
        
        # Consider all unit cells within specified range
        for da in range(-range_a, range_a + 1):
            for db in range(-range_b, range_b + 1):
                for dc in range(-range_c, range_c + 1):
                    # Calculate cell displacement vector
                    cell_shift = da * self.a_vec + db * self.b_vec + dc * self.c_vec
                    
                    # Add all atoms from this cell
                    for i, atom in enumerate(atoms):
                        shifted_atom = atom.copy()
                        shifted_atom['cart_x'] += cell_shift[0]
                        shifted_atom['cart_y'] += cell_shift[1]
                        shifted_atom['cart_z'] += cell_shift[2]
                        
                        # Add periodic image information
                        shifted_atom['original_index'] = i
                        shifted_atom['cell_shift'] = (da, db, dc)
                        shifted_atom['is_original'] = (da == 0 and db == 0 and dc == 0)
                        
                        extended_atoms.append(shifted_atom)
        
        return extended_atoms
    
    def create_neighbor_shell(self, atoms, max_distance=10.0):
        """
        Generate all periodic images within specified maximum distance.
        
        Args:
            atoms (list): Original atom list
            max_distance (float): Maximum distance to consider (Å)
        
        Returns:
            list: Extended atom list
        """
        if not self.cell_vectors:
            return atoms
        
        # Calculate required number of cells in each direction
        a_length = np.linalg.norm(self.a_vec)
        b_length = np.linalg.norm(self.b_vec)
        c_length = np.linalg.norm(self.c_vec)
        
        range_a = int(np.ceil(max_distance / a_length)) + 1
        range_b = int(np.ceil(max_distance / b_length)) + 1
        range_c = int(np.ceil(max_distance / c_length)) + 1
        
        return self.create_periodic_images(atoms, range_a, range_b, range_c)
    
    def find_minimum_image_distance(self, atom1, atom2):
        """
        Calculate minimum image distance between two atoms. (Fractional coordinate based - intuitive and efficient)
        
        Args:
            atom1, atom2 (dict): Atomic information
        
        Returns:
            tuple: (minimum distance, minimum distance vector)
        """
        if not self.cell_vectors:
            # General distance calculation when no periodicity
            dx = atom2['cart_x'] - atom1['cart_x']
            dy = atom2['cart_y'] - atom1['cart_y']
            dz = atom2['cart_z'] - atom1['cart_z']
            distance = np.sqrt(dx*dx + dy*dy + dz*dz)
            return distance, (dx, dy, dz)
        
        # 1. Convert Cartesian coordinates to fractional coordinates
        pos1 = np.array([atom1['cart_x'], atom1['cart_y'], atom1['cart_z']])
        pos2 = np.array([atom2['cart_x'], atom2['cart_y'], atom2['cart_z']])
        
        # Create cell matrix (arranged as column vectors)
        cell_matrix = np.column_stack([self.a_vec, self.b_vec, self.c_vec])
        
        # Convert to fractional coordinates (using inverse matrix)
        try:
            cell_inv = np.linalg.inv(cell_matrix)
            frac1 = cell_inv @ pos1
            frac2 = cell_inv @ pos2
        except np.linalg.LinAlgError:
            # Fallback to existing method if inverse matrix calculation fails
            return self._fallback_minimum_image_distance(atom1, atom2)
        
        # 2. Calculate shortest distance vector in fractional coordinates (Key: wrap to ±0.5 range)
        frac_diff = frac2 - frac1
        
        # Consider periodic boundaries in each direction (normalize to -0.5 ~ +0.5 range)
        for i in range(3):
            if frac_diff[i] > 0.5:
                frac_diff[i] -= 1.0
            elif frac_diff[i] < -0.5:
                frac_diff[i] += 1.0
        
        # 3. Convert fractional coordinate difference back to Cartesian coordinates
        cart_diff = cell_matrix @ frac_diff
        
        # 4. Calculate distance
        distance = np.linalg.norm(cart_diff)
        
        return distance, tuple(cart_diff)
    
    def _fallback_minimum_image_distance(self, atom1, atom2):
        """
        Fallback method used when fractional coordinate conversion fails (3x3x3 check)
        """
        min_distance = float('inf')
        min_vector = None
        
        # Basic 3x3x3 check
        for da in [-1, 0, 1]:
            for db in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    cell_shift = da * self.a_vec + db * self.b_vec + dc * self.c_vec
                    
                    shifted_x = atom2['cart_x'] + cell_shift[0]
                    shifted_y = atom2['cart_y'] + cell_shift[1]
                    shifted_z = atom2['cart_z'] + cell_shift[2]
                    
                    dx = shifted_x - atom1['cart_x']
                    dy = shifted_y - atom1['cart_y']
                    dz = shifted_z - atom1['cart_z']
                    distance = np.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    if distance < min_distance:
                        min_distance = distance
                        min_vector = (dx, dy, dz)
        
        return min_distance, min_vector
    
    def get_cell_info(self):
        """
        Return cell information.
        
        Returns:
            dict: Cell information
        """
        if not self.cell_vectors:
            return None
        
        return {
            'vectors': self.cell_vectors,
            'lengths': {
                'a': np.linalg.norm(self.a_vec),
                'b': np.linalg.norm(self.b_vec),
                'c': np.linalg.norm(self.c_vec)
            },
            'volume': np.abs(np.dot(self.a_vec, np.cross(self.b_vec, self.c_vec)))
        } 