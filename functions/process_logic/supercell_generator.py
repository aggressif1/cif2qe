import numpy as np

class SupercellGenerator:
    """Class for generating supercells"""
    
    @staticmethod
    def calculate_supercell_vectors(cell_vectors, nx, ny, nz):
        """
        Calculate lattice vectors of supercell.
        
        Args:
            cell_vectors (dict): Unit cell vectors {'a': [...], 'b': [...], 'c': [...]}
            nx, ny, nz (int): Number of repetitions in each direction
            
        Returns:
            dict: Supercell vectors
        """
        return {
            'a': [v * nx for v in cell_vectors['a']],
            'b': [v * ny for v in cell_vectors['b']],
            'c': [v * nz for v in cell_vectors['c']]
        }
    
    @staticmethod
    def generate_supercell_atoms(atoms, nx, ny, nz, cell_vectors=None):
        """
        Generate atomic positions in supercell.
        
        Args:
            atoms (list): List of dictionaries containing atomic information
            nx, ny, nz (int): Number of repetitions in each direction
            cell_vectors (dict): Unit cell vector information
            
        Returns:
            list: Atomic information in supercell
        """
        supercell_atoms = []
        atom_index = 1  # Global index for all atoms
        
        # Convert cell vectors to numpy arrays
        if cell_vectors:
            a_vec = np.array(cell_vectors['a'])
            b_vec = np.array(cell_vectors['b'])
            c_vec = np.array(cell_vectors['c'])
        
        # Replicate unit cell in each direction
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for atom in atoms:
                        # Create new atom information
                        new_atom = atom.copy()
                        
                        # Add atom_id (element + global index)
                        new_atom['atom_id'] = f"{atom['element']}{atom_index}"
                        new_atom['index'] = atom_index
                        atom_index += 1
                        
                        # Calculate fractional coordinates (new fractional coordinates in supercell)
                        new_atom['x'] = (atom['x'] + ix) / nx
                        new_atom['y'] = (atom['y'] + iy) / ny
                        new_atom['z'] = (atom['z'] + iz) / nz
                        
                        # Calculate Cartesian coordinates
                        if cell_vectors:
                            # Calculate Cartesian coordinates using precise vector operations
                            fractional_pos = np.array([new_atom['x'], new_atom['y'], new_atom['z']])
                            cartesian_pos = (fractional_pos[0] * a_vec * nx + 
                                           fractional_pos[1] * b_vec * ny + 
                                           fractional_pos[2] * c_vec * nz)
                            
                            new_atom['cart_x'] = cartesian_pos[0]
                            new_atom['cart_y'] = cartesian_pos[1]
                            new_atom['cart_z'] = cartesian_pos[2]
                        else:
                            # Simple calculation based on unit cell translation
                            # Add lattice vectors to original atom's Cartesian coordinates in unit cell
                            unit_a = atom['cart_x'] / atom['x'] if atom['x'] != 0 else 0
                            unit_b = atom['cart_y'] / atom['y'] if atom['y'] != 0 else 0
                            unit_c = atom['cart_z'] / atom['z'] if atom['z'] != 0 else 0
                            
                            new_atom['cart_x'] = atom['cart_x'] + ix * unit_a
                            new_atom['cart_y'] = atom['cart_y'] + iy * unit_b
                            new_atom['cart_z'] = atom['cart_z'] + iz * unit_c
                        
                        supercell_atoms.append(new_atom)
        
        # Remove duplicate atoms at periodic boundaries (e.g., coordinate = 1.0 when 0.0 exists)
        supercell_atoms = SupercellGenerator._remove_duplicate_boundary_atoms(supercell_atoms)
        
        return supercell_atoms
    
    @staticmethod
    def _remove_duplicate_boundary_atoms(atoms):
        """
        Remove duplicate atoms at periodic boundaries.
        
        When fractional coordinates include 1.0, they are equivalent to 0.0 due to 
        periodic boundary conditions. This method removes atoms with coordinate 1.0
        when corresponding atoms with 0.0 exist at the same position.
        
        Args:
            atoms (list): List of atomic information
            
        Returns:
            list: Atomic information with duplicates removed
        """
        tolerance = 1e-8  # Tolerance for floating point comparison
        atoms_to_keep = []
        
        for i, atom in enumerate(atoms):
            x, y, z = atom['x'], atom['y'], atom['z']
            
            # Check if any coordinate is approximately 1.0
            has_boundary_coord = (
                abs(x - 1.0) < tolerance or 
                abs(y - 1.0) < tolerance or 
                abs(z - 1.0) < tolerance
            )
            
            if has_boundary_coord:
                # Convert boundary coordinates to 0.0 equivalent
                equiv_x = 0.0 if abs(x - 1.0) < tolerance else x
                equiv_y = 0.0 if abs(y - 1.0) < tolerance else y  
                equiv_z = 0.0 if abs(z - 1.0) < tolerance else z
                
                # Check if equivalent atom exists
                equivalent_exists = False
                for j, other_atom in enumerate(atoms):
                    if i != j and atom['element'] == other_atom['element']:
                        if (abs(other_atom['x'] - equiv_x) < tolerance and
                            abs(other_atom['y'] - equiv_y) < tolerance and  
                            abs(other_atom['z'] - equiv_z) < tolerance):
                            equivalent_exists = True
                            break
                
                # Keep atom only if no equivalent exists
                if not equivalent_exists:
                    atoms_to_keep.append(atom)
            else:
                # Keep atoms that don't have boundary coordinates
                atoms_to_keep.append(atom)
        
        original_count = len(atoms)
        final_count = len(atoms_to_keep)
        if original_count != final_count:
            print(f"   Removed {original_count - final_count} duplicate boundary atoms")
            print(f"   Atoms: {original_count} → {final_count}")
        
        return atoms_to_keep
    
    @staticmethod
    def add_vacuum(atoms, cell_vectors, vacuum_size, direction='z'):
        """
        Add vacuum region in specific direction.
        
        Args:
            atoms (list): List of dictionaries containing atomic information
            cell_vectors (dict): Cell vectors
            vacuum_size (float): Size of vacuum region to add (Å)
            direction (str): Direction to add vacuum ('x', 'y', or 'z')
            
        Returns:
            tuple: (modified atomic information, modified cell vectors)
        """
        if vacuum_size <= 0:
            return atoms, cell_vectors
        
        # Calculate current cell size
        cell_size = {
            'x': np.linalg.norm(cell_vectors['a']),
            'y': np.linalg.norm(cell_vectors['b']),
            'z': np.linalg.norm(cell_vectors['c'])
        }
        
        # Modify cell vectors for vacuum addition
        new_vectors = cell_vectors.copy()
        if direction == 'x':
            scale = (cell_size['x'] + vacuum_size) / cell_size['x']
            new_vectors['a'] = [v * scale for v in cell_vectors['a']]
        elif direction == 'y':
            scale = (cell_size['y'] + vacuum_size) / cell_size['y']
            new_vectors['b'] = [v * scale for v in cell_vectors['b']]
        else:  # z direction
            scale = (cell_size['z'] + vacuum_size) / cell_size['z']
            new_vectors['c'] = [v * scale for v in cell_vectors['c']]
        
        # Don't change atomic positions (use fractional coordinates)
        return atoms, new_vectors
    
    @staticmethod
    def center_atoms(atoms, cell_vectors):
        """
        Move atoms to center of cell.
        
        Args:
            atoms (list): List of dictionaries containing atomic information
            cell_vectors (dict): Cell vectors
            
        Returns:
            list: Atomic information moved to center
        """
        # Calculate center of mass of atoms
        center = {'x': 0, 'y': 0, 'z': 0}
        total_atoms = len(atoms)
        
        for atom in atoms:
            center['x'] += atom['cart_x']
            center['y'] += atom['cart_y']
            center['z'] += atom['cart_z']
        
        center = {k: v/total_atoms for k, v in center.items()}
        
        # Calculate center of cell
        cell_center = {
            'x': sum(cell_vectors['a'])/2,
            'y': sum(cell_vectors['b'])/2,
            'z': sum(cell_vectors['c'])/2
        }
        
        # Calculate translation vector
        shift = {k: cell_center[k] - center[k] for k in ['x', 'y', 'z']}
        
        # Translate atomic positions
        centered_atoms = []
        for atom in atoms:
            new_atom = atom.copy()
            new_atom['cart_x'] += shift['x']
            new_atom['cart_y'] += shift['y']
            new_atom['cart_z'] += shift['z']
            
            # Update fractional coordinates
            new_atom['x'] = new_atom['cart_x'] / np.linalg.norm(cell_vectors['a'])
            new_atom['y'] = new_atom['cart_y'] / np.linalg.norm(cell_vectors['b'])
            new_atom['z'] = new_atom['cart_z'] / np.linalg.norm(cell_vectors['c'])
            
            centered_atoms.append(new_atom)
        
        return centered_atoms 