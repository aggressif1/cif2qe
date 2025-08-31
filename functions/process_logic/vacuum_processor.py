"""
Phase4 vacuum layer processing related functions
"""
import numpy as np
import copy

class VacuumProcessor:
    """Phase4 vacuum layer processing class"""
    
    def __init__(self):
        """Initialize Phase4VacuumProcessor"""
        pass
    
    def add_vacuum_layers(self, expanded_supercell, vacuum_config):
        """
        Add vacuum layers to supercell.
        
        Args:
            expanded_supercell (dict): Expanded supercell data
            vacuum_config (dict): Vacuum layer configuration
            
        Returns:
            dict: Structure data with added vacuum layers
        """
        try:
            print(f"Vacuum layer addition information:")
            print(f"   - Coordinate system: {vacuum_config['coordinate_system']}")
            print(f"   - Vacuum layer settings: {vacuum_config['vacuum_layers']}")
            
            # Copy original data
            result_structure = copy.deepcopy(expanded_supercell)
            
            # Process according to coordinate system
            if vacuum_config['coordinate_system'] == 'cartesian':
                result_structure = self._add_cartesian_vacuum(result_structure, vacuum_config)
            else:
                result_structure = self._add_user_defined_vacuum(result_structure, vacuum_config)
            
            if result_structure:
                print("SUCCESS: Vacuum layer addition completed")
                return result_structure
            else:
                print("ERROR: Vacuum layer addition failed")
                return None
                
        except Exception as e:
            print(f"ERROR: Error during vacuum layer addition: {str(e)}")
            return None
    
    def _add_cartesian_vacuum(self, structure, vacuum_config):
        """Add vacuum layers in Cartesian coordinate system"""
        try:
            vacuum_layers = vacuum_config['vacuum_layers']
            
            # Calculate current structure boundaries
            atoms = structure.get('atoms', [])
            if not atoms:
                print("ERROR: No atomic data available.")
                return None
            
            # Extract atomic coordinates
            coords = []
            for atom in atoms:
                x = atom.get('x', 0.0)
                y = atom.get('y', 0.0)
                z = atom.get('z', 0.0)
                coords.append([x, y, z])
            
            coords = np.array(coords)
            
            # Calculate current boundaries
            min_coords = np.min(coords, axis=0)
            max_coords = np.max(coords, axis=0)
            
            print(f"   Current structure boundaries:")
            print(f"     X: {min_coords[0]:.4f} ~ {max_coords[0]:.4f} Å")
            print(f"     Y: {min_coords[1]:.4f} ~ {max_coords[1]:.4f} Å")
            print(f"     Z: {min_coords[2]:.4f} ~ {max_coords[2]:.4f} Å")
            
            # Add vacuum layers (positive direction only)
            vacuum_x = vacuum_layers.get('x', 0.0)
            vacuum_y = vacuum_layers.get('y', 0.0)
            vacuum_z = vacuum_layers.get('z', 0.0)
            
            # Calculate new cell size
            new_cell_size = [
                max_coords[0] - min_coords[0] + vacuum_x,
                max_coords[1] - min_coords[1] + vacuum_y,
                max_coords[2] - min_coords[2] + vacuum_z
            ]
            
            print(f"   Cell size after vacuum layer addition:")
            print(f"     X: {new_cell_size[0]:.4f} Å (+{vacuum_x:.2f} Å)")
            print(f"     Y: {new_cell_size[1]:.4f} Å (+{vacuum_y:.2f} Å)")
            print(f"     Z: {new_cell_size[2]:.4f} Å (+{vacuum_z:.2f} Å)")
            
            # Update cell vectors (Cartesian coordinate system)
            structure['cell_vectors'] = {
                'a': [new_cell_size[0], 0.0, 0.0],
                'b': [0.0, new_cell_size[1], 0.0],
                'c': [0.0, 0.0, new_cell_size[2]]
            }
            
            # Keep atomic coordinates as is (vacuum layers added in positive direction only)
            structure['vacuum_info'] = {
                'coordinate_system': 'cartesian',
                'vacuum_layers': vacuum_layers,
                'original_cell_size': [
                    max_coords[0] - min_coords[0],
                    max_coords[1] - min_coords[1],
                    max_coords[2] - min_coords[2]
                ],
                'new_cell_size': new_cell_size
            }
            
            return structure
            
        except Exception as e:
            print(f"ERROR: Error during Cartesian vacuum layer addition: {str(e)}")
            return None
    
    def _add_user_defined_vacuum(self, structure, vacuum_config):
        """Add vacuum layers in user-defined cell vectors"""
        try:
            vacuum_layers = vacuum_config['vacuum_layers']
            
            # Get user-defined supercell vectors from new supercell data structure
            user_supercell_vectors = structure.get('user_defined_supercell_vectors')
            
            if user_supercell_vectors:
                # Use newly generated user-defined supercell vectors from Phase 3
                a_user_extended = user_supercell_vectors.get('a_user_extended')
                b_user_extended = user_supercell_vectors.get('b_user_extended')
                c_user_extended = user_supercell_vectors.get('c_user_extended')
                
                print(f"   Using user-defined supercell vectors:")
                print(f"     a_user_extended: {a_user_extended}")
                print(f"     b_user_extended: {b_user_extended}")
                print(f"     c_user_extended: {c_user_extended}")
                
                # Set working vectors
                a_work = np.array(a_user_extended)
                b_work = np.array(b_user_extended)
                c_work = np.array(c_user_extended)
                
            else:
                # Legacy method: Use actual supercell vectors (compatibility)
                actual_supercell = structure.get('actual_supercell')
                if actual_supercell:
                    a_work = np.array(actual_supercell.get('a_supercell'))
                    b_work = np.array(actual_supercell.get('b_supercell'))
                    c_work = np.array(actual_supercell.get('c_supercell'))
                    
                    print(f"   Using actual supercell vectors (compatibility mode):")
                    print(f"     a_supercell: {a_work}")
                    print(f"     b_supercell: {b_work}")
                    print(f"     c_supercell: {c_work}")
                else:
                    # Legacy method: Use a_user, b_user, c_user directly
                    a_user = structure.get('a_user')
                    b_user = structure.get('b_user')
                    c_user = structure.get('c_user')
                    
                    # Safe check for numpy arrays
                    vectors_valid = []
                    for vector in [a_user, b_user, c_user]:
                        if vector is not None:
                            # Check if numpy array or list and has appropriate size
                            if hasattr(vector, '__len__') and len(vector) >= 3:
                                vectors_valid.append(True)
                            else:
                                vectors_valid.append(False)
                        else:
                            vectors_valid.append(False)
                    
                    if not all(vectors_valid):
                        print("ERROR: No user-defined cell vectors available.")
                        return None
                    
                    a_work = np.array(a_user)
                    b_work = np.array(b_user)
                    c_work = np.array(c_user)
                    
                    print(f"   Using legacy user-defined cell vectors:")
                    print(f"     a_user: [{a_work[0]:.4f}, {a_work[1]:.4f}, {a_work[2]:.4f}]")
                    print(f"     b_user: [{b_work[0]:.4f}, {b_work[1]:.4f}, {b_work[2]:.4f}]")
                    print(f"     c_user: [{c_work[0]:.4f}, {c_work[1]:.4f}, {c_work[2]:.4f}]")
            
            # Add vacuum layers
            vacuum_a = vacuum_layers.get('a_user', 0.0)
            vacuum_b = vacuum_layers.get('b_user', 0.0)
            vacuum_c = vacuum_layers.get('c_user', 0.0)
            
            # Calculate new cell vectors (existing vectors + vacuum layers)
            # Extend in each direction by vacuum layer amount
            new_a = a_work + (a_work / np.linalg.norm(a_work)) * vacuum_a
            new_b = b_work + (b_work / np.linalg.norm(b_work)) * vacuum_b
            new_c = c_work + (c_work / np.linalg.norm(c_work)) * vacuum_c
            
            print(f"   Cell vectors after vacuum layer addition:")
            print(f"     new_a: [{new_a[0]:.4f}, {new_a[1]:.4f}, {new_a[2]:.4f}] (+{vacuum_a:.2f} Å)")
            print(f"     new_b: [{new_b[0]:.4f}, {new_b[1]:.4f}, {new_b[2]:.4f}] (+{vacuum_b:.2f} Å)")
            print(f"     new_c: [{new_c[0]:.4f}, {new_c[1]:.4f}, {new_c[2]:.4f}] (+{vacuum_c:.2f} Å)")
            
            # Update cell vectors
            structure['cell_vectors'] = {
                'a': new_a.tolist(),
                'b': new_b.tolist(),
                'c': new_c.tolist()
            }
            
            # Update user-defined vectors as well (for compatibility)
            structure['a_user'] = new_a.tolist()
            structure['b_user'] = new_b.tolist()
            structure['c_user'] = new_c.tolist()
            
            # Keep atomic coordinates as is (vacuum layers added in positive direction only)
            structure['vacuum_info'] = {
                'coordinate_system': 'user_defined',
                'vacuum_layers': vacuum_layers,
                'original_vectors': {
                    'a_work': a_work.tolist(),
                    'b_work': b_work.tolist(),
                    'c_work': c_work.tolist()
                },
                'new_vectors': {
                    'a_user': new_a.tolist(),
                    'b_user': new_b.tolist(),
                    'c_user': new_c.tolist()
                }
            }
            
            return structure
            
        except Exception as e:
            print(f"ERROR: Error during user-defined vacuum layer addition: {str(e)}")
            return None 