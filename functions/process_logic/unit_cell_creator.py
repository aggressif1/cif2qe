"""
Unit cell creation functions separated from Phase3
"""
import numpy as np
import csv
import os

class UnitCellCreator:
    """Phase3 dedicated unit cell creation class"""
    
    def __init__(self):
        """Initialize Phase3UnitCellCreator"""
        pass
    
    def create_new_unit_cell_group(self, a_user, b_user, c_user):
        """
        Create new unit cell group.
        
        Args:
            a_user: a_user vector
            b_user: b_user vector  
            c_user: c_user vector
            
        Returns:
            dict: Unit cell group information
        """
        print("\nStep 4: Creating new unit cell group")
        print("=" * 80)
        
        try:
            # 1. Find center point of first reference plane
            first_plane_center = self.find_first_reference_plane_center()
            if first_plane_center is None:
                print("ERROR: Cannot find first reference plane center point")
                return None
            
            # 2. Determine new origin
            new_origin = self.determine_new_origin(first_plane_center, a_user, b_user, c_user)
            if new_origin is None:
                print("ERROR: Cannot determine new origin")
                return None
            
            # 3. Find atoms inside unit cell
            unit_cell_atoms = self.find_atoms_in_unit_cell(new_origin, a_user, b_user, c_user)
            if not unit_cell_atoms:
                print("ERROR: Cannot find atoms inside unit cell")
                return None
            
            # 4. Exclude boundary atoms
            filtered_atoms = self.exclude_boundary_atoms(unit_cell_atoms, new_origin, a_user, b_user, c_user)
            
            unit_cell_group = {
                'origin': new_origin,
                'a_user': a_user,
                'b_user': b_user,
                'c_user': c_user,
                'atoms': filtered_atoms,
                'first_plane_center': first_plane_center
            }
            
            # 5. Display result summary
            self.display_unit_cell_group_summary(unit_cell_group)
            
            return unit_cell_group
            
        except Exception as e:
            print(f"ERROR: Unit cell group creation failed: {str(e)}")
            return None
    
    def find_first_reference_plane_center(self):
        """
        Find center point of first reference plane.
        
        Returns:
            numpy.ndarray: Center point or None
        """
        try:
            # Find atoms of first reference plane from CSV
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("ERROR: CSV file not found")
                return None
            
            first_plane_atoms = []
            
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    if row.get('reference_plane') == 'First reference plane':
                        try:
                            x = float(row.get('cart_x', 0))
                            y = float(row.get('cart_y', 0))
                            z = float(row.get('cart_z', 0))
                            first_plane_atoms.append([x, y, z])
                        except ValueError:
                            continue
            
            if not first_plane_atoms:
                print("ERROR: Cannot find atoms of first reference plane")
                return None
            
            # Calculate center point (average)
            center = np.mean(first_plane_atoms, axis=0)
            
            print(f"First reference plane center point calculation: {len(first_plane_atoms)} atoms in first reference plane | Center point: ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f})")
            
            return center
            
        except Exception as e:
            print(f"ERROR: First reference plane center point calculation failed: {str(e)}")
            return None
    
    def determine_new_origin(self, first_plane_center, a_user, b_user, c_user):
        """
        Determine new origin.
        
        Args:
            first_plane_center: First reference plane center point
            a_user, b_user, c_user: User-defined cell vectors
            
        Returns:
            numpy.ndarray: New origin
        """
        try:
            print(f"\nDetermining new origin:")
            print(f"   First reference plane center point: ({first_plane_center[0]:.4f}, {first_plane_center[1]:.4f}, {first_plane_center[2]:.4f})")
            
            # Find atom closest to center point in first reference plane
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("ERROR: CSV file not found")
                return None
            
            closest_atom_pos = None
            min_distance = float('inf')
            
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    if row.get('reference_plane') == 'First reference plane':
                        try:
                            x = float(row.get('cart_x', 0))
                            y = float(row.get('cart_y', 0))
                            z = float(row.get('cart_z', 0))
                            atom_pos = np.array([x, y, z])
                            
                            # Calculate distance to center point
                            distance = np.linalg.norm(atom_pos - first_plane_center)
                            
                            if distance < min_distance:
                                min_distance = distance
                                closest_atom_pos = atom_pos
                                
                        except ValueError:
                            continue
            
            if closest_atom_pos is None:
                print("ERROR: Cannot find closest atom")
                return None
            
            print(f"   Closest atom position: ({closest_atom_pos[0]:.4f}, {closest_atom_pos[1]:.4f}, {closest_atom_pos[2]:.4f})")
            print(f"   Distance to center point: {min_distance:.4f} Å")
            print(f"   New origin: ({closest_atom_pos[0]:.4f}, {closest_atom_pos[1]:.4f}, {closest_atom_pos[2]:.4f})")
            
            return closest_atom_pos
            
        except Exception as e:
            print(f"ERROR: New origin determination failed: {str(e)}")
            return None
    
    def find_atoms_in_unit_cell(self, origin, a_user, b_user, c_user):
        """
        Find atoms inside unit cell.
        
        Args:
            origin: Unit cell origin
            a_user, b_user, c_user: Unit cell vectors
            
        Returns:
            list: List of atoms inside unit cell
        """
        try:
            print(f"\nSearching for atoms inside unit cell:")
            
            # Read all atom information from CSV
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("ERROR: CSV file not found")
                return []
            
            all_atoms = []
            
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    try:
                        atom_info = {
                            'element': row.get('element', 'Unknown'),
                            'cart_x': float(row.get('cart_x', 0)),  # Preserve original coordinates
                            'cart_y': float(row.get('cart_y', 0)),
                            'cart_z': float(row.get('cart_z', 0)),
                            'plane_id': row.get('plane_id', ''),
                            'reference_plane': row.get('reference_plane', ''),
                            'plane_type': row.get('plane_type', '')
                        }
                        all_atoms.append(atom_info)
                    except ValueError:
                        continue
            
            print(f"   Total atoms loaded: {len(all_atoms)}")
            
            # Find atoms inside unit cell using coordinate transformation
            unit_cell_atoms = []
            tolerance = 1e-6  # Tolerance for boundary judgment
            
            # Convert to numpy arrays
            a_vec = np.array(a_user)
            b_vec = np.array(b_user)
            c_vec = np.array(c_user)
            
            # Create transformation matrix (cell vectors as columns)
            cell_matrix = np.column_stack([a_vec, b_vec, c_vec])
            
            # Check if matrix is invertible
            if np.abs(np.linalg.det(cell_matrix)) < 1e-10:
                print("ERROR: Cell vectors are linearly dependent")
                return []
            
            inv_cell_matrix = np.linalg.inv(cell_matrix)
            
            for atom in all_atoms:
                # Relative position from origin
                position = np.array([atom['cart_x'], atom['cart_y'], atom['cart_z']])
                relative_pos = position - origin
                
                # Transform to fractional coordinates
                fractional_coords = inv_cell_matrix @ relative_pos
                
                # Check if inside unit cell (0 <= coord < 1)
                inside_cell = all(
                    -tolerance <= coord < 1.0 - tolerance 
                    for coord in fractional_coords
                )
                
                if inside_cell:
                    # Add fractional coordinates to atom information
                    atom['frac_x'] = fractional_coords[0]
                    atom['frac_y'] = fractional_coords[1]
                    atom['frac_z'] = fractional_coords[2]
                    unit_cell_atoms.append(atom)
            
            print(f"   Atoms found inside unit cell: {len(unit_cell_atoms)}")
            
            return unit_cell_atoms
            
        except Exception as e:
            print(f"ERROR: Unit cell atom search failed: {str(e)}")
            return []
    
    def exclude_boundary_atoms(self, unit_cell_atoms, origin, a_user, b_user, c_user):
        """
        Exclude boundary atoms.
        
        Args:
            unit_cell_atoms: Atoms inside unit cell
            origin: Origin
            a_user, b_user, c_user: Cell vectors
            
        Returns:
            list: Atoms with boundary atoms excluded
        """
        try:
            print(f"\nExcluding boundary atoms:")
            
            boundary_tolerance = 1e-4  # Tolerance for boundary determination
            filtered_atoms = []
            
            for atom in unit_cell_atoms:
                frac_coords = [atom['frac_x'], atom['frac_y'], atom['frac_z']]
                
                # Check if it's a boundary atom
                is_boundary = any(
                    abs(coord) < boundary_tolerance or abs(coord - 1.0) < boundary_tolerance
                    for coord in frac_coords
                )
                
                if not is_boundary:
                    filtered_atoms.append(atom)
            
            excluded_count = len(unit_cell_atoms) - len(filtered_atoms)
            print(f"   Original atoms: {len(unit_cell_atoms)}")
            print(f"   Boundary atoms excluded: {excluded_count}")
            print(f"   Final atoms: {len(filtered_atoms)}")
            
            return filtered_atoms
            
        except Exception as e:
            print(f"ERROR: Boundary atom exclusion failed: {str(e)}")
            return unit_cell_atoms  # Return original if failed
    
    def display_unit_cell_group_summary(self, unit_cell_group):
        """
        Display unit cell group summary.
        
        Args:
            unit_cell_group: Unit cell group information
        """
        try:
            print(f"\nUnit cell group creation completed:")
            print("=" * 60)
            
            origin = unit_cell_group['origin']
            atoms = unit_cell_group['atoms']
            
            print(f"Origin: ({origin[0]:.4f}, {origin[1]:.4f}, {origin[2]:.4f})")
            print(f"Number of atoms: {len(atoms)}")
            
            # Element statistics
            from ..utils.element_counter import ElementCounter
            element_counts = ElementCounter.count_elements(atoms)
            print(f"Element composition: {element_counts}")
            
            # Cell vector information
            a_user = unit_cell_group['a_user']
            b_user = unit_cell_group['b_user']
            c_user = unit_cell_group['c_user']
            
            print(f"Cell vectors:")
            print(f"   a_user: [{a_user[0]:.4f}, {a_user[1]:.4f}, {a_user[2]:.4f}]")
            print(f"   b_user: [{b_user[0]:.4f}, {b_user[1]:.4f}, {b_user[2]:.4f}]")
            print(f"   c_user: [{c_user[0]:.4f}, {c_user[1]:.4f}, {c_user[2]:.4f}]")
            
            # Calculate unit cell volume
            a_vec = np.array(a_user)
            b_vec = np.array(b_user)
            c_vec = np.array(c_user)
            volume = abs(np.dot(a_vec, np.cross(b_vec, c_vec)))
            print(f"Unit cell volume: {volume:.4f} Å³")
            
        except Exception as e:
            print(f"ERROR: Unit cell group summary display failed: {str(e)}") 