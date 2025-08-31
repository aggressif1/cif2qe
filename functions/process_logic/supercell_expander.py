"""
Functions related to supercell expansion separated from Phase3
"""
import numpy as np
import os
import json

class SupercellExpander:
    """Supercell expansion class for Phase3"""
    
    def __init__(self):
        """Initialize Phase3SupercellExpander"""
        pass
    
    def create_supercell_expansion(self, unit_cell_group):
        """
        Creates supercell expansion based on unit cell group.
        
        Args:
            unit_cell_group: Unit cell group information
            
        Returns:
            dict: Expanded supercell information
        """
        print("\n🔄 Step 5: Create supercell expansion")
        print("=" * 80)
        
        try:
            # 1. Get supercell size input
            supercell_size = self.get_supercell_size_input()
            if not supercell_size:
                print("❌ Failed to input supercell size")
                return None
            
            # 2. Generate expanded supercell
            expanded_supercell = self.generate_expanded_supercell(unit_cell_group, supercell_size)
            if not expanded_supercell:
                print("❌ Failed to create supercell expansion")
                return None
            
            # 3. Display result summary
            self.display_supercell_expansion_summary(expanded_supercell)
            
            # 4. Save results
            self.save_expanded_supercell(expanded_supercell)
            
            return expanded_supercell
            
        except Exception as e:
            print(f"❌ Failed to create supercell expansion: {str(e)}")
            return None
    
    def get_supercell_size_input(self):
        """
        Gets supercell size input from the user.
        
        Returns:
            dict: Supercell size information or None
        """
        try:
            print("\n📏 Supercell size configuration:")
            print("   Format: 'nx ny nz' (e.g., '2 2 2' or '3 3 1')")
            print("   Specifies how many unit cells to replicate in each direction.")
            
            while True:
                try:
                    user_input = input("Enter supercell size: ").strip()
                    
                    if not user_input:
                        print("❌ Input is empty. Please enter again.")
                        continue
                    
                    supercell_size = self.parse_supercell_size(user_input)
                    if supercell_size:
                        print(f"✅ Set supercell size: {supercell_size['nx']} × {supercell_size['ny']} × {supercell_size['nz']}")
                        return supercell_size
                    else:
                        print("❌ Please enter in correct format. (e.g., '2 2 2')")
                        
                except KeyboardInterrupt:
                    print("\n❌ User cancelled input.")
                    return None
                except Exception as e:
                    print(f"❌ Input processing error: {str(e)}")
                    continue
                    
        except Exception as e:
            print(f"❌ Failed to input supercell size: {str(e)}")
            return None
    
    def parse_supercell_size(self, user_input):
        """
        Parses user input to extract supercell size.
        
        Args:
            user_input (str): User input string
            
        Returns:
            dict: Parsed supercell size or None
        """
        try:
            # Split by whitespace
            parts = user_input.split()
            
            if len(parts) == 1:
                # If single number, apply equally to all directions
                size = int(parts[0])
                if size > 0:
                    return {'nx': size, 'ny': size, 'nz': size}
            elif len(parts) == 3:
                # If 3 numbers
                nx, ny, nz = int(parts[0]), int(parts[1]), int(parts[2])
                if nx > 0 and ny > 0 and nz > 0:
                    return {'nx': nx, 'ny': ny, 'nz': nz}
            
            return None
            
        except ValueError:
            return None
        except Exception as e:
            print(f"⚠️ Parsing error: {str(e)}")
            return None
    
    def generate_expanded_supercell(self, unit_cell_group, supercell_size):
        """
        Expands unit cell to generate supercell.
        Actual supercell is generated with user input + 1,
        and user-defined supercell cell vectors are extended only by user input amount.
        Expands based on fractional coordinates to solve boundary issues.
        
        Args:
            unit_cell_group: Unit cell group information
            supercell_size: User input supercell size {'nx': int, 'ny': int, 'nz': int}
            
        Returns:
            dict: Expanded supercell information (actual supercell + user-defined supercell cell vectors)
        """
        try:
            # User input size
            user_nx = supercell_size['nx']
            user_ny = supercell_size['ny'] 
            user_nz = supercell_size['nz']
            
            # Actual supercell size (user input + 1)
            actual_nx = user_nx + 1
            actual_ny = user_ny + 1
            actual_nz = user_nz + 1
            
            print(f"\n🔧 Starting supercell expansion:")
            print(f"   User input size: {user_nx}×{user_ny}×{user_nz}")
            print(f"   Actual generation size: {actual_nx}×{actual_ny}×{actual_nz}")
            
            # Unit cell vectors
            a_user = np.array(unit_cell_group['a_user'])
            b_user = np.array(unit_cell_group['b_user'])
            c_user = np.array(unit_cell_group['c_user'])
            
            # User-defined supercell cell vectors (extended by user input amount)
            a_user_extended = user_nx * a_user
            b_user_extended = user_ny * b_user
            c_user_extended = user_nz * c_user
            
            # Actual supercell vectors (actual generation size)
            actual_expanded_a = actual_nx * a_user
            actual_expanded_b = actual_ny * b_user
            actual_expanded_c = actual_nz * c_user
            
            print(f"📐 User-defined supercell cell vectors:")
            print(f"   a_user_extended = {user_nx} × a_user = ({a_user_extended[0]:.4f}, {a_user_extended[1]:.4f}, {a_user_extended[2]:.4f})")
            print(f"   b_user_extended = {user_ny} × b_user = ({b_user_extended[0]:.4f}, {b_user_extended[1]:.4f}, {b_user_extended[2]:.4f})")
            print(f"   c_user_extended = {user_nz} × c_user = ({c_user_extended[0]:.4f}, {c_user_extended[1]:.4f}, {c_user_extended[2]:.4f})")
            
            expanded_atoms = []
            atom_count = 0
            
            # Fractional coordinate-based expansion (to actual supercell size)
            for i in range(actual_nx):
                for j in range(actual_ny):
                    for k in range(actual_nz):
                        for atom in unit_cell_group['atoms']:
                            # Read fractional coordinates (using uvw key)
                            if 'uvw' in atom:
                                frac_u, frac_v, frac_w = atom['uvw']
                            else:
                                # Use default values if uvw not available
                                frac_u, frac_v, frac_w = 0.0, 0.0, 0.0
                            
                            # Calculate new fractional coordinates (based on actual supercell)
                            new_frac_u = (frac_u + i) / actual_nx
                            new_frac_v = (frac_v + j) / actual_ny
                            new_frac_w = (frac_w + k) / actual_nz
                            
                            # Convert to Cartesian coordinates (using actual supercell vectors)
                            new_pos = (new_frac_u * actual_expanded_a + 
                                     new_frac_v * actual_expanded_b + 
                                     new_frac_w * actual_expanded_c)
                            
                            # Create new atom
                            new_atom = atom.copy()
                            new_atom['x'] = new_pos[0]
                            new_atom['y'] = new_pos[1]
                            new_atom['z'] = new_pos[2]
                            
                            # Add position information within supercell
                            new_atom['supercell_index'] = (i, j, k)
                            new_atom['original_frac'] = (frac_u, frac_v, frac_w)
                            new_atom['new_frac'] = (new_frac_u, new_frac_v, new_frac_w)
                            
                            expanded_atoms.append(new_atom)
                            atom_count += 1
            
            # Compose supercell information
            expanded_supercell = {
                'origin': unit_cell_group['origin'],
                
                # Actually generated supercell information
                'actual_supercell': {
                    'a_supercell': actual_expanded_a.tolist(),
                    'b_supercell': actual_expanded_b.tolist(), 
                    'c_supercell': actual_expanded_c.tolist(),
                    'size': {'nx': actual_nx, 'ny': actual_ny, 'nz': actual_nz}
                },
                
                # User-defined supercell cell vectors (passed to Phase 4)
                'user_defined_supercell_vectors': {
                    'a_user_extended': a_user_extended.tolist(),
                    'b_user_extended': b_user_extended.tolist(),
                    'c_user_extended': c_user_extended.tolist(),
                    'user_input_size': {'nx': user_nx, 'ny': user_ny, 'nz': user_nz}
                },
                
                # Original unit cell vectors (for reference)
                'unit_cell_vectors': {
                    'a_user': a_user.tolist(),
                    'b_user': b_user.tolist(),
                    'c_user': c_user.tolist()
                },
                
                'atoms': expanded_atoms,
                'supercell_size': supercell_size,  # Maintain user input size
                'actual_supercell_size': {'nx': actual_nx, 'ny': actual_ny, 'nz': actual_nz},
                'total_atoms': atom_count,
                'unit_cell_atoms': len(unit_cell_group['atoms'])
            }
            
            # Atom group classification (for Phase 4)
            atom_groups = self._classify_atoms_by_fractional_coordinates(
                expanded_atoms, a_user_extended, b_user_extended, c_user_extended
            )
            
            if atom_groups:
                expanded_supercell['atom_groups'] = atom_groups
                print(f"📊 Atom group classification completed: bulk_group: {len(atom_groups['bulk_group'])} atoms | boundary_group: {len(atom_groups['boundary_group'])} atoms | Not_used_group: {len(atom_groups['Not_used_group'])} atoms (deleted)")
            
            print(f"✅ Supercell expansion completed:")
            print(f"   - Unit cell atom count: {len(unit_cell_group['atoms'])}")
            print(f"   - Actual supercell atom count: {atom_count}")
            print(f"   - Actual expansion ratio: {actual_nx}×{actual_ny}×{actual_nz} = {actual_nx*actual_ny*actual_nz}")
            print(f"   - User-defined expansion ratio: {user_nx}×{user_ny}×{user_nz} = {user_nx*user_ny*user_nz}")
            
            return expanded_supercell
            
        except Exception as e:
            print(f"❌ Error during supercell expansion: {e}")
            return None
    
    def _classify_atoms_by_fractional_coordinates(self, atoms, a_user_extended, b_user_extended, c_user_extended):
        """
        Classifies atoms according to fractional coordinates based on user-defined supercell cell vectors.
        
        Args:
            atoms: Atom list
            a_user_extended, b_user_extended, c_user_extended: User-defined supercell cell vectors
            
        Returns:
            dict: Classified atom groups
        """
        try:
            print(f"\n🔍 Starting atom group classification...")
            
            # Convert user-defined supercell cell vectors to numpy arrays
            a_vec = np.array(a_user_extended)
            b_vec = np.array(b_user_extended)
            c_vec = np.array(c_user_extended)
            
            # Create transformation matrix (Cartesian → fractional coordinates)
            cell_matrix = np.column_stack([a_vec, b_vec, c_vec])
            
            # Calculate inverse matrix (for fractional coordinate transformation)
            try:
                inv_cell_matrix = np.linalg.inv(cell_matrix)
            except np.linalg.LinAlgError:
                print("❌ Cell vector matrix is singular. Skipping classification.")
                return None
            
            # Initialize atom groups
            bulk_group = []
            boundary_group = []
            not_used_group = []
            
            tolerance = 1e-6  # Numerical error tolerance
            
            for atom in atoms:
                # Extract Cartesian coordinates
                cart_pos = np.array([atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)])
                
                # Convert to fractional coordinates
                frac_pos = inv_cell_matrix @ cart_pos
                
                # Store fractional coordinates (for debugging)
                atom['frac_coords_extended'] = frac_pos.tolist()
                
                # Check classification criteria
                is_bulk = True
                is_boundary = False
                is_not_used = False
                
                for coord in frac_pos:
                    if coord > 1.0 + tolerance:
                        is_not_used = True
                        is_bulk = False
                        break
                    elif abs(coord - 1.0) <= tolerance:
                        is_boundary = True
                        is_bulk = False
                    elif coord < -tolerance:
                        is_not_used = True
                        is_bulk = False
                        break
                
                # Group classification
                if is_not_used:
                    not_used_group.append(atom)
                elif is_boundary:
                    boundary_group.append(atom)
                else:
                    bulk_group.append(atom)
            
            print(f"   Classification results:")
            print(f"     - bulk_group: {len(bulk_group)} atoms (0 ≤ coordinates < 1)")
            print(f"     - boundary_group: {len(boundary_group)} atoms (coordinates = 1)")
            print(f"     - Not_used_group: {len(not_used_group)} atoms (coordinates > 1, deleted)")
            
            # Debug: Print fractional coordinates of a few atoms
            print(f"   Sample fractional coordinates (first 3 atoms):")
            for i, atom in enumerate(atoms[:3]):
                if 'frac_coords_extended' in atom:
                    frac = atom['frac_coords_extended']
                    element = atom.get('element', 'X')
                    print(f"     {element}: ({frac[0]:.6f}, {frac[1]:.6f}, {frac[2]:.6f})")
            
            return {
                'bulk_group': bulk_group,
                'boundary_group': boundary_group,
                'Not_used_group': not_used_group,  # Not passed to Phase 4
                'classification_info': {
                    'total_atoms': len(atoms),
                    'bulk_count': len(bulk_group),
                    'boundary_count': len(boundary_group),
                    'not_used_count': len(not_used_group),
                    'tolerance': tolerance
                }
            }
            
        except Exception as e:
            print(f"❌ Error during atom group classification: {str(e)}")
            return None
    
    def display_supercell_expansion_summary(self, expanded_supercell):
        """
        Displays summary of supercell expansion results.
        
        Args:
            expanded_supercell: Expanded supercell information
        """
        try:
            print(f"\n📊 Supercell expansion results summary:")
            
            size = expanded_supercell['supercell_size']
            print(f"   Expansion size: {size['nx']} × {size['ny']} × {size['nz']}")
            print(f"   Unit cell atom count: {expanded_supercell['unit_cell_atoms']}")
            print(f"   Total atom count: {expanded_supercell['total_atoms']}")
            print(f"   Expansion ratio: {expanded_supercell['total_atoms'] / expanded_supercell['unit_cell_atoms']:.1f}x")
            
            # Calculate count by element
            element_counts = {}
            for atom in expanded_supercell['atoms']:
                element = atom['element']
                element_counts[element] = element_counts.get(element, 0) + 1
            
            print(f"   Count by element:")
            for element, count in sorted(element_counts.items()):
                print(f"     {element}: {count}")
            
            # Supercell vector information
            supercell_vectors = expanded_supercell['supercell_vectors']
            print(f"   Supercell vectors:")
            print(f"     a_expanded: ({supercell_vectors['a_expanded'][0]:.4f}, {supercell_vectors['a_expanded'][1]:.4f}, {supercell_vectors['a_expanded'][2]:.4f})")
            print(f"     b_expanded: ({supercell_vectors['b_expanded'][0]:.4f}, {supercell_vectors['b_expanded'][1]:.4f}, {supercell_vectors['b_expanded'][2]:.4f})")
            print(f"     c_expanded: ({supercell_vectors['c_expanded'][0]:.4f}, {supercell_vectors['c_expanded'][1]:.4f}, {supercell_vectors['c_expanded'][2]:.4f})")
            
        except Exception as e:
            print(f"⚠️ Failed to display summary information: {str(e)}")
    
    def save_expanded_supercell(self, expanded_supercell):
        """
        Saves the expanded supercell to a file.
        
        Args:
            expanded_supercell: Expanded supercell information
        """
        try:
            print(f"\n💾 Saving supercell results:")
            
            # Check output directory
            output_dir = 'output'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # Save as JSON file (convert numpy arrays to lists)
            json_data = {}
            for key, value in expanded_supercell.items():
                if key == 'unit_cell_vectors' or key == 'supercell_vectors':
                    json_data[key] = {}
                    for sub_key, sub_value in value.items():
                        if isinstance(sub_value, np.ndarray):
                            json_data[key][sub_key] = sub_value.tolist()
                        else:
                            json_data[key][sub_key] = sub_value
                elif isinstance(value, np.ndarray):
                    json_data[key] = value.tolist()
                else:
                    json_data[key] = value
            
            json_path = os.path.join(output_dir, 'expanded_supercell.json')
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(json_data, f, indent=2, ensure_ascii=False)
            
            print(f"   JSON file saved: {json_path}")
            
            # Also save as CSV file (atom coordinates)
            csv_path = os.path.join(output_dir, 'expanded_supercell_coordinates.csv')
            with open(csv_path, 'w', encoding='utf-8') as f:
                f.write("element,x,y,z,cell_index_i,cell_index_j,cell_index_k\n")
                
                for atom in expanded_supercell['atoms']:
                    cell_idx = atom.get('cell_index', (0, 0, 0))
                    f.write(f"{atom['element']},{atom['x']:.6f},{atom['y']:.6f},{atom['z']:.6f},{cell_idx[0]},{cell_idx[1]},{cell_idx[2]}\n")
            
            print(f"   CSV file saved: {csv_path}")
            print(f"✅ Supercell results saved successfully")
            
        except Exception as e:
            print(f"❌ Failed to save supercell results: {str(e)}") 