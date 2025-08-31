"""
Phase5 hydrogen atom placement related functions
"""
import numpy as np
import random
import copy

class HydrogenPlacer:
    """Phase5 hydrogen atom placement class"""
    
    def __init__(self):
        """Initialize Phase5HydrogenPlacer"""
        pass
    
    def define_placement_region(self, structure, hydrogen_config):
        """
        Define hydrogen placement regions.
        
        Args:
            structure (dict): Structure data
            hydrogen_config (dict): Hydrogenation configuration
            
        Returns:
            dict: Placement region information
        """
        try:
            atoms = structure.get('atoms', [])
            if not atoms:
                print("ERROR: No atomic data available.")
                return None
            
            # Extract atomic coordinates (prioritize moved final coordinates after vacuum layer insertion in slab structure)
            coords = []
            for atom in atoms:
                # Prioritize final coordinates moved due to vacuum layer insertion in slab structure
                if 'final_x' in atom and 'final_y' in atom and 'final_z' in atom:
                    x = atom['final_x']
                    y = atom['final_y']
                    z = atom['final_z']
                else:
                    # fallback: use existing coordinates
                    x = atom.get('x', 0.0)
                    y = atom.get('y', 0.0)
                    z = atom.get('z', 0.0)
                coords.append([x, y, z])
            
            coords = np.array(coords)
            
            # Calculate structure boundaries (based on coordinates after vacuum layer addition)
            min_coords = np.min(coords, axis=0)
            max_coords = np.max(coords, axis=0)
            
            print(f"   Structure boundaries (after vacuum layer addition):")
            print(f"     X: {min_coords[0]:.4f} ~ {max_coords[0]:.4f} Å")
            print(f"     Y: {min_coords[1]:.4f} ~ {max_coords[1]:.4f} Å")
            print(f"     Z: {min_coords[2]:.4f} ~ {max_coords[2]:.4f} Å")
            
            # Get user-defined cell vector information (use as received from Phase 4)
            cell_vectors = structure.get('cell_vectors', {})
            if not cell_vectors:
                print("ERROR: No user-defined cell vector information available.")
                return None
            
            print(f"   Available cell vectors (final version with vacuum layers added):")
            # Separate basic and extended cell vectors
            basic_vectors = {k: v for k, v in cell_vectors.items() if k in ['a', 'b', 'c']}
            extended_vectors = {k: v for k, v in cell_vectors.items() if k not in ['a', 'b', 'c']}
            
            if basic_vectors:
                basic_info = " | ".join([f"{name}: length {(vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5:.2f} Å" for name, vector in basic_vectors.items()])
                print(f"     Basic: {basic_info}")
            
            if extended_vectors:
                extended_info = " | ".join([f"{name}: length {(vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5:.2f} Å" for name, vector in extended_vectors.items()])
                print(f"     Extended: {extended_info}")
            
            # Define hydrogenation regions (based on new configuration structure)
            direction_config = hydrogen_config['direction_config']
            
            placement_regions = {}
            
            # Process each user-defined cell vector direction
            for vector_name, directions in direction_config.items():
                # Get corresponding cell vector information
                if vector_name not in cell_vectors:
                    print(f"WARNING: Cannot find {vector_name} cell vector information.")
                    continue
                
                cell_vector = np.array(cell_vectors[vector_name])
                
                # Process positive direction
                pos_config = directions['positive']
                if pos_config['enabled']:
                    region_name = f"{vector_name}_positive"
                    region_data = self._define_vector_direction_region(
                        coords, cell_vector, pos_config, True, region_name
                    )
                    if region_data:
                        placement_regions[region_name] = region_data
                
                # Process negative direction
                neg_config = directions['negative']
                if neg_config['enabled']:
                    region_name = f"{vector_name}_negative"
                    region_data = self._define_vector_direction_region(
                        coords, cell_vector, neg_config, False, region_name
                    )
                    if region_data:
                        placement_regions[region_name] = region_data
            
            if not placement_regions:
                print("ERROR: No enabled hydrogenation directions.")
                return None
            
            print(f"   SUCCESS: {len(placement_regions)} hydrogenation regions defined")
            
            return {
                'regions': placement_regions,
                'structure_bounds': {
                    'min_coords': min_coords.tolist(),
                    'max_coords': max_coords.tolist()
                }
            }
            
        except Exception as e:
            print(f"ERROR: Error during hydrogen placement region definition: {str(e)}")
            return None
    
    def _define_vector_direction_region(self, coords, cell_vector, range_config, is_positive, region_name):
        """Define hydrogenation region for specific cell vector direction"""
        try:
            min_range = range_config['min_range']
            max_range = range_config['max_range']
            
            # Calculate cell vector direction unit vector
            vector_length = np.linalg.norm(cell_vector)
            if vector_length == 0:
                print(f"WARNING: {region_name}: Cannot process zero vector.")
                return None
            
            unit_vector = cell_vector / vector_length
            
            # Calculate structure center point
            center = np.mean(coords, axis=0)
            
            # Define region according to direction
            if is_positive:
                # Positive direction: from center towards cell vector direction
                region_start = center + unit_vector * min_range
                region_end = center + unit_vector * max_range
            else:
                # Negative direction: from center towards opposite of cell vector direction
                region_start = center - unit_vector * max_range
                region_end = center - unit_vector * min_range
            
            # Calculate region boundaries (extend structure boundaries in other directions)
            min_coords = np.min(coords, axis=0)
            max_coords = np.max(coords, axis=0)
            
            # Determine main direction axis (largest component)
            main_axis = np.argmax(np.abs(unit_vector))
            
            # Set region boundaries
            region_min = min_coords.copy()
            region_max = max_coords.copy()
            
            # Use configured range for main direction
            if is_positive:
                region_min[main_axis] = center[main_axis] + min_range
                region_max[main_axis] = center[main_axis] + max_range
            else:
                region_min[main_axis] = center[main_axis] - max_range
                region_max[main_axis] = center[main_axis] - min_range
            
            # Slightly extend structure boundaries in other directions
            for i in range(3):
                if i != main_axis:
                    margin = (max_coords[i] - min_coords[i]) * 0.1  # 10% margin
                    region_min[i] -= margin
                    region_max[i] += margin
            
            print(f"     {region_name}: {min_range:.1f}~{max_range:.1f} Å region defined")
            
            return {
                'min_bound': region_min.tolist(),
                'max_bound': region_max.tolist(),
                'vector_direction': unit_vector.tolist(),
                'is_positive': is_positive,
                'range_config': range_config
            }
            
        except Exception as e:
            print(f"ERROR: Error during {region_name} region definition: {str(e)}")
            return None
    
    def place_hydrogen_atoms(self, structure, placement_region, hydrogen_config):
        """
        Place hydrogen atoms at random positions.
        
        Args:
            structure (dict): Structure data
            placement_region (dict): Placement region information
            hydrogen_config (dict): Hydrogenation configuration
            
        Returns:
            dict: Structure with placed hydrogen atoms
        """
        try:
            print("\nStarting hydrogen atom placement...")
            
            # Deep copy original structure
            result_structure = copy.deepcopy(structure)
            
            # Initialize placement count
            total_placed = 0
            placement_details = {}
            
            # Process each region
            for region_name, region_data in placement_region['regions'].items():
                # Extract region configuration
                range_config = region_data['range_config']
                atom_count = range_config.get('count', 0)
                
                if atom_count <= 0:
                    print(f"   {region_name}: Skipping (count: {atom_count})")
                    continue
                
                # Place hydrogen atoms in this region
                placed_atoms = self._place_in_region(region_data, atom_count, region_name)
                
                if placed_atoms:
                    # Add to main atom list
                    if 'atoms' not in result_structure:
                        result_structure['atoms'] = []
                    
                    result_structure['atoms'].extend(placed_atoms)
                    total_placed += len(placed_atoms)
                    placement_details[region_name] = len(placed_atoms)
                    
                    print(f"   {region_name}: {len(placed_atoms)} atoms placed")
                else:
                    print(f"   {region_name}: Failed to place atoms")
                    placement_details[region_name] = 0
            
            # Update hydrogen placement information
            result_structure['hydrogen_placement'] = {
                'total_placed': total_placed,
                'placement_details': placement_details,
                'total_atoms': len(result_structure.get('atoms', []))
            }
            
            print(f"\nSUCCESS: Hydrogen placement completed")
            print(f"   - Total hydrogen atoms placed: {total_placed}")
            print(f"   - Total atoms (including hydrogen): {len(result_structure.get('atoms', []))}")
            
            return result_structure
            
        except Exception as e:
            print(f"ERROR: Error during hydrogen atom placement: {str(e)}")
            return None
    
    def _place_in_region(self, region_data, count, region_name):
        """
        Place specified number of hydrogen atoms in given region.
        
        Args:
            region_data (dict): Region data
            count (int): Number of atoms to place
            region_name (str): Region name
            
        Returns:
            list: List of placed hydrogen atoms
        """
        try:
            min_bound = np.array(region_data['min_bound'])
            max_bound = np.array(region_data['max_bound'])
            
            placed_atoms = []
            
            for i in range(count):
                # Generate random position within region boundaries
                x = random.uniform(min_bound[0], max_bound[0])
                y = random.uniform(min_bound[1], max_bound[1])
                z = random.uniform(min_bound[2], max_bound[2])
                
                # Create hydrogen atom data
                hydrogen_atom = {
                    'element': 'H',
                    'x': x,
                    'y': y,
                    'z': z,
                    'cart_x': x,  # Same as fractional coordinates
                    'cart_y': y,
                    'cart_z': z,
                    'atom_id': f"H_{region_name}_{i+1}",
                    'placement_region': region_name,
                    'is_hydrogen': True
                }
                
                placed_atoms.append(hydrogen_atom)
            
            return placed_atoms
            
        except Exception as e:
            print(f"ERROR: Error placing atoms in {region_name}: {str(e)}")
            return [] 