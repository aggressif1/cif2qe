"""
Functions related to cell vector calculation separated from Phase3
"""
import numpy as np
import csv
import os

class CellVectorCalculator:
    """Cell vector calculation class"""

    def __init__(self):
        """Initialize CellVectorCalculator"""
        pass
    
    def initialize_user_cell_vectors(self):
        """
        Initializes user-defined cell vectors.
        
        Returns:
            dict: Initialized cell vector information
        """
        print("\n🔧 Step 3-1: Initialize user-defined cell vectors")
        print("-" * 50)
        
        user_cell_vectors = {
            'a_user': None,
            'b_user': None, 
            'c_user': None
        }
        
        print("✅ User-defined cell vector initialization completed")
        return user_cell_vectors
    
    def calculate_miller_transformation_matrix(self, miller_indices):
        """
        Calculates transformation matrix from initial Miller indices (001) to user-defined Miller indices.
        
        Args:
            miller_indices (dict): User-defined Miller indices
            
        Returns:
            numpy.ndarray: Transformation matrix or None
        """
        print("\n🔧 Step 3-2: Calculate transformation matrix")
        print("-" * 50)
        
        try:
            # Get user-defined Miller indices
            h_user = miller_indices['h']
            k_user = miller_indices['k'] 
            l_user = miller_indices['l']
            
            print(f"📐 Initial Miller indices: (0 0 1)")
            print(f"📐 User Miller indices: ({h_user} {k_user} {l_user})")
            
            # Simple transformation matrix (start with identity matrix)
            transformation_matrix = np.eye(3)
            
            print(f"📊 Transformation matrix:")
            for i, row in enumerate(transformation_matrix):
                print(f"   [{row[0]:8.4f} {row[1]:8.4f} {row[2]:8.4f}]")
            
            return transformation_matrix
            
        except Exception as e:
            print(f"❌ Transformation matrix calculation failed: {str(e)}")
            return None
    
    def calculate_in_plane_vectors(self, transformation_matrix):
        """
        Defines a_user, b_user vectors based on actual atom arrangement in the plane.
        
        Args:
            transformation_matrix: Not used (kept for compatibility)
            
        Returns:
            tuple: (a_user, b_user) vectors or (None, None)
        """
        print("\n🔲 Step 3-2-1 & 3-2-2: Define vectors based on atom arrangement in plane")
        print("-" * 50)
        
        try:
            # 1. Get atoms in the first reference plane
            plane_atoms = self.get_first_reference_plane_atoms()
            if not plane_atoms or len(plane_atoms) < 2:
                print("❌ Insufficient atoms in the first reference plane")
                return None, None
            
            print(f"📊 Number of atoms in first reference plane: {len(plane_atoms)}")
            
            # 2. Calculate plane center and select reference atom
            plane_center = self.calculate_plane_center(plane_atoms)
            print(f"📍 Plane center point: ({plane_center[0]:.4f}, {plane_center[1]:.4f}, {plane_center[2]:.4f})")
            
            # Select the atom closest to the center as the reference point
            reference_atom = self.find_reference_atom(plane_center, plane_atoms)
            if reference_atom is None:
                print("❌ Cannot find reference atom")
                return None, None
            
            print(f"📍 Reference atom: ({reference_atom[0]:.4f}, {reference_atom[1]:.4f}, {reference_atom[2]:.4f})")
            
            # 3. Define a_user: vector from reference atom to periodic repetition unit
            a_user = self.define_a_user_vector(reference_atom, plane_atoms)
            if a_user is None:
                print("❌ Failed to define a_user vector")
                return None, None
            
            print(f"📐 a_user = ({a_user[0]:.4f}, {a_user[1]:.4f}, {a_user[2]:.4f}), magnitude: {np.linalg.norm(a_user):.4f} Å")
            
            # 4. Define b_user: vector to the second closest atom
            b_user = self.define_b_user_vector(reference_atom, plane_atoms, a_user)
            if b_user is None:
                print("❌ Failed to define b_user vector")
                return None, None
            
            print(f"📐 b_user = ({b_user[0]:.4f}, {b_user[1]:.4f}, {b_user[2]:.4f}), magnitude: {np.linalg.norm(b_user):.4f} Å")
            
            # 5. Calculate angle between vectors
            cos_angle = np.dot(a_user, b_user) / (np.linalg.norm(a_user) * np.linalg.norm(b_user))
            cos_angle = np.clip(cos_angle, -1, 1)
            angle_deg = np.degrees(np.arccos(abs(cos_angle)))
            print(f"📐 Angle between a_user and b_user: {angle_deg:.2f}°")
            
            # 6. Verify periodicity
            periodicity_verified = self.verify_in_plane_periodicity(a_user, b_user, None)
            
            if periodicity_verified:
                print("✅ Vector definition based on atom arrangement in plane completed")
                return a_user, b_user
            else:
                print("⚠️ Periodicity verification failed, but continuing")
                return a_user, b_user
                
        except Exception as e:
            print(f"❌ Failed to define vectors in plane: {str(e)}")
            return None, None
    
    def get_first_reference_plane_atoms(self):
        """Gets atoms belonging to the first reference plane."""
        try:
            import csv
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            
            plane_atoms = []
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    if row.get('reference_plane') == '제1기준평면':
                        try:
                            x = float(row.get('cart_x', 0))
                            y = float(row.get('cart_y', 0))
                            z = float(row.get('cart_z', 0))
                            plane_atoms.append(np.array([x, y, z]))
                        except (ValueError, TypeError):
                            continue
            
            return plane_atoms
            
        except Exception as e:
            print(f"⚠️ Failed to get first reference plane atom information: {e}")
            return []
    
    def calculate_plane_center(self, plane_atoms):
        """Calculates the center point of atoms in the plane."""
        if not plane_atoms:
            return np.array([0.0, 0.0, 0.0])
        
        center = np.mean(plane_atoms, axis=0)
        return center
    
    def find_reference_atom(self, center, plane_atoms):
        """Finds the atom closest to the center point."""
        min_distance = float('inf')
        closest_atom = None
        
        for atom in plane_atoms:
            vector = atom - center
            distance = np.linalg.norm(vector)
            
            if distance < min_distance:
                min_distance = distance
                closest_atom = atom
        
        return closest_atom
    
    def define_a_user_vector(self, reference_atom, plane_atoms):
        """
        Defines a_user as the vector from reference atom to the closest atom.
        Ensures it has the minimum length for periodicity.
        
        Args:
            reference_atom: Reference atom coordinates
            plane_atoms: Atoms in the plane
            
        Returns:
            numpy.ndarray: a_user vector or None
        """
        try:
            # Calculate distances and vectors from reference point to other atoms
            distances = []
            vectors = []
            
            for atom in plane_atoms:
                vector = atom - reference_atom
                distance = np.linalg.norm(vector)
                
                # Exclude self (very small distance case)
                if distance > 1e-6:
                    distances.append(distance)
                    vectors.append(vector)
            
            if not distances:
                print("⚠️ Cannot find appropriate atom, using default vector")
                return np.array([1.0, 0.0, 0.0])  # Default value
            
            # Sort distance and vector pairs (ascending order by distance)
            distance_vector_pairs = list(zip(distances, vectors))
            distance_vector_pairs.sort(key=lambda x: x[0])
            
            # Select vector to the closest atom
            closest_distance, closest_vector = distance_vector_pairs[0]
            
            print(f"   📏 Distance to closest atom: {closest_distance:.4f} Å")
            
            # Check minimum length for periodicity
# If the closest distance is the interatomic bond distance (half of lattice constant),
# consider expanding by 2x for periodicity
            original_cell_vectors = self.get_original_cell_vectors()
            if original_cell_vectors:
                lattice_constant = np.linalg.norm(original_cell_vectors['a'])
                half_lattice = lattice_constant / 2
                
                # If the closest distance is about half of the lattice constant (interatomic bond distance)
                if abs(closest_distance - half_lattice) < 0.2:
                    # Find vector expanded by 2x for periodicity
                    target_distance = lattice_constant
                    for dist, vec in distance_vector_pairs:
                        if abs(dist - target_distance) < 0.2:
                            print(f"   ✅ 주기성을 위해 확장된 거리 선택: {dist:.4f} Å")
                            return vec
                    
                    # If expanded vector not found, directly expand by 2x
                    expanded_vector = closest_vector * 2
                    expanded_distance = np.linalg.norm(expanded_vector)
                    print(f"   ✅ 주기성을 위해 벡터 2배 확장: {closest_distance:.4f} Å → {expanded_distance:.4f} Å")
                    return expanded_vector
            
            # Return the closest vector by default
            return closest_vector
            
        except Exception as e:
            print(f"❌ Error during a_user vector definition: {e}")
            return None
    
    def define_b_user_vector(self, reference_atom, plane_atoms, a_user):
        """
        Defines b_user as the vector from reference atom to the second closest atom.
        Ensures it is independent from a_user and has minimum length for periodicity.
        
        Args:
            reference_atom: Reference atom coordinates
            plane_atoms: Atoms in the plane
            a_user: Already defined a_user vector
            
        Returns:
            numpy.ndarray: b_user vector or None
        """
        try:
            a_user_normalized = a_user / np.linalg.norm(a_user)
            
            # Calculate distances and vectors from reference point to other atoms
            distances = []
            vectors = []
            
            for atom in plane_atoms:
                vector = atom - reference_atom
                distance = np.linalg.norm(vector)
                
                # Exclude self (very small distance case)
                if distance <= 1e-6:
                    continue
                
                # Check independence from a_user
                vector_normalized = vector / distance
                cos_angle = abs(np.dot(a_user_normalized, vector_normalized))
                angle_deg = np.degrees(np.arccos(np.clip(cos_angle, 0, 1)))
                
                # Only consider vectors with angles between 30° and 150° (ensuring independence)
                if 30.0 <= angle_deg <= 150.0:
                    distances.append(distance)
                    vectors.append(vector)
            
            if not distances:
                print("⚠️ Cannot find independent atom, generating with cross product")
                # Calculate plane normal vector (temporary)
                normal = np.array([0, 0, 1])  # Temporary normal vector
                b_user = np.cross(normal, a_user)
                if np.linalg.norm(b_user) < 1e-6:
                    b_user = np.cross(np.array([1, 0, 0]), a_user)
                return b_user
            
            # Sort distance and vector pairs (ascending order by distance)
            distance_vector_pairs = list(zip(distances, vectors))
            distance_vector_pairs.sort(key=lambda x: x[0])
            
            # Select vector to the closest independent atom (second closest overall atom)
            closest_distance, closest_vector = distance_vector_pairs[0]
            
            print(f"   📏 Distance to second closest atom: {closest_distance:.4f} Å")
            
            # Check minimum length for periodicity (same logic as a_user)
            original_cell_vectors = self.get_original_cell_vectors()
            if original_cell_vectors:
                lattice_constant = np.linalg.norm(original_cell_vectors['a'])
                half_lattice = lattice_constant / 2
                
                # If the closest distance is about half of the lattice constant (interatomic bond distance)
                if abs(closest_distance - half_lattice) < 0.2:
                    # Find vector expanded by 2x for periodicity
                    target_distance = lattice_constant
                    for dist, vec in distance_vector_pairs:
                        if abs(dist - target_distance) < 0.2:
                            print(f"   ✅ 주기성을 위해 확장된 거리 선택: {dist:.4f} Å")
                            return vec
                    
                    # If expanded vector not found, directly expand by 2x
                    expanded_vector = closest_vector * 2
                    expanded_distance = np.linalg.norm(expanded_vector)
                    print(f"   ✅ 주기성을 위해 벡터 2배 확장: {closest_distance:.4f} Å → {expanded_distance:.4f} Å")
                    return expanded_vector
            
            # Return the closest independent vector by default
            return closest_vector
            
        except Exception as e:
            print(f"❌ Error during b_user vector definition: {e}")
            return None
    
    def verify_in_plane_periodicity(self, a_user, b_user, first_reference_plane):
        """
        Verifies that a_user and b_user are properly defined within the first reference plane.
        Checks if there are actual atoms at positions moved by each vector from the reference point.
        
        Args:
            a_user: a_user vector
            b_user: b_user vector
            first_reference_plane: First reference plane information
            
        Returns:
            bool: Verification result
        """
        print(f"\n🔍 In-plane periodicity verification:")
        print(f"   a_user = ({a_user[0]:.4f}, {a_user[1]:.4f}, {a_user[2]:.4f}), magnitude: {np.linalg.norm(a_user):.4f} Å")
        print(f"   b_user = ({a_user[0]:.4f}, {b_user[1]:.4f}, {b_user[2]:.4f}), magnitude: {np.linalg.norm(b_user):.4f} Å")
        print(f"   a_user · b_user: {np.dot(a_user, b_user):.4f}")
        
        try:
            # Get reference atom position
            reference_atom = self.get_reference_atom_position()
            if reference_atom is None:
                print("❌ Cannot find reference atom position")
                return False
            
            # Get all atom coordinates in the first reference plane
            plane_atoms = self.get_first_reference_plane_atoms()
            if not plane_atoms:
                print("❌ Cannot find atoms in first reference plane")
                return False
            
            print(f"   Reference atom: ({reference_atom[0]:.4f}, {reference_atom[1]:.4f}, {reference_atom[2]:.4f})")
            print(f"   Number of atoms in first reference plane: {len(plane_atoms)}")
            
            # Define vectors to verify
            test_vectors = [
                ("a_user", a_user),
                ("b_user", b_user),
                ("a_user + b_user", a_user + b_user),
                ("-a_user", -a_user),
                ("-b_user", -b_user)
            ]
            
            verification_results = []
            
            for vector_name, test_vector in test_vectors:
                target_position = reference_atom + test_vector
                
                # Check if there are atoms near the target position (0.1Å tolerance)
                found_atom = False
                closest_distance = float('inf')
                closest_atom = None
                
                for atom in plane_atoms:
                    distance = np.linalg.norm(target_position - atom)
                    if distance < closest_distance:
                        closest_distance = distance
                        closest_atom = atom
                    
                    if distance < 0.1:  # Consider atom exists if within 0.1Å
                        found_atom = True
                        break
                
                if found_atom:
                    print(f"   ✅ Atom exists at {vector_name} movement position (distance: {closest_distance:.4f} Å)")
                    verification_results.append(True)
                else:
                    print(f"   ❌ No atom at {vector_name} movement position (distance to closest atom: {closest_distance:.4f} Å)")
                    print(f"      Target position: ({target_position[0]:.4f}, {target_position[1]:.4f}, {target_position[2]:.4f})")
                    if closest_atom is not None:
                        print(f"      Closest atom: ({closest_atom[0]:.4f}, {closest_atom[1]:.4f}, {closest_atom[2]:.4f})")
                    verification_results.append(False)
            
            # Evaluate results
            passed_tests = sum(verification_results)
            total_tests = len(verification_results)
            
            print(f"\n📊 Periodicity verification results:")
            print(f"   Tests passed: {passed_tests}/{total_tests}")
            
            if passed_tests >= 3:  # Success if at least 3 tests pass
                print("✅ In-plane periodicity verification passed")
                return True
            else:
                print("❌ In-plane periodicity verification failed")
                return False
                
        except Exception as e:
            print(f"❌ Error during periodicity verification: {e}")
            return False
            
    def calculate_out_of_plane_vector(self, sequence_info, transformation_matrix):
        """
        Defines c_user as a vector to the next equivalent plane by moving in the Miller index direction (h*a + k*b + l*c).
        Uses a new origin on the first reference plane as the reference point.
        
        Args:
            sequence_info: Plane sequence information
            transformation_matrix: Not used (kept for compatibility)
            
        Returns:
            numpy.ndarray: c_user vector or None
        """
        print("\n🔲 Step 3-2-3: Calculate out-of-plane vector")
        print("-" * 50)
        
        try:
            # 1. Get Miller indices and original cell vector information
            miller_indices = self.get_miller_indices()
            original_cell_vectors = self.get_original_cell_vectors()
            
            if not miller_indices or not original_cell_vectors:
                print("❌ Cannot find Miller indices or cell vector information")
                return None
            
            h, k, l = miller_indices
            a_vec = np.array(original_cell_vectors['a'])
            b_vec = np.array(original_cell_vectors['b']) 
            c_vec = np.array(original_cell_vectors['c'])
            
            print(f"📐 Miller indices: ({h}, {k}, {l})")
            print(f"📐 Original cell vectors:")
            print(f"   a = ({a_vec[0]:.4f}, {a_vec[1]:.4f}, {a_vec[2]:.4f})")
            print(f"   b = ({b_vec[0]:.4f}, {b_vec[1]:.4f}, {b_vec[2]:.4f})")
            print(f"   c = ({c_vec[0]:.4f}, {c_vec[1]:.4f}, {c_vec[2]:.4f})")
            
            # 2. Calculate Miller index direction vector: c_direction = h*a + k*b + l*c
            c_direction = h * a_vec + k * b_vec + l * c_vec
            c_direction_normalized = c_direction / np.linalg.norm(c_direction)
            
            print(f"📐 Miller index direction vector: ({c_direction[0]:.4f}, {c_direction[1]:.4f}, {c_direction[2]:.4f})")
            print(f"📐 Normalized direction: ({c_direction_normalized[0]:.4f}, {c_direction_normalized[1]:.4f}, {c_direction_normalized[2]:.4f})")
            
            # 3. Get new origin coordinates (reference point on the first reference plane)
            reference_atom = self.get_reference_atom_position()
            if reference_atom is None:
                print("❌ Cannot find new origin coordinates")
                return None
            
            print(f"📍 New origin (first reference plane): ({reference_atom[0]:.4f}, {reference_atom[1]:.4f}, {reference_atom[2]:.4f})")
            
            # 4. Find next equivalent planes to the first reference plane
            equivalent_planes = self.get_equivalent_planes()
            if not equivalent_planes:
                print("⚠️ Cannot find equivalent plane information")
                return None
            
            # 5. Calculate distance to next equivalent plane in Miller index direction
            next_plane_distance = self.find_next_equivalent_plane_in_direction(
                reference_atom, c_direction_normalized, equivalent_planes
            )
            
            if next_plane_distance is None:
                print("❌ Cannot find next equivalent plane in Miller index direction")
                return None
            
            print(f"📏 Distance to next equivalent plane in Miller index direction: {next_plane_distance:.4f} Å")
            
            # 6. Define c_user vector: direction (Miller indices) × length (distance to next equivalent plane)
            c_user = c_direction_normalized * next_plane_distance
            
            print(f"📐 c_user = ({c_user[0]:.4f}, {c_user[1]:.4f}, {c_user[2]:.4f}), magnitude: {np.linalg.norm(c_user):.4f} Å")
            
            # 7. Verification: Check if c_user actually reaches the equivalent plane
            verification_passed = self.verify_c_user_reaches_equivalent_plane(
                reference_atom, c_user, equivalent_planes
            )
            
            if verification_passed:
                print("✅ Out-of-plane vector calculation completed")
                return c_user
            else:
                print("⚠️ Verification failed, but continuing")
                return c_user
                
        except Exception as e:
            print(f"❌ Failed to calculate out-of-plane vector: {str(e)}")
            return None
    
    def get_miller_indices(self):
        """Gets Miller indices."""
        try:
            # Load Miller indices information from Phase3 temporary file
            from ..utils.temporary_file_handler import TemporaryFileHandler
            temp_handler = TemporaryFileHandler()
            phase3_data = temp_handler.load_phase_results("phase3")
            
            if phase3_data and 'miller_indices' in phase3_data:
                miller_indices = phase3_data['miller_indices']
                
                # Handle various forms of Miller indices
                if isinstance(miller_indices, dict):
                    # Case: {'h': 1, 'k': 1, 'l': 1} format
                    return (miller_indices.get('h', 0), miller_indices.get('k', 0), miller_indices.get('l', 0))
                elif isinstance(miller_indices, (list, tuple)) and len(miller_indices) == 3:
                    # Case: (1, 1, 1) or [1, 1, 1] format
                    return tuple(miller_indices)
                elif isinstance(miller_indices, str):
                    # Case: "(1, 1, 1)" string format
                    import ast
                    parsed_indices = ast.literal_eval(miller_indices)
                    return tuple(parsed_indices)
            
            print("❌ Cannot find Miller indices in Phase3 temporary file")
            return None
            
        except Exception as e:
            print(f"❌ Failed to get Miller indices information: {e}")
            return None
    
    def get_original_cell_vectors(self):
        """Gets original cell vector information."""
        try:
            # Load cell vector information from Phase1 temporary file
            from ..utils.temporary_file_handler import TemporaryFileHandler
            temp_handler = TemporaryFileHandler()
            phase1_data = temp_handler.load_phase_results("phase1")
            
            if phase1_data and 'cell_vectors' in phase1_data:
                cell_vectors_str = phase1_data['cell_vectors']
                
                # Parse vectors stored as strings
                import ast
                cell_vectors = {}
                for key, value_str in cell_vectors_str.items():
                    cell_vectors[key] = ast.literal_eval(value_str)
                
                return cell_vectors
            
            return None
        except Exception as e:
            print(f"⚠️ Failed to get original cell vector information: {e}")
            return None
    
    def get_reference_atom_position(self):
        """Gets the coordinates of the reference atom used in a_user, b_user calculations."""
        try:
            # Get atoms in the first reference plane
            plane_atoms = self.get_first_reference_plane_atoms()
            if not plane_atoms:
                return None
            
            # Calculate plane center
            plane_center = self.calculate_plane_center(plane_atoms)
            
            # Find atom closest to center (same logic used in a_user, b_user)
            reference_atom = self.find_reference_atom(plane_center, plane_atoms)
            
            return reference_atom
        except Exception as e:
            print(f"⚠️ Failed to get reference atom coordinates: {e}")
            return None
    
    def find_next_equivalent_plane_in_direction(self, reference_atom, direction, equivalent_planes):
        """
        Finds the distance to the next equivalent plane by moving in the Miller index direction from the new origin.
        
        Args:
            reference_atom: New origin coordinates (on first reference plane)
            direction: Miller index direction (normalized vector)
            equivalent_planes: List of equivalent plane information
            
        Returns:
            float: Distance to next equivalent plane or None
        """
        try:
            if not equivalent_planes:
                print("⚠️ No equivalent plane information")
                return None
            
            # Find the first equivalent plane encountered while moving in the Miller index direction
            positive_distances = []
            
            for plane_info in equivalent_planes:
                # Representative point of the plane (first atom position)
                plane_point = np.array([
                    plane_info['cart_x'],
                    plane_info['cart_y'], 
                    plane_info['cart_z']
                ])
                
                # Vector from new origin to plane point
                vector_to_plane = plane_point - reference_atom
                
                # Projected distance in the Miller index direction
                projected_distance = np.dot(vector_to_plane, direction)
                
                # Only consider positive direction (Miller index direction)
                if projected_distance > 0.1:  # At least 0.1Å
                    positive_distances.append({
                        'distance': projected_distance,
                        'plane_id': plane_info['plane_id'],
                        'plane_point': plane_point
                    })
                    print(f"   📏 {plane_info['plane_id']}: Distance in Miller index direction {projected_distance:.4f} Å")
            
            if not positive_distances:
                print("⚠️ Cannot find equivalent plane in Miller index direction")
                return None
            
            # Sort by distance and select the closest plane
            positive_distances.sort(key=lambda x: x['distance'])
            next_plane = positive_distances[0]
            
            print(f"📍 Selected next equivalent plane: {next_plane['plane_id']}")
            print(f"📏 Distance in Miller index direction: {next_plane['distance']:.4f} Å")
            
            return next_plane['distance']
            
        except Exception as e:
            print(f"❌ Failed to calculate distance to equivalent plane in Miller index direction: {e}")
            return None
    
    def get_equivalent_planes(self):
        """Gets information about planes that passed the equivalence test."""
        try:
            # Load first reference plane information from temporary file
            from ..utils.temporary_file_handler import TemporaryFileHandler
            temp_handler = TemporaryFileHandler()
            phase3_data = temp_handler.load_phase_results("phase3")
            
            if not phase3_data or 'first_reference_plane' not in phase3_data:
                print("❌ Cannot find first reference plane information")
                return []
            
            first_plane_id = phase3_data['first_reference_plane']['plane_id']
            print(f"📍 First reference plane ID: {first_plane_id}")
            
            import csv
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            
            equivalent_planes = []
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                processed_planes = set()  # Prevent duplicates
                
                for row in reader:
                    plane_id = row.get('plane_id')
                    # Find Type A planes different from the first reference plane
                    if (plane_id and 
                        plane_id != first_plane_id and
                        plane_id not in processed_planes and
                        row.get('plane_type') == 'Type A'):
                        
                        plane_info = {
                            'plane_id': plane_id,
                            'reference_plane': row.get('reference_plane', ''),
                            'plane_type': row.get('plane_type', ''),
                            'plane_equation': row.get('plane_equation', ''),
                            'cart_x': float(row.get('cart_x', 0)),
                            'cart_y': float(row.get('cart_y', 0)),
                            'cart_z': float(row.get('cart_z', 0))
                        }
                        equivalent_planes.append(plane_info)
                        processed_planes.add(plane_id)
            
            print(f"📊 Number of equivalent planes: {len(processed_planes)}")
            for plane_id in sorted(processed_planes):
                print(f"   - {plane_id}")
            
            return equivalent_planes
            
        except Exception as e:
            print(f"❌ Failed to load equivalent plane information: {str(e)}")
            return []
    
    def calculate_distance_to_plane(self, reference_atom, direction, plane_info):
        """Calculates the Miller index direction distance from the reference point to a specific plane."""
        try:
            plane_point = np.array([
                plane_info['cart_x'],
                plane_info['cart_y'], 
                plane_info['cart_z']
            ])
            
            # Vector from reference point to plane point
            vector_to_plane = plane_point - reference_atom
            
            # Projected distance in the Miller index direction
            distance = np.dot(vector_to_plane, direction)
            
            return distance if distance > 0 else None
            
        except Exception as e:
            print(f"⚠️ Failed to calculate distance to plane: {e}")
            return None
    
    def verify_c_user_reaches_equivalent_plane(self, reference_atom, c_user, equivalent_planes):
        """Verifies if the c_user vector actually reaches the equivalent plane."""
        try:
            target_position = reference_atom + c_user
            print(f"🎯 Position moved by c_user: ({target_position[0]:.4f}, {target_position[1]:.4f}, {target_position[2]:.4f})")
            
            # Compare with all atoms in equivalent planes
            import csv
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            
            # Collect equivalent plane IDs
            equivalent_plane_ids = set()
            for plane_info in equivalent_planes:
                equivalent_plane_ids.add(plane_info['plane_id'])
            
            # Check all atoms in equivalent planes from CSV
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                
                for row in reader:
                    plane_id = row.get('plane_id')
                    if plane_id in equivalent_plane_ids:
                        atom_position = np.array([
                            float(row.get('cart_x', 0)),
                            float(row.get('cart_y', 0)),
                            float(row.get('cart_z', 0))
                        ])
                        
                        distance = np.linalg.norm(target_position - atom_position)
                        if distance < 0.1:  # Consider near equivalent plane if within 0.1Å
                            atom_name = row.get('atom_name', '')
                            print(f"✅ Confirmed reaching near equivalent plane - {atom_name}: ({atom_position[0]:.6f}, {atom_position[1]:.6f}, {atom_position[2]:.6f})")
                            print(f"   Distance: {distance:.6f} Å")
                            return True
            
            print("⚠️ Failed to reach near equivalent plane")
            return False
            
        except Exception as e:
            print(f"⚠️ c_user verification failed: {e}")
            return False
            
    def display_final_cell_vectors(self, a_user, b_user, c_user, sequence_info, transformation_matrix):
        """
        Displays final cell vector information.
        
        Args:
            a_user: a_user vector
            b_user: b_user vector
            c_user: c_user vector
            sequence_info: Sequence analysis information
            transformation_matrix: Transformation matrix
        """
        print("\n📊 Step 3-3: Final cell vector summary")
        print("=" * 50)
        
        print(f"📐 User-defined cell vectors:")
        print(f"   a_user = ({a_user[0]:.4f}, {a_user[1]:.4f}, {a_user[2]:.4f}), magnitude: {np.linalg.norm(a_user):.4f} Å")
        print(f"   b_user = ({b_user[0]:.4f}, {b_user[1]:.4f}, {b_user[2]:.4f}), magnitude: {np.linalg.norm(b_user):.4f} Å")
        print(f"   c_user = ({c_user[0]:.4f}, {c_user[1]:.4f}, {c_user[2]:.4f}), magnitude: {np.linalg.norm(c_user):.4f} Å")
        
        # Calculate cell volume
        cell_volume = abs(np.dot(a_user, np.cross(b_user, c_user)))
        print(f"   Cell volume: {cell_volume:.4f} Å³")
        
        # Calculate vector angles
        def calculate_angle(v1, v2):
            cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            # Handle numerical errors
            cos_angle = min(1.0, max(-1.0, cos_angle))
            return np.degrees(np.arccos(cos_angle))
        
        alpha = calculate_angle(b_user, c_user)
        beta = calculate_angle(a_user, c_user)
        gamma = calculate_angle(a_user, b_user)
        
        print(f"   Cell angles: α={alpha:.2f}°, β={beta:.2f}°, γ={gamma:.2f}°")
        
        # Summary of sequence information
        if sequence_info:
            if 'is_default' in sequence_info:
                print(f"\n⚠️ Default pattern used: {sequence_info['pattern']}")
            else:
                pattern = sequence_info.get('pattern', [])
                pattern_length = sequence_info.get('pattern_length', 0)
                repetitions = sequence_info.get('repetitions', 0)
                coverage = sequence_info.get('coverage', 0)
                
                print(f"\n📋 Sequence information:")
                print(f"   Pattern: {pattern}")
                print(f"   Pattern length: {pattern_length}")
                print(f"   Repetitions: {repetitions}")
                print(f"   Coverage: {coverage:.2f}")
        
        print("\n✅ Cell vector calculation completed") 