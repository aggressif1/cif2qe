import numpy as np
import math
from .periodicity_calculator import PeriodicityCalculator
import csv
import re
import os

class PlaneComparisonAnalyzer:
    """Class for performing plane comparison analysis"""
    
    def __init__(self):
        """Initialize PlaneComparisonAnalyzer"""
        self.bond_criteria = None  # Bond criteria received from Step2
        self.periodicity_calc = None  # Periodicity calculator
        
    def set_bond_criteria(self, bond_criteria):
        """
        Set bond criteria received from Step2
        
        Args:
            bond_criteria (dict): Bond type, minimum/maximum distance information
        """
        self.bond_criteria = bond_criteria
        print(f"Bond analysis criteria set: {bond_criteria}")
    
    def set_cell_vectors(self, cell_vectors):
        """
        Set cell vector information to enable periodicity calculation.
        
        Args:
            cell_vectors (dict): Cell vector information {'a': [x,y,z], 'b': [x,y,z], 'c': [x,y,z]}
        """
        self.cell_vectors = cell_vectors  # Store cell vectors directly
        self.periodicity_calc = PeriodicityCalculator(cell_vectors)
        print(f"Periodicity calculator setup complete")
    
    def calculate_plane_center(self, plane_atoms):
        """
        Calculate the center point of the comparison plane.
        
        Args:
            plane_atoms (list): List of atoms belonging to the plane
            
        Returns:
            tuple: Center point coordinates (x, y, z)
        """
        if not plane_atoms:
            print("No atoms in the plane.")
            return None
        
        try:
            # Calculate average coordinates of all atoms
            total_x = sum(atom['cart_x'] for atom in plane_atoms)
            total_y = sum(atom['cart_y'] for atom in plane_atoms)
            total_z = sum(atom['cart_z'] for atom in plane_atoms)
            
            count = len(plane_atoms)
            center_x = total_x / count
            center_y = total_y / count
            center_z = total_z / count
            
            center = (center_x, center_y, center_z)
            
            print(f"Plane center point: ({center_x:.4f}, {center_y:.4f}, {center_z:.4f})")
            print(f"Plane atom count: {count}")
            
            return center
            
        except Exception as e:
            print(f"Plane center point calculation failed: {str(e)}")
            return None
    
    def find_closest_atoms_to_center(self, plane_atoms, center, num_atoms=4):
        """
        Select the atoms closest to the center point.
        
        Args:
            plane_atoms (list): List of atoms belonging to the plane
            center (tuple): Center point coordinates (x, y, z)
            num_atoms (int): Number of atoms to select (default: 4)
            
        Returns:
            list: List of atoms closest to the center point
        """
        if not plane_atoms or not center:
            print("Invalid input.")
            return []
        
        if len(plane_atoms) < num_atoms:
            print(f"Plane atom count ({len(plane_atoms)}) is less than requested atom count ({num_atoms}).")
            num_atoms = len(plane_atoms)
        
        try:
            center_x, center_y, center_z = center
            
            # Calculate distance between each atom and center point
            atom_distances = []
            
            for i, atom in enumerate(plane_atoms):
                distance = math.sqrt(
                    (atom['cart_x'] - center_x) ** 2 +
                    (atom['cart_y'] - center_y) ** 2 +
                    (atom['cart_z'] - center_z) ** 2
                )
                
                atom_distances.append((distance, i, atom))
            
            # Sort by distance
            atom_distances.sort(key=lambda x: x[0])
            
            # Select closest atoms
            closest_atoms = []
            for i in range(num_atoms):
                distance, atom_index, atom = atom_distances[i]
                
                # Add unique ID to atom
                atom_with_id = atom.copy()
                
                # Use existing atom_id if available, otherwise create one
                if 'atom_id' in atom:
                    atom_with_id['atom_id'] = atom['atom_id']
                else:
                    # Use legacy logic (backward compatibility)
                    atom_with_id['atom_id'] = atom.get('label', f"{atom.get('element', 'Unknown')}{i+1}")
                
                atom_with_id['distance_to_center'] = distance
                
                closest_atoms.append(atom_with_id)
            
            print(f"Successfully selected {num_atoms} atoms closest to center:")
            for i, atom in enumerate(closest_atoms):
                element = atom.get('element', 'Unknown')
                x, y, z = atom['cart_x'], atom['cart_y'], atom['cart_z']
                distance = atom['distance_to_center']
                atom_id = atom['atom_id']
                print(f"   {i+1}. {atom_id} ({element}): Position({x:.4f}, {y:.4f}, {z:.4f}), Center distance: {distance:.4f}Å")
            
            return closest_atoms
            
        except Exception as e:
            print(f"Selection of atoms closest to center failed: {str(e)}")
            return []
    
    def analyze_comparison_plane_atoms(self, plane_assignment):
        """
        Analyze atoms of the comparison plane and select the 4 atoms closest to the center.
        
        Args:
            plane_assignment (dict): Plane assignment information
            
        Returns:
            dict: Analysis result {'center': tuple, 'selected_atoms': list} or None
        """
        print("\nStarting comparison plane atom analysis...")
        
        try:
            plane_atoms = plane_assignment.get('atoms', [])
            plane_id = plane_assignment.get('plane_info', {}).get('plane_id', 'Unknown')
            
            print(f"Target plane for analysis: {plane_id}")
            
            if len(plane_atoms) < 4:
                print(f"Insufficient plane atoms (minimum 4 required, current {len(plane_atoms)})")
                return None
            
            # Step 1: Calculate center point
            center = self.calculate_plane_center(plane_atoms)
            if not center:
                return None
            
            # Step 2: Select 4 atoms closest to center
            selected_atoms = self.find_closest_atoms_to_center(plane_atoms, center, 4)
            if len(selected_atoms) < 4:
                print(f"Failed to select sufficient atoms.")
                return None
            
            result = {
                'plane_id': plane_id,
                'center': center,
                'selected_atoms': selected_atoms,
                'total_atoms': len(plane_atoms)
            }
            
            print(f"Comparison plane atom analysis completed")
            return result
            
        except Exception as e:
            print(f"Comparison plane atom analysis failed: {str(e)}")
            return None
    
    def display_analysis_summary(self, analysis_result):
        """
        Display analysis result summary.
        
        Args:
            analysis_result (dict): Analysis result
        """
        if not analysis_result:
            print("No analysis result to display.")
            return
        
        print("\n" + "=" * 60)
        print("Comparison Plane Analysis Result Summary")
        print("=" * 60)
        
        plane_id = analysis_result.get('plane_id', 'Unknown')
        center = analysis_result.get('center')
        selected_atoms = analysis_result.get('selected_atoms', [])
        total_atoms = analysis_result.get('total_atoms', 0)
        
        print(f"Plane ID: {plane_id}")
        print(f"Total atoms: {total_atoms}")
        print(f"Center point: ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f})")
        print(f"Selected atoms: {len(selected_atoms)} | Selected atom list:")
        print(f"{'No.':<4} {'Atom ID':<10} {'Element':<6} {'Coordinates (X, Y, Z)':<25} {'Center Dist(Å)':<12}")
        print("-" * 65)
        
        for i, atom in enumerate(selected_atoms, 1):
            atom_id = atom.get('atom_id', 'Unknown')
            element = atom.get('element', 'Unknown')
            x, y, z = atom['cart_x'], atom['cart_y'], atom['cart_z']
            distance = atom.get('distance_to_center', 0)
            
            coords_str = f"({x:.3f}, {y:.3f}, {z:.3f})"
            print(f"{i:<4} {atom_id:<10} {element:<6} {coords_str:<25} {distance:<12.4f}")
        
        print("=" * 60)
    
    def analyze_bond_information(self, selected_atoms, all_supercell_atoms):
        """
        Analyze bond information of selected atoms.
        
        Args:
            selected_atoms (list): Selected 4 atoms
            all_supercell_atoms (list): All supercell atom information
            
        Returns:
            dict: Bond information for each atom or None
        """
        print("\nStarting bond analysis of selected atoms...")
        
        if not self.bond_criteria:
            print("Bond analysis criteria not set.")
            return None
        
        if not selected_atoms or len(selected_atoms) < 4:
            print("Insufficient selected atoms.")
            return None
        
        try:
            bond_analysis_results = {}
            min_distance = self.bond_criteria['min_distance']
            max_distance = self.bond_criteria['max_distance']
            
            print(f"Bond analysis criteria:")
            print(f"   - Bond type: {self.bond_criteria['bond_type']}")
            print(f"   - Minimum distance: {min_distance}Å")
            print(f"   - Maximum distance: {max_distance}Å")
            
            for i, atom in enumerate(selected_atoms):
                atom_id = atom['atom_id']
                atom_element = atom.get('element', 'Unknown')
                atom_pos = (atom['cart_x'], atom['cart_y'], atom['cart_z'])
                
                print(f"\nAtom {i+1}: {atom_id} ({atom_element})")
                print(f"   Position: ({atom_pos[0]:.4f}, {atom_pos[1]:.4f}, {atom_pos[2]:.4f})")
                
                # Find bond partners for this atom
                bond_partners = self._find_bond_partners(
                    atom, all_supercell_atoms, min_distance, max_distance
                )
                
                if bond_partners:
                    print(f"   Bonds found: {len(bond_partners)}")
                    
                    bond_info = {
                        'atom_id': atom_id,
                        'atom_element': atom_element,
                        'atom_position': atom_pos,
                        'bond_count': len(bond_partners),
                        'bond_partners': bond_partners
                    }
                    
                    # Display detailed bond information
                    for j, partner in enumerate(bond_partners):
                        partner_element = partner['element']
                        distance = partner['distance']
                        vector = partner['bond_vector']
                        print(f"      {j+1}. {partner_element}: {distance:.4f}Å, Vector({vector[0]:.3f}, {vector[1]:.3f}, {vector[2]:.3f})")
                    
                    bond_analysis_results[atom_id] = bond_info
                else:
                    print("   No bond partners found.")
                    bond_analysis_results[atom_id] = {
                        'atom_id': atom_id,
                        'atom_element': atom_element,
                        'atom_position': atom_pos,
                        'bond_count': 0,
                        'bond_partners': []
                    }
            
            print(f"\nBond analysis completed: {len(bond_analysis_results)} atoms processed")
            return bond_analysis_results
            
        except Exception as e:
            print(f"Bond analysis failed: {str(e)}")
            return None
    
    def _find_bond_partners(self, target_atom, all_atoms, min_distance, max_distance):
        """
        Find bond partners of target atom. (Apply periodicity + remove duplicates)
        
        Args:
            target_atom (dict): Target atom information
            all_atoms (list): All atom list
            min_distance (float): Minimum bond distance
            max_distance (float): Maximum bond distance
            
        Returns:
            list: List of bond partners (duplicates removed)
        """
        target_pos = (target_atom['cart_x'], target_atom['cart_y'], target_atom['cart_z'])
        target_element = target_atom.get('element', 'Unknown')
        
        bond_partners = []
        seen_bonds = set()  # Set for duplicate removal
        
        for atom in all_atoms:
            # Exclude self
            atom_pos = (atom['cart_x'], atom['cart_y'], atom['cart_z'])
            if abs(atom_pos[0] - target_pos[0]) < 1e-6 and \
               abs(atom_pos[1] - target_pos[1]) < 1e-6 and \
               abs(atom_pos[2] - target_pos[2]) < 1e-6:
                continue
            
            # Calculate distance considering periodicity
            if self.periodicity_calc:
                # Calculate minimum image distance applying periodic boundary conditions
                # Vector from central atom (target_atom) to bonding atom (atom)
                distance, bond_vector = self.periodicity_calc.find_minimum_image_distance(target_atom, atom)
            else:
                # Regular distance calculation (no periodicity)
                # Vector from central atom (target_atom) to bonding atom (atom)
                dx = atom['cart_x'] - target_atom['cart_x']
                dy = atom['cart_y'] - target_atom['cart_y']
                dz = atom['cart_z'] - target_atom['cart_z']
                distance = math.sqrt(dx*dx + dy*dy + dz*dz)
                bond_vector = (dx, dy, dz)
            
            # Check if within bond distance range
            if min_distance <= distance <= max_distance:
                # Generate unique key for duplicate removal (using rounded distance + element + vector)
                rounded_distance = round(distance, 4)
                rounded_vector = (round(bond_vector[0], 3), round(bond_vector[1], 3), round(bond_vector[2], 3))
                bond_key = (atom.get('element', 'Unknown'), rounded_distance, rounded_vector)
                
                # Check if bond already seen
                if bond_key not in seen_bonds:
                    seen_bonds.add(bond_key)
                    
                    partner_info = {
                        'element': atom.get('element', 'Unknown'),
                        'atom_id': atom.get('atom_id', 'Unknown'),
                        'position': atom_pos,
                        'distance': distance,
                        'bond_vector': bond_vector,
                        'periodic': self.periodicity_calc is not None,
                        'partner_atom': atom  # Preserve complete atom information
                    }
                    
                    bond_partners.append(partner_info)
        
        # Sort by distance
        bond_partners.sort(key=lambda x: x['distance'])
        
        return bond_partners
    
    def display_bond_analysis_summary(self, bond_analysis_results):
        """
        Display bond analysis result summary.
        
        Args:
            bond_analysis_results (dict): Bond analysis results
        """
        if not bond_analysis_results:
            print("No bond analysis results to display.")
            return
        
        print("\n" + "=" * 80)
        print("Bond Analysis Result Summary")
        print("=" * 80)
        
        total_bonds = 0
        
        print(f"{'Atom ID':<12} {'Element':<6} {'Bond Count':<8} {'Bond Partners':<40}")
        print("-" * 80)
        
        for atom_id, bond_info in bond_analysis_results.items():
            element = bond_info['atom_element']
            bond_count = bond_info['bond_count']
            total_bonds += bond_count
            
            # Generate bond partner summary
            partners_summary = []
            for partner in bond_info['bond_partners']:
                partner_element = partner['element']
                distance = partner['distance']
                partners_summary.append(f"{partner_element}({distance:.3f}Å)")
            
            partners_str = ", ".join(partners_summary[:3])  # Display maximum 3
            if len(bond_info['bond_partners']) > 3:
                partners_str += f" +{len(bond_info['bond_partners']) - 3} more"
            
            print(f"{atom_id:<12} {element:<6} {bond_count:<8} {partners_str}")
        
        print("-" * 80)
        print(f"Total bonds: {total_bonds}")
        print("=" * 80)
    
    def verify_atoms_on_plane(self, selected_atoms, plane_info, tolerance=0.01):
        """
        Verify that selected atoms are actually on the corresponding plane.
        
        Args:
            selected_atoms (list): Selected atoms
            plane_info (dict): Plane information (including plane equation)
            tolerance (float): Tolerance error
            
        Returns:
            dict: Verification results
        """
        print("\nStarting plane position verification...")
        
        if not selected_atoms or not plane_info:
            print("No input data.")
            return None
        
        try:
            # Extract plane equation information
            plane_equation = plane_info.get('plane_equation', '')
            coefficients = plane_info.get('coefficients', {})
            
            if not coefficients:
                print("No plane equation coefficient information.")
                return None
            
            a = coefficients.get('a', 0)
            b = coefficients.get('b', 0) 
            c = coefficients.get('c', 0)
            d = coefficients.get('d', 0)
            
            print(f"Plane equation: Coefficients: a={a:.4f}, b={b:.4f}, c={c:.4f}, d={d:.4f} | Tolerance: ±{tolerance}")
            
            verification_results = []
            all_on_plane = True
            
            for i, atom in enumerate(selected_atoms):
                atom_id = atom.get('atom_id', f'Atom_{i+1}')
                element = atom.get('element', 'Unknown')
                x = atom['cart_x']
                y = atom['cart_y'] 
                z = atom['cart_z']
                
                # Substitute coordinates into plane equation: ax + by + cz + d = 0
                plane_value = a * x + b * y + c * z + d
                distance_to_plane = abs(plane_value)
                
                is_on_plane = distance_to_plane <= tolerance
                
                verification_results.append({
                    'atom_id': atom_id,
                    'element': element,
                    'coordinates': (x, y, z),
                    'plane_equation_value': plane_value,
                    'distance_to_plane': distance_to_plane,
                    'is_on_plane': is_on_plane
                })
                
                if not is_on_plane:
                    all_on_plane = False
            
            # Display plane position verification results in table format
            print(f"\nPlane Position Verification Results:")
            print("=" * 110)
            print(f"{'No.':<4} {'Atom ID':<12} {'Element':<4} {'Coordinates (X, Y, Z)':<30} {'Plane Eq Value':<12} {'Plane Dist':<10} {'On Plane':<10}")
            print("-" * 110)
            
            for i, result in enumerate(verification_results, 1):
                atom_id = result['atom_id']
                element = result['element']
                x, y, z = result['coordinates']
                plane_val = result['plane_equation_value']
                distance = result['distance_to_plane']
                on_plane = result['is_on_plane']
                
                coord_str = f"({x:.3f}, {y:.3f}, {z:.3f})"
                plane_val_str = f"{plane_val:.6f}"
                distance_str = f"{distance:.6f}Å"
                status_icon = "Yes" if on_plane else "No"
                
                print(f"{i:<4} {atom_id:<12} {element:<4} {coord_str:<30} {plane_val_str:<12} {distance_str:<10} {status_icon:<10}")
            
            print("=" * 110)
            
            result = {
                'plane_equation': plane_equation,
                'tolerance': tolerance,
                'total_atoms': len(selected_atoms),
                'atoms_on_plane': sum(1 for r in verification_results if r['is_on_plane']),
                'all_on_plane': all_on_plane,
                'verification_details': verification_results
            }
            
            print(f"\nVerification Result Summary:")
            print(f"   - Total atoms: {result['total_atoms']}")
            print(f"   - Atoms on plane: {result['atoms_on_plane']}")
            print(f"   - All on plane: {'Yes' if all_on_plane else 'No'}")
            
            if not all_on_plane:
                print("Some atoms are off the plane. Please check the plane assignment logic.")
            else:
                print("All atoms are on the correct plane.")
            
            return result
            
        except Exception as e:
            print(f"Plane position verification failed: {str(e)}")
            return None
    
    def clone_selected_atoms(self, selected_atoms):
        """
        Clone selected atoms and group them.
        Cloned atoms have 'A' appended to their atom index.
        
        Args:
            selected_atoms (list): Selected 4 atoms
            
        Returns:
            dict: Cloned atom group information or None
        """
        print("\nStarting cloning and grouping of selected atoms...")
        
        if not selected_atoms or len(selected_atoms) < 4:
            print("Insufficient atoms to clone.")
            return None
        
        try:
            cloned_atoms = []
            clone_list = []
            
            for i, atom in enumerate(selected_atoms):
                # Copy original atom information
                cloned_atom = atom.copy()
                
                # Generate new ID (original ID + 'A' suffix)
                original_id = atom['atom_id']
                cloned_id = f"{original_id}A"
                cloned_atom['atom_id'] = cloned_id
                cloned_atom['original_atom_id'] = original_id
                
                cloned_atoms.append(cloned_atom)
                clone_list.append(f"{original_id}→{cloned_id}")
            
            print(f"   Cloning completed: {' | '.join(clone_list)}")
            
            # Generate group information
            group_info = {
                'group_id': f"ClonedGroup_{len(selected_atoms)}atoms",
                'original_atoms': selected_atoms,
                'cloned_atoms': cloned_atoms,
                'atom_count': len(cloned_atoms),
                'creation_timestamp': None  # Timestamp can be added if needed
            }
            
            print(f"Atom cloning and grouping completed: Group ID: {group_info['group_id']}")
            
            return group_info
            
        except Exception as e:
            print(f"Atom cloning and grouping failed: {str(e)}")
            return None
    
    def move_cloned_atoms_to_reference_plane(self, cloned_group, reference_plane_info, miller_indices):
        """
        Move cloned atoms to the reference plane.
        The direction of the movement vector is calculated using cell vectors and Miller indices.
        
        Args:
            cloned_group (dict): Cloned atom group information
            reference_plane_info (dict): Reference plane information
            miller_indices (tuple): Miller indices (h, k, l)
            
        Returns:
            dict: Moved atom group information or None
        """
        print("\nStarting movement of cloned atoms to reference plane...")
        
        if not cloned_group or not reference_plane_info or not miller_indices:
            print("Insufficient input data.")
            return None
        
        if not self.cell_vectors:
            print("Cell vector information not set.")
            return None
        
        try:
            cloned_atoms = cloned_group['cloned_atoms']
            
            # Process Miller indices
            if isinstance(miller_indices, (list, tuple)) and len(miller_indices) == 3:
                h, k, l = miller_indices
                print(f"Miller indices: ({h} {k} {l})")
            else:
                print(f"Miller indices format error: {miller_indices}")
                return None
            # Display only simple summary information (excluding large atom information)
            summary_info = {
                'plane_id': reference_plane_info.get('plane_id', 'Unknown'),
                'd_value': reference_plane_info.get('d_value', 0.0),  # Set default to 0.0
                'atom_count': reference_plane_info.get('atom_count', 'Unknown')
            }
            
            print(f"Reference plane info: Plane ID: {summary_info['plane_id']} | D-value: {float(summary_info['d_value']):.4f} | Atom count: {summary_info['atom_count']}")
            
            # Extract reference plane equation coefficients - Read normalized plane equation from CSV file
            print(f"Reading normalized equation of reference plane from CSV file...")
            
            # Find reference plane information in CSV file
            import csv
            import re
            
            csv_file_path = "output/supercell_coordinates.csv"
            
            # Extract plane_id accurately
            raw_plane_id = reference_plane_info.get('plane_id', 'Unknown')
            
            # CSV stores in 'Plane 1', 'Plane 2' format
            if isinstance(raw_plane_id, str):
                if raw_plane_id.startswith('Plane '):
                    # Use as-is if already in 'Plane X' format
                    target_plane_id = raw_plane_id
                elif raw_plane_id.isdigit():
                    # Convert to 'Plane X' format if only digits
                    target_plane_id = f"Plane {raw_plane_id}"
                else:
                    # Use as-is for other cases
                    target_plane_id = raw_plane_id
            else:
                # Convert to 'Plane X' format if numeric type
                target_plane_id = f"Plane {raw_plane_id}"
            
            print(f"CSV file path: {csv_file_path}")
            
            plane_equation_text = None
            found_plane_ids = []  # For debugging
            
            try:
                with open(csv_file_path, 'r', encoding='utf-8') as file:
                    csv_reader = csv.DictReader(file)
                    
                    for row in csv_reader:
                        current_plane_id = row['plane_id']
                        found_plane_ids.append(current_plane_id)
                        
                        if current_plane_id == target_plane_id:
                            plane_equation_text = row['plane_equation']
                            break
                
                if not plane_equation_text:
                    print(f"Could not find plane equation for {target_plane_id} in CSV file.")
                    print(f"Found plane_ids in CSV (first 10): {found_plane_ids[:10]}")
                    print(f"Total unique plane_ids found: {len(set(found_plane_ids))}")
                    
                    # Find similar IDs
                    similar_ids = [pid for pid in set(found_plane_ids) if target_plane_id in pid or pid in target_plane_id]
                    if similar_ids:
                        print(f"Similar plane_ids: {similar_ids}")
                    
                    return None
                
                # Parse plane equation: "0.7071x + 0.7071y - 23.7084 = 0" format
                # First check if z term exists
                if 'z' in plane_equation_text:
                    # Case with z term: ax + by + cz + d = 0
                    pattern = r'([+-]?\d*\.?\d+)x\s*([+-]\s*\d*\.?\d+)y\s*([+-]\s*\d*\.?\d+)z\s*([+-]\s*\d*\.?\d+)\s*=\s*0'
                    match = re.search(pattern, plane_equation_text)
                    if match:
                        plane_h = float(match.group(1))
                        plane_k = float(match.group(2).replace(' ', ''))
                        plane_l = float(match.group(3).replace(' ', ''))
                        plane_d = float(match.group(4).replace(' ', ''))
                    else:
                        print(f"Failed to parse plane equation with z term: {plane_equation_text}")
                        return None
                else:
                    # Case without z term - use more general parsing method
                    print(f"Parsing plane equation without z term: {plane_equation_text}")
                    
                    # Extract coefficients by finding each term individually
                    plane_h = 0.0  # x coefficient
                    plane_k = 0.0  # y coefficient 
                    plane_l = 0.0  # z coefficient (always 0)
                    plane_d = 0.0  # constant term
                    
                    # Find x term
                    x_pattern = r'([+-]?\d*\.?\d*)x'
                    x_match = re.search(x_pattern, plane_equation_text)
                    if x_match:
                        x_coeff = x_match.group(1)
                        if x_coeff == '' or x_coeff == '+':
                            plane_h = 1.0
                        elif x_coeff == '-':
                            plane_h = -1.0
                        else:
                            plane_h = float(x_coeff)
                    
                    # Find y term
                    y_pattern = r'([+-]?\d*\.?\d*)y'
                    y_match = re.search(y_pattern, plane_equation_text)
                    if y_match:
                        y_coeff = y_match.group(1)
                        if y_coeff == '' or y_coeff == '+':
                            plane_k = 1.0
                        elif y_coeff == '-':
                            plane_k = -1.0
                        else:
                            plane_k = float(y_coeff)
                    
                    # Find constant term - use direct method
                    # Find constant term on the left side of equation
                    equation_left = plane_equation_text.split('=')[0].strip()
                    
                    # Find terms consisting only of numbers (no variables)
                    # Use more precise pattern to distinguish from variable coefficients
                    # Find only independent constant terms not immediately before variables (x,y,z)
                    const_pattern = r'([+-]\s*\d+\.?\d*)\s*(?![xyz]\w*)'
                    
                    # Additionally, remove variable coefficient parts first, then find constant terms
                    temp_equation = equation_left
                    # Remove x, y, z coefficient parts (e.g., remove "0.7071x", "0.7071y")
                    temp_equation = re.sub(r'[+-]?\s*\d*\.?\d*[xyz]', '', temp_equation)
                    temp_equation = temp_equation.strip()
                    
                    # Find constant terms in remaining part
                    const_matches = re.findall(r'([+-]?\s*\d+\.?\d*)', temp_equation)
                    
                    for const in const_matches:
                        const = const.replace(' ', '').strip()  # Remove spaces
                        if const:
                            try:
                                plane_d += float(const)
                                print(f"   - Constant term found: {const}")
                            except ValueError:
                                pass
                    
                    print(f"Parsing results:")
                    print(f"   - x coefficient (h): {plane_h}")
                    print(f"   - y coefficient (k): {plane_k}")
                    print(f"   - z coefficient (l): {plane_l}")
                    print(f"   - constant term (d): {plane_d}")
                    
                    # Check if parsing was successful (at least one coefficient should not be 0)
                    if abs(plane_h) < 1e-10 and abs(plane_k) < 1e-10 and abs(plane_l) < 1e-10:
                        print(f"Failed to parse plane equation without z term: all coefficients are 0")
                        return None
                    
                    print(f"Successfully parsed plane equation without z term")
                
                print(f"Plane found in CSV: {target_plane_id} → {plane_equation_text}")
                
            except Exception as e:
                print(f"Failed to read CSV file: {str(e)}")
                return None
            
            # Calculate movement vector using cell vectors and Miller indices
            # Movement vector = h*a + k*b + l*c (user requirement)
            # cell_vectors is dict format: {'a': [x,y,z], 'b': [x,y,z], 'c': [x,y,z]}
            a_vec = self.cell_vectors['a']  # [5.5881, 0, 0]
            b_vec = self.cell_vectors['b']  # [0, 5.5881, 0] 
            c_vec = self.cell_vectors['c']  # [0, 0, 5.5881]
            
            # Convert vectors to lists for safe processing
            a_vec = [float(a_vec[0]), float(a_vec[1]), float(a_vec[2])]
            b_vec = [float(b_vec[0]), float(b_vec[1]), float(b_vec[2])]
            c_vec = [float(c_vec[0]), float(c_vec[1]), float(c_vec[2])]
            
            # Convert Miller indices to float
            h, k, l = float(h), float(k), float(l)
            
            # User requirement: direction of movement vector = h*a + k*b + l*c
            direction_vector = [
                h * a_vec[0] + k * b_vec[0] + l * c_vec[0],
                h * a_vec[1] + k * b_vec[1] + l * c_vec[1], 
                h * a_vec[2] + k * b_vec[2] + l * c_vec[2]
            ]
            
            print(f"Movement vector information:")
            print(f"   - Cell vectors: a=({a_vec[0]:.4f}, {a_vec[1]:.4f}, {a_vec[2]:.4f}) | b=({b_vec[0]:.4f}, {b_vec[1]:.4f}, {b_vec[2]:.4f}) | c=({c_vec[0]:.4f}, {c_vec[1]:.4f}, {c_vec[2]:.4f})")
            print(f"   - Miller indices: ({h:.0f}, {k:.0f}, {l:.0f}) | Calculated direction vector (h*a + k*b + l*c): ({direction_vector[0]:.4f}, {direction_vector[1]:.4f}, {direction_vector[2]:.4f})")
            
            # User requirement: precise movement vector calculation method
            print(f"Precise movement calculation:")
            
            # 1. Calculate center point of current plane (comparison plane)
            current_center = [0, 0, 0]
            if len(cloned_atoms) > 0:
                for atom in cloned_atoms:
                    current_center[0] += atom['cart_x']
                    current_center[1] += atom['cart_y'] 
                    current_center[2] += atom['cart_z']
                current_center = [coord / len(cloned_atoms) for coord in current_center]
            
            print(f"   1. Comparison plane center: ({current_center[0]:.4f}, {current_center[1]:.4f}, {current_center[2]:.4f})")
            
            # 2. Normalize direction vector (make unit vector)
            import math
            direction_magnitude = math.sqrt(sum(comp**2 for comp in direction_vector))
            if direction_magnitude == 0:
                print("Movement direction vector is 0.")
                return None
            
            unit_direction = [comp / direction_magnitude for comp in direction_vector]
            print(f"   2. Normalized direction vector: ({unit_direction[0]:.4f}, {unit_direction[1]:.4f}, {unit_direction[2]:.4f})")
            
            # 3. Get reference plane equation coefficients
            # Use normalized coefficients read from CSV
            print(f"   3. Normalized reference plane equation: {plane_h:.4f}x + {plane_k:.4f}y + {plane_l:.4f}z + ({plane_d:.4f}) = 0")
            
            # 4. Move comparison plane center point in direction vector direction to find point that satisfies reference plane
            # Use parameter t: new_point = current_center + t * unit_direction_vector
            # Substitute into reference plane equation: h(x0 + t*dx) + k(y0 + t*dy) + l(z0 + t*dz) + d = 0
            # Rearranged: (hx0 + ky0 + lz0 + d) + t(h*dx + k*dy + l*dz) = 0
            # Therefore: t = -(hx0 + ky0 + lz0 + d) / (h*dx + k*dy + l*dz)
            
            x0, y0, z0 = current_center
            dx, dy, dz = unit_direction
            
            # Numerator: value when current center point is substituted into reference plane equation
            numerator = plane_h * x0 + plane_k * y0 + plane_l * z0 + plane_d
            
            # Denominator: dot product of direction vector and plane normal vector
            denominator = plane_h * dx + plane_k * dy + plane_l * dz
            
            print(f"   4. Parameter calculation:")
            print(f"      - Numerator (current point plane equation value): {numerator:.4f}")
            print(f"      - Denominator (direction vector · plane normal): {denominator:.4f}")
            
            if abs(denominator) < 1e-10:
                print("Direction vector and plane are parallel. Cannot move.")
                return None
            
            # Calculate parameter t
            t = -numerator / denominator
            print(f"      - Parameter t: {t:.4f}")
            
            # 5. Calculate new center point
            new_center = [
                x0 + t * dx,
                y0 + t * dy,
                z0 + t * dz
            ]
            
            # 6. Verify reference plane equation
            verification_value = plane_h * new_center[0] + plane_k * new_center[1] + plane_l * new_center[2] + plane_d
            print(f"   5. New center point: ({new_center[0]:.4f}, {new_center[1]:.4f}, {new_center[2]:.4f})")
            print(f"   6. Reference plane equation verification: {verification_value:.8f} (should be close to 0)")
            
            # 7. Calculate final movement vector (new center point - current center point)
            move_vector = [
                new_center[0] - current_center[0],
                new_center[1] - current_center[1],
                new_center[2] - current_center[2]
            ]
            
            # 8. Calculate movement distance
            move_distance = math.sqrt(sum(comp**2 for comp in move_vector))
            
            print(f"   7. Final movement vector: ({move_vector[0]:.4f}, {move_vector[1]:.4f}, {move_vector[2]:.4f})")
            print(f"   8. Movement distance: {move_distance:.4f}Å")
            
            # 9. Move all cloned atoms by the movement vector
            print(f"   9. Moving cloned atoms:")
            for atom in cloned_atoms:
                old_x, old_y, old_z = atom['cart_x'], atom['cart_y'], atom['cart_z']
                new_x = old_x + move_vector[0]
                new_y = old_y + move_vector[1]
                new_z = old_z + move_vector[2]
                
                atom['cart_x'] = new_x
                atom['cart_y'] = new_y
                atom['cart_z'] = new_z
                
                atom_id = atom.get('atom_id', 'Unknown')
                element = atom.get('element', 'Unknown')
                print(f"   {atom_id} ({element}): ({old_x:.4f}, {old_y:.4f}, {old_z:.4f}) → ({new_x:.4f}, {new_y:.4f}, {new_z:.4f})")
            
            # Generate moved group information
            moved_group = cloned_group.copy()
            moved_group['moved_atoms'] = cloned_atoms
            moved_group['move_vector'] = move_vector
            moved_group['target_plane'] = reference_plane_info.get('plane_id', 'Unknown')
            moved_group['miller_indices'] = miller_indices
            
            return moved_group
            
        except Exception as e:
            print(f"Atom movement failed: {str(e)}")
            return None
    
    def perform_equivalence_test_step1(self, moved_group, reference_plane_atoms, tolerance=0.01):
        """
        Equivalence test step 1: Check distances between cloned atoms and reference plane atoms
        
        Args:
            moved_group (dict): Moved atom group information
            reference_plane_atoms (list): Atoms of the reference plane
            tolerance (float): Allowed tolerance (default: 0.01Å)
            
        Returns:
            dict: Test results or None
        """

        
        if not moved_group or not reference_plane_atoms:
            print("Insufficient input data.")
            return None
        
        try:
            moved_atoms = moved_group.get('moved_atoms', [])
            target_plane = moved_group.get('target_plane', 'Unknown')
            
            print(f"Reference plane: {target_plane}")
            print(f"Number of moved atoms: {len(moved_atoms)}")
            print(f"Number of reference plane atoms: {len(reference_plane_atoms)}")
            print(f"Tolerance: ±{tolerance}Å")
            
            distance_results = []
            all_passed = True
            
            # Pre-read normalized normal vector of reference plane from CSV file
            import csv
            import re
            
            csv_file_path = "output/supercell_coordinates.csv"
            target_plane_id = moved_group.get('target_plane', 'Unknown')
            
            reference_normal = (1.0, 0.0, 0.0)  # Default value (Miller (1,0,0))
            
            try:
                with open(csv_file_path, 'r', encoding='utf-8') as file:
                    csv_reader = csv.DictReader(file)
                    for row in csv_reader:
                        if row['plane_id'] == target_plane_id:
                            plane_equation_text = row['plane_equation']
                            
                            # Extract normal vector coefficients from plane equation
                            if 'z' in plane_equation_text:
                                # Case with z term
                                pattern = r'([+-]?\d*\.?\d+)x\s*([+-]\s*\d*\.?\d+)y\s*([+-]\s*\d*\.?\d+)z\s*([+-]\s*\d*\.?\d+)\s*=\s*0'
                                match = re.search(pattern, plane_equation_text)
                                if match:
                                    norm_h = float(match.group(1))
                                    norm_k = float(match.group(2).replace(' ', ''))
                                    norm_l = float(match.group(3).replace(' ', ''))
                                    reference_normal = (norm_h, norm_k, norm_l)
                            else:
                                # Case without z term - use general parsing method
                                norm_h = 0.0  # x coefficient
                                norm_k = 0.0  # y coefficient 
                                norm_l = 0.0  # z coefficient (always 0)
                                
                                # Find x term
                                x_pattern = r'([+-]?\d*\.?\d*)x'
                                x_match = re.search(x_pattern, plane_equation_text)
                                if x_match:
                                    x_coeff = x_match.group(1)
                                    if x_coeff == '' or x_coeff == '+':
                                        norm_h = 1.0
                                    elif x_coeff == '-':
                                        norm_h = -1.0
                                    else:
                                        norm_h = float(x_coeff)
                                
                                # Find y term
                                y_pattern = r'([+-]?\d*\.?\d*)y'
                                y_match = re.search(y_pattern, plane_equation_text)
                                if y_match:
                                    y_coeff = y_match.group(1)
                                    if y_coeff == '' or y_coeff == '+':
                                        norm_k = 1.0
                                    elif y_coeff == '-':
                                        norm_k = -1.0
                                    else:
                                        norm_k = float(y_coeff)
                                
                                reference_normal = (norm_h, norm_k, norm_l)
                            break
            except Exception as e:
                print(f"   Warning: Failed to read normal vector from CSV file, using default: {str(e)}")
            
            # Output table header
            print("─" * 125)
            print("Moved Atoms and Coordinates                   | Reference Plane Atoms and Coordinates | Normal Vector Dot Product and In-Plane Vector Match | Position Match")
            print("─" * 125)
            
            for i, moved_atom in enumerate(moved_atoms):
                moved_id = moved_atom.get('atom_id', f'MovedAtom_{i+1}')
                moved_element = moved_atom.get('element', 'Unknown')
                moved_pos = (
                    float(moved_atom['cart_x']), 
                    float(moved_atom['cart_y']), 
                    float(moved_atom['cart_z'])
                )
                
                # Find closest atom in reference plane
                min_distance = float('inf')
                closest_ref_atom = None
                closest_distance = None
                
                for ref_atom in reference_plane_atoms:
                    ref_pos = (
                        float(ref_atom['cart_x']), 
                        float(ref_atom['cart_y']), 
                        float(ref_atom['cart_z'])
                    )
                    
                    # Calculate Euclidean distance
                    distance = math.sqrt(
                        (moved_pos[0] - ref_pos[0])**2 + 
                        (moved_pos[1] - ref_pos[1])**2 + 
                        (moved_pos[2] - ref_pos[2])**2
                    )
                    
                    if distance < min_distance:
                        min_distance = distance
                        closest_ref_atom = ref_atom
                        closest_distance = distance
                
                if closest_ref_atom:
                    # Use key that matches actual key structure
                    ref_id = (closest_ref_atom.get('atom_id') or 
                             closest_ref_atom.get('label', f"{closest_ref_atom.get('element', 'Unknown')}{len(reference_plane_atoms)}"))
                    # Use element name key that matches actual key structure
                    ref_element = closest_ref_atom.get('element', 'Unknown')
                    ref_pos = (
                        float(closest_ref_atom['cart_x']), 
                        float(closest_ref_atom['cart_y']), 
                        float(closest_ref_atom['cart_z'])
                    )
                    
                    # Calculate distance vector (replicated atom → reference plane atom)
                    distance_vector = (
                        ref_pos[0] - moved_pos[0],
                        ref_pos[1] - moved_pos[1],
                        ref_pos[2] - moved_pos[2]
                    )
                    
                    # Calculate dot product: if vector is in-plane, dot product should be close to 0
                    dot_product = (distance_vector[0] * reference_normal[0] + 
                                 distance_vector[1] * reference_normal[1] + 
                                 distance_vector[2] * reference_normal[2])
                    
                    # in-plane tolerance (consider in-plane vector if within 0.01Å)
                    in_plane_tolerance = 0.01
                    is_in_plane = abs(dot_product) <= in_plane_tolerance
                    
                    is_at_same_position = closest_distance <= tolerance
                    
                    # Output in table format
                    moved_info = f"{moved_id} ({moved_element}) @ ({moved_pos[0]:.4f}, {moved_pos[1]:.4f}, {moved_pos[2]:.4f})"
                    ref_info = f"{ref_id} ({ref_element}) @ ({ref_pos[0]:.4f}, {ref_pos[1]:.4f}, {ref_pos[2]:.4f})"
                    normal_info = f"{dot_product:.6f} (tolerance ±{in_plane_tolerance}):  {'Yes' if is_in_plane else 'No'}"
                    position_info = f"{'Pass' if is_at_same_position else 'Fail'}"
                    
                    print(f"{moved_info:<42} | {ref_info} | {normal_info} | {position_info}")
                    
                    # Output warning if not in-plane
                    if not is_in_plane:
                        print(f"   Warning: Distance vector is not perpendicular to reference plane! (atoms may be on different planes)")
                    
                    distance_results.append({
                        'moved_atom': {
                            'id': moved_id,
                            'element': moved_element,
                            'position': moved_pos
                        },
                        'closest_reference_atom': {
                            'id': ref_id,
                            'element': ref_element,
                            'position': ref_pos
                        },
                        'distance': closest_distance,
                        'distance_vector': distance_vector,
                        'dot_product_with_normal': dot_product,
                        'is_in_plane': is_in_plane,
                        'position_match': is_at_same_position
                    })
                    
                    if not is_at_same_position:
                        all_passed = False
                else:
                    print(f"   Cannot find close atom in reference plane.")
                    all_passed = False
            
            print("─" * 125)
            
            # Result summary
            result = {
                'test_step': 1,
                'test_name': 'Position Distance Check',
                'target_plane': target_plane,
                'tolerance': tolerance,
                'total_atoms': len(moved_atoms),
                'passed_atoms': sum(1 for r in distance_results if r['position_match']),
                'all_passed': all_passed,
                'distance_results': distance_results
            }
            
            print(f"Step 1 test result summary: Total atoms {result['total_atoms']} | Passed {result['passed_atoms']} | {'Equivalence test step 1 passed: All atoms are at the same position.' if all_passed else 'Equivalence test step 1 failed: Some atoms are at different positions.'}")
            
            return result
            
        except Exception as e:
            print(f"Failed to execute equivalence test step 1: {str(e)}")
            return None

    def perform_equivalence_test_step2(self, moved_group, reference_plane_atoms, step1_result=None, tolerance=0.01):
        """
        Equivalence test step 2: Element type verification
        Use atom pairs matched in step 1 to check if element types match.
        
        Args:
            moved_group (dict): Moved atom group information
            reference_plane_atoms (list): Atoms of the reference plane
            step1_result (dict): Step 1 test results (including matched atom pair information)
            tolerance (float): Allowed tolerance (not used, kept for compatibility)
            
        Returns:
            dict: Test results or None
        """

        
        if not moved_group or not reference_plane_atoms:
            print("Insufficient input data.")
            return None
        
        try:
            moved_atoms = moved_group.get('moved_atoms', [])
            if not moved_atoms:
                print("No moved atoms.")
                return None
            
            print(f"Reference plane: {moved_group.get('target_plane', 'Unknown')}")
            print(f"Number of moved atoms: {len(moved_atoms)}")
            print(f"Number of reference plane atoms: {len(reference_plane_atoms)}")
            
            results = []
            all_passed = True
            
            # Use matched atom pairs from step 1 if available, otherwise use existing method
            if step1_result and 'distance_results' in step1_result:
                print("Using matched atom pairs from step 1.")
                
                for i, distance_result in enumerate(step1_result['distance_results'], 1):
                    moved_atom_info = distance_result['moved_atom']
                    ref_atom_info = distance_result['closest_reference_atom']
                    
                    # Find moved_atom object
                    moved_atom = None
                    for atom in moved_atoms:
                        if atom.get('atom_id') == moved_atom_info['id']:
                            moved_atom = atom
                            break
                    
                    if not moved_atom:
                        print(f"   Cannot find moved atom {moved_atom_info['id']}.")
                        all_passed = False
                        continue
                    
                    # Find reference plane atom
                    closest_ref_atom = None
                    for atom in reference_plane_atoms:
                        if atom.get('atom_id') == ref_atom_info['id']:
                            closest_ref_atom = atom
                            break
                    
                    if not closest_ref_atom:
                        print(f"   Cannot find reference plane atom {ref_atom_info['id']}.")
                        all_passed = False
                        continue
                    
                    moved_element = moved_atom.get('element', 'Unknown')
                    ref_element = closest_ref_atom.get('element', 'Unknown')
                    element_match = moved_element == ref_element
                    passed = element_match
                    
                    print(f"Moved atom {i}: {moved_atom_info['id']} ({moved_element}) @ ({moved_atom_info['position'][0]:.4f}, {moved_atom_info['position'][1]:.4f}, {moved_atom_info['position'][2]:.4f}) ↔ Reference plane {ref_atom_info['id']} ({ref_element}) @ ({ref_atom_info['position'][0]:.4f}, {ref_atom_info['position'][1]:.4f}, {ref_atom_info['position'][2]:.4f}) | Distance: {distance_result['distance']:.6f}Å | Element match: {'Pass' if element_match else 'Fail'}")
                    
                    results.append({
                        'moved_atom': moved_atom,
                        'closest_reference_atom': closest_ref_atom,
                        'distance': distance_result['distance'],
                        'element_match': element_match,
                        'passed': passed
                    })
                    
                    if not passed:
                        all_passed = False
            
            else:
                print("No step 1 results available, using existing matching method.")
                # Existing method: Find closest same element atom in reference plane for each moved atom
                for i, moved_atom in enumerate(moved_atoms, 1):
                    moved_element = moved_atom.get('element', 'Unknown')
                    moved_id = moved_atom.get('atom_id', f'Unknown_{i}')
                    moved_pos = (moved_atom['cart_x'], moved_atom['cart_y'], moved_atom['cart_z'])
                    
                    print(f"Moved atom {i}: {moved_id} ({moved_element}) @ ({moved_pos[0]:.4f}, {moved_pos[1]:.4f}, {moved_pos[2]:.4f})", end="")
                    
                    # Find closest atom among reference plane atoms with same element
                    same_element_refs = [atom for atom in reference_plane_atoms 
                                       if atom.get('element', '') == moved_element]
                    
                    if not same_element_refs:
                        print(f" | Reference plane has no {moved_element} element.")
                        results.append({
                            'moved_atom': moved_atom,
                            'closest_reference_atom': None,
                            'element_match': False,
                            'passed': False
                        })
                        all_passed = False
                        continue
                    
                    # Find closest same element atom
                    min_distance = float('inf')
                    closest_ref_atom = None
                    
                    for ref_atom in same_element_refs:
                        ref_pos = (ref_atom['cart_x'], ref_atom['cart_y'], ref_atom['cart_z'])
                        
                        # Calculate Euclidean distance
                        import math
                        distance = math.sqrt(
                            (moved_pos[0] - ref_pos[0])**2 + 
                            (moved_pos[1] - ref_pos[1])**2 + 
                            (moved_pos[2] - ref_pos[2])**2
                        )
                        
                        if distance < min_distance:
                            min_distance = distance
                            closest_ref_atom = ref_atom
                    
                    # Store result
                    if closest_ref_atom:
                        ref_id = closest_ref_atom.get('atom_id', 'Unknown')
                        ref_element = closest_ref_atom.get('element', 'Unknown')
                        ref_pos = (closest_ref_atom['cart_x'], closest_ref_atom['cart_y'], closest_ref_atom['cart_z'])
                        
                        element_match = moved_element == ref_element
                        passed = element_match
                        
                        print(f" ↔ Reference plane {ref_id} ({ref_element}) @ ({ref_pos[0]:.4f}, {ref_pos[1]:.4f}, {ref_pos[2]:.4f}) | Distance: {min_distance:.6f}Å | Element match: {'Pass' if element_match else 'Fail'}")
                        
                        results.append({
                            'moved_atom': moved_atom,
                            'closest_reference_atom': closest_ref_atom,
                            'distance': min_distance,
                            'element_match': element_match,
                            'passed': passed
                        })
                        
                        if not passed:
                            all_passed = False
                    else:
                        print(f" | Cannot find reference plane atom with same element.")
                        results.append({
                            'moved_atom': moved_atom,
                            'closest_reference_atom': None,
                            'element_match': False,
                            'passed': False
                        })
                        all_passed = False
            
            # Output result summary
            passed_count = sum(1 for r in results if r['passed'])
            print(f"\nStep 2 test result summary:")
            print(f"   - Total atoms: {len(results)}")
            print(f"   - Passed atoms: {passed_count}")
            print(f"   - Overall pass: {'Pass' if all_passed else 'Fail'}")
            
            if all_passed:
                print(f"Equivalence test step 2 passed: All atom elements match.")
            else:
                print(f"Equivalence test step 2 failed: Some atom elements do not match.")
            
            # Generate atom pair matching information for use in Step 3
            atom_pairs = []
            for result in results:
                if result['closest_reference_atom']:
                    atom_pairs.append({
                        'comparison_atom': result['moved_atom'],
                        'reference_atom': result['closest_reference_atom']
                    })
            
            return {
                'all_passed': all_passed,
                'passed_count': passed_count,
                'total_count': len(results),
                'detailed_results': results,
                'atom_pairs': atom_pairs  # Add matching information for use in Step 3
            }
            
        except Exception as e:
            print(f"Failed to execute equivalence test step 2: {str(e)}")
            return None

    def perform_equivalence_test_step3(self, comparison_atoms, reference_atoms, reference_plane_info, all_supercell_atoms, step2_result=None):
        """
        Equivalence test step 3: Bond analysis comparison
        Compare bond information between comparison plane atoms and reference plane atoms.
        
        Args:
            comparison_atoms (list): Atoms of the comparison plane
            reference_atoms (list): Matched reference plane atoms
            reference_plane_info (dict): Reference plane information
            all_supercell_atoms (list): Complete supercell atom information
            
        Returns:
            dict: Test results or None
        """

        
        # Input data validation
        if not comparison_atoms or not reference_atoms:
            print("Insufficient input atom data.")
            return None
            
        if not all_supercell_atoms:
            print("Insufficient complete supercell atom data.")
            return None
        
        if not hasattr(self, 'bond_criteria') or not self.bond_criteria:
            print("Bond analysis criteria not set.")
            return None
        
        try:
            # Output basic information
            plane_id = reference_plane_info.get('plane_id', 'Unknown') if reference_plane_info else 'Unknown'
            print(f"Reference plane: {plane_id}")
            print(f"Number of comparison plane atoms: {len(comparison_atoms)}")
            print(f"Number of reference plane atoms: {len(reference_atoms)}")
            
            bond_type = self.bond_criteria.get('bond_type', 'Unknown')
            min_dist = self.bond_criteria.get('min_distance', 0.0)
            max_dist = self.bond_criteria.get('max_distance', 10.0)
            print(f"Bond analysis criteria: {bond_type} ({min_dist:.1f}Å ~ {max_dist:.1f}Å)")
            
            # Use atom pair matching information established in Step 2
            if step2_result and step2_result.get('atom_pairs'):
                print("Using matched atom pairs from step 2.")
                atom_pairs = step2_result['atom_pairs']
            else:
                print("Warning: No step 2 matching information available, matching in order.")
                atom_pairs = []
                min_count = min(len(comparison_atoms), len(reference_atoms))
                for i in range(min_count):
                    atom_pairs.append({
                        'comparison_atom': comparison_atoms[i],
                        'reference_atom': reference_atoms[i]
                    })
            
            # 1. Bond analysis of comparison plane atoms (based on matched pairs)
            print(f"\nBond analysis of comparison plane atoms:")
            comparison_bond_info = []
            
            for i, pair in enumerate(atom_pairs):
                atom = pair['comparison_atom']
                atom_id = atom.get('atom_id', f'CompAtom_{i+1}')
                element = atom.get('element', 'Unknown')
                atom_pos = (atom['cart_x'], atom['cart_y'], atom['cart_z'])
                
                # For cloned atoms, find original atom and calculate bond vectors
                original_atom = None
                original_atom_id = None
                if atom_id.endswith('A'):  # If it's a cloned atom
                    original_atom_id = atom_id[:-1]  # Remove 'A' suffix
                    
                    # Find original atom in the complete atom list
                    for orig_atom in all_supercell_atoms:
                        if orig_atom.get('atom_id') == original_atom_id:
                            original_atom = orig_atom
                            break
                
                # Find bond partners - based on original atom if available, otherwise current atom
                target_atom_for_bonds = original_atom if original_atom else atom
                bonds = self._find_bond_partners(target_atom_for_bonds, all_supercell_atoms, min_dist, max_dist)
                bond_count = len(bonds)
                
                if original_atom:
                    print(f"Atom {i+1}: {atom_id} ({element}) @ ({original_atom['cart_x']:.4f}, {original_atom['cart_y']:.4f}, {original_atom['cart_z']:.4f}) → ({atom_pos[0]:.4f}, {atom_pos[1]:.4f}, {atom_pos[2]:.4f}) [Original: {original_atom_id}]")
                else:
                    print(f"Atom {i+1}: {atom_id} ({element}) @ ({atom_pos[0]:.4f}, {atom_pos[1]:.4f}, {atom_pos[2]:.4f})")
                print(f"   Found bonds: {bond_count}")
                
                # Detailed output of all bond information
                for j, bond in enumerate(bonds):
                    partner_atom_id = bond.get('atom_id', 'Unknown')
                    partner_element = bond.get('element', 'Unknown')
                    partner_pos = bond.get('position', (0, 0, 0))
                    distance = bond.get('distance', 0.0)
                    vector = bond.get('bond_vector', (0, 0, 0))
                    
                    print(f"      {j+1}. {partner_atom_id} ({partner_element}) - Position({partner_pos[0]:.4f}, {partner_pos[1]:.4f}, {partner_pos[2]:.4f}) - {distance:.3f}Å - Vector({vector[0]:.3f}, {vector[1]:.3f}, {vector[2]:.3f})")
                
                if bond_count == 0:
                    print(f"      Warning: No bonds found!")
                
                comparison_bond_info.append({
                    'atom': atom,
                    'atom_id': atom_id,
                    'element': element,
                    'bonds': bonds,
                    'bond_count': bond_count
                })
            
            # 2. Bond analysis of reference plane atoms (based on matched pairs)
            print(f"\nBond analysis of reference plane atoms:")
            reference_bond_info = []
            
            for i, pair in enumerate(atom_pairs):
                atom = pair['reference_atom']
                atom_id = atom.get('atom_id', f'RefAtom_{i+1}')
                element = atom.get('element', 'Unknown')
                atom_pos = (atom['cart_x'], atom['cart_y'], atom['cart_z'])
                
                # Find bond partners
                bonds = self._find_bond_partners(atom, all_supercell_atoms, min_dist, max_dist)
                bond_count = len(bonds)
                
                print(f"Atom {i+1}: {atom_id} ({element}) @ ({atom_pos[0]:.4f}, {atom_pos[1]:.4f}, {atom_pos[2]:.4f})")
                print(f"   Found bonds: {bond_count}")
                
                # Detailed output of all bond information
                for j, bond in enumerate(bonds):
                    partner_atom_id = bond.get('atom_id', 'Unknown')
                    partner_element = bond.get('element', 'Unknown')
                    partner_pos = bond.get('position', (0, 0, 0))
                    distance = bond.get('distance', 0.0)
                    vector = bond.get('bond_vector', (0, 0, 0))
                    
                    print(f"      {j+1}. {partner_atom_id} ({partner_element}) - Position({partner_pos[0]:.4f}, {partner_pos[1]:.4f}, {partner_pos[2]:.4f}) - {distance:.3f}Å - Vector({vector[0]:.3f}, {vector[1]:.3f}, {vector[2]:.3f})")
                
                if bond_count == 0:
                    print(f"      Warning: No bonds found!")
                
                reference_bond_info.append({
                    'atom': atom,
                    'atom_id': atom_id,
                    'element': element,
                    'bonds': bonds,
                    'bond_count': bond_count
                })
            
            # 3. Bond pattern comparison
            print(f"\nBond pattern comparison analysis:")
            comparison_results = []
            overall_passed = True
            
            # Check atom count
            if len(comparison_bond_info) != len(reference_bond_info):
                print(f"Warning: Different atom counts: Comparison plane({len(comparison_bond_info)}) vs Reference plane({len(reference_bond_info)})")
                overall_passed = False
            
            # Compare each atom pair (based on matched pairs)
            for i in range(len(atom_pairs)):
                comp_info = comparison_bond_info[i]
                ref_info = reference_bond_info[i]
                
                comp_atom_id = comp_info['atom_id']
                ref_atom_id = ref_info['atom_id']
                comp_element = comp_info['element']
                ref_element = ref_info['element']
                comp_bond_count = comp_info['bond_count']
                ref_bond_count = ref_info['bond_count']
                
                print(f"Atom {i+1} comparison:")
                print(f"   Comparison plane: {comp_atom_id} ({comp_element}) - {comp_bond_count} bonds | Reference plane: {ref_atom_id} ({ref_element}) - {ref_bond_count} bonds")
                
                # Element type comparison
                element_match = (comp_element == ref_element)
                
                # Bond count comparison
                bond_count_match = (comp_bond_count == ref_bond_count)
                
                # Bond distance pattern comparison
                distance_pattern_match = True
                if bond_count_match and comp_bond_count > 0:
                    # Sort distance lists for comparison
                    comp_distances = sorted([bond['distance'] for bond in comp_info['bonds']])
                    ref_distances = sorted([bond['distance'] for bond in ref_info['bonds']])
                    
                    for comp_dist, ref_dist in zip(comp_distances, ref_distances):
                        if abs(comp_dist - ref_dist) > 0.01:  # 0.01Å tolerance
                            distance_pattern_match = False
                            break
                
                print(f"   Element match: {'Yes' if element_match else 'No'} | Bond count match: {'Yes' if bond_count_match else 'No'} | Distance pattern match: {'Yes' if distance_pattern_match else 'No'}")
                
                # Bond vector pattern comparison (distance + direction)
                vector_pattern_match = True
                if bond_count_match and comp_bond_count > 0:
                    # Extract bond vectors from comparison plane and reference plane
                    comp_vectors = []
                    ref_vectors = []
                    
                    for bond in comp_info['bonds']:
                        bond_vector = bond.get('bond_vector', [0, 0, 0])
                        # Convert numpy type to regular float and limit to 3 decimal places
                        clean_vector = [round(float(x), 3) for x in bond_vector]
                        comp_vectors.append(tuple(clean_vector))
                    
                    for bond in ref_info['bonds']:
                        bond_vector = bond.get('bond_vector', [0, 0, 0])
                        # Convert numpy type to regular float and limit to 3 decimal places
                        clean_vector = [round(float(x), 3) for x in bond_vector]
                        ref_vectors.append(tuple(clean_vector))
                    
                    # Sort vectors for comparison (considering both direction and distance)
                    comp_vectors_sorted = sorted(comp_vectors)
                    ref_vectors_sorted = sorted(ref_vectors)
                    
                    # Compare each vector pair
                    for idx, (comp_vec, ref_vec) in enumerate(zip(comp_vectors_sorted, ref_vectors_sorted)):
                        # Compare by vector components (tolerance relaxed from 0.005 to 0.02)
                        max_diff = 0.0
                        component_diffs = []
                        for i in range(3):
                            diff = abs(comp_vec[i] - ref_vec[i])
                            # Ignore very small differences (floating point errors)
                            if diff < 1e-10:
                                diff = 0.0
                            component_diffs.append(diff)
                            max_diff = max(max_diff, diff)
                            if diff > 0.02:  # Tolerance relaxed from 0.01 to 0.02
                                vector_pattern_match = False
                        
                        # Display result in one line
                        status = "Yes" if max_diff <= 0.02 else "No"
                        print(f"      Vector pair {idx+1}: {comp_vec} vs {ref_vec} → Difference: ({component_diffs[0]:.6f}, {component_diffs[1]:.6f}, {component_diffs[2]:.6f}) {status}")
                        
                        if not vector_pattern_match:
                            print(f"         Mismatch found in vector pair {idx+1} (max difference: {max_diff:.6f} > 0.02), stopping")
                            break
                    
                    print(f"   Bond vector match: {'Yes' if vector_pattern_match else 'No'}")
                    
                    # Output detailed comparison information (only when there's a mismatch)
                    if not vector_pattern_match:
                        print(f"   Vector mismatch details:")
                        print(f"      Comparison plane vectors: {comp_vectors_sorted}")
                        print(f"      Reference plane vectors: {ref_vectors_sorted}")
                
                # Individual atom pass/fail (including vector pattern match)
                atom_passed = element_match and bond_count_match and vector_pattern_match
                
                comparison_results.append({
                    'comparison_atom': comp_info['atom'],
                    'reference_atom': ref_info['atom'],
                    'element_match': element_match,
                    'bond_count_match': bond_count_match,
                    'distance_pattern_match': distance_pattern_match,
                    'vector_pattern_match': vector_pattern_match,
                    'passed': atom_passed
                })
                
                if not atom_passed:
                    overall_passed = False
            
            # Result summary
            passed_count = sum(1 for result in comparison_results if result['passed'])
            total_count = len(comparison_results)
            
            print(f"\nStep 3 test result summary:")
            print(f"   - Total atoms: {total_count}")
            print(f"   - Passed atoms: {passed_count}")
            print(f"   - Overall pass: {'Pass' if overall_passed else 'Fail'}")
            
            if overall_passed:
                print(f"Equivalence test step 3 passed: All atom bond patterns match.")
            else:
                print(f"Equivalence test step 3 failed: Some atom bond patterns do not match.")
            
            return {
                'all_passed': overall_passed,
                'passed_count': passed_count,
                'total_count': total_count,
                'detailed_results': comparison_results,
                'comparison_bond_results': comparison_bond_info,
                'reference_bond_results': reference_bond_info
            }
            
        except Exception as e:
            print(f"Failed to execute equivalence test step 3: {str(e)}")
            import traceback
            traceback.print_exc()
            return None 

    def update_csv_for_equivalent_plane(self, csv_file_path, plane_id, reference_plane_id, plane_type):
        """
        Update CSV file for planes that passed all equivalence tests.
        
        Args:
            csv_file_path (str): CSV file path
            plane_id (str): Plane ID that passed equivalence test
            reference_plane_id (str): Reference plane ID
            plane_type (str): Plane type (e.g., "Type A")
            
        Returns:
            bool: Update success status
        """
        print(f"\nUpdating CSV for equivalent plane: {plane_id}")
        
        try:
            if not os.path.exists(csv_file_path):
                print(f"CSV file does not exist: {csv_file_path}")
                return False
            
            # Read CSV file
            rows = []
            updated_count = 0
            
            with open(csv_file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)
                rows.append(header)
                
                # Find column indices in header
                plane_id_idx = header.index('plane_id')
                ref_plane_idx = header.index('reference_plane') 
                plane_type_idx = header.index('plane_type')
                
                for row in reader:
                    try:
                        # Update atoms belonging to the plane
                        if row[plane_id_idx] == plane_id:
                            # reference_plane_id is passed in "Reference Plane X" format, so use as is
                            row[ref_plane_idx] = reference_plane_id
                            row[plane_type_idx] = plane_type
                            updated_count += 1
                        
                        rows.append(row)
                        
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Error processing row: {e}")
                        rows.append(row)
            
            # Save updated CSV file
            with open(csv_file_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerows(rows)
            
            print(f"CSV update completed: Plane={plane_id} | Reference plane={reference_plane_id} | Type={plane_type} | Updated atoms={updated_count}")
            
            return True
            
        except Exception as e:
            print(f"CSV update failed: {str(e)}")
            return False
    
    def _get_reference_plane_name(self, csv_file_path, reference_plane_id):
        """
        Dynamically find reference plane mapping in CSV file and generate correct reference name.
        
        Args:
            csv_file_path (str): CSV file path
            reference_plane_id (str): Reference plane ID (e.g., "Plane 18")
            
        Returns:
            str: Correct reference name (e.g., "Equivalent to Reference Plane 1")
        """
        try:
            # Find current reference plane mapping in CSV file
            with open(csv_file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)
                
                plane_id_idx = header.index('plane_id')
                ref_plane_idx = header.index('reference_plane')
                
                # Find which reference plane the reference_plane_id is set to  
                for row in reader:
                    if row[plane_id_idx] == reference_plane_id:
                        reference_plane_name = row[ref_plane_idx]  # e.g., "Reference Plane 1"
                        if reference_plane_name.endswith("Reference Plane"):
                            return f"Equivalent to {reference_plane_name}"
                        break
            
            # Use default value if not found
            return f"Equivalent to {reference_plane_id}"
            
        except Exception as e:
            print(f"⚠️  Failed to find reference plane mapping: {e}")
            return f"Equivalent to {reference_plane_id}"
    
 