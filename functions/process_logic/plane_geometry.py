"""
Plane geometry functions separated from Phase3
"""
import numpy as np
import re

class PlaneGeometry:
    """Plane geometry calculation class"""
    
    def __init__(self):
        """Initialize PlaneGeometry"""
        pass
    
    def get_point_on_plane(self, plane_equation):
        """
        Find a point on the plane.
        
        Args:
            plane_equation (str): Plane equation string
            
        Returns:
            numpy.ndarray: Point on plane or None
        """
        try:
            print(f"   Plane equation: {plane_equation}")
            
            # Parse as ax + by + cz = d form
            a, b, c, d = 0, 0, 0, 0
            
            # Separate left and right sides by equals sign
            if '=' in plane_equation:
                left_side, right_side = plane_equation.split('=')
                left_side = left_side.strip()
                right_side = right_side.strip()
                
                # Constant on right side
                d = float(right_side) if right_side != '0' else 0
                
                # Extract x, y, z coefficients from left side
                pattern = r'([+-]?\s*\d*\.?\d*)\s*\*?\s*([xyz])'
                matches = re.findall(pattern, left_side)
                
                for coeff_str, var in matches:
                    coeff_str = coeff_str.replace(' ', '')  # Remove spaces
                    
                    if coeff_str == '' or coeff_str == '+':
                        coeff = 1.0
                    elif coeff_str == '-':
                        coeff = -1.0
                    else:
                        coeff = float(coeff_str)
                    
                    if var == 'x':
                        a = coeff
                    elif var == 'y':
                        b = coeff
                    elif var == 'z':
                        c = coeff
                
                # Extract constant term from left side (numbers without variables)
                const_pattern = r'([+-]\s*\d+\.?\d*)\s*(?![xyz])'
                const_matches = re.findall(const_pattern, left_side)
                
                if const_matches:
                    for const_str in const_matches:
                        const_str = const_str.replace(' ', '')
                        d = d - float(const_str)  # Move left side constant to right side
                
                print(f"   Coefficients: a={a:.4f}, b={b:.4f}, c={c:.4f}, d={d:.4f}")
                
                # Find a point on the plane
                if abs(c) > 1e-6:  # If z coefficient is not 0
                    point = np.array([0.0, 0.0, d/c])
                elif abs(b) > 1e-6:  # If y coefficient is not 0
                    point = np.array([0.0, d/b, 0.0])
                elif abs(a) > 1e-6:  # If x coefficient is not 0
                    point = np.array([d/a, 0.0, 0.0])
                else:
                    print("ERROR: Invalid plane equation")
                    return None
                
                print(f"   Calculated point: ({point[0]:.4f}, {point[1]:.4f}, {point[2]:.4f})")
                return point
            
            return None
                
        except Exception as e:
            print(f"ERROR: Failed to calculate point on plane: {str(e)}")
            return None
    
    def calculate_line_plane_intersection(self, start_point, direction, plane_equation):
        """
        Calculate intersection of line and plane.
        
        Args:
            start_point: Starting point of line
            direction: Direction vector of line (unit vector)
            plane_equation: Plane equation
            
        Returns:
            numpy.ndarray: Intersection point or None
        """
        try:
            print(f"   Line-plane intersection calculation:")
            print(f"     Start point: ({start_point[0]:.4f}, {start_point[1]:.4f}, {start_point[2]:.4f})")
            print(f"     Direction: ({direction[0]:.4f}, {direction[1]:.4f}, {direction[2]:.4f})")
            print(f"     Plane: {plane_equation}")
            
            # Parse as ax + by + cz = d form
            a, b, c, d = 0, 0, 0, 0
            
            # Separate left and right sides by equals sign
            if '=' in plane_equation:
                left_side, right_side = plane_equation.split('=')
                left_side = left_side.strip()
                right_side = right_side.strip()
                
                # Constant on right side
                d = float(right_side) if right_side != '0' else 0
                
                # Extract x, y, z coefficients from left side
                pattern = r'([+-]?\s*\d*\.?\d*)\s*\*?\s*([xyz])'
                matches = re.findall(pattern, left_side)
                
                for coeff_str, var in matches:
                    coeff_str = coeff_str.replace(' ', '')  # Remove spaces
                    
                    if coeff_str == '' or coeff_str == '+':
                        coeff = 1.0
                    elif coeff_str == '-':
                        coeff = -1.0
                    else:
                        coeff = float(coeff_str)
                    
                    if var == 'x':
                        a = coeff
                    elif var == 'y':
                        b = coeff
                    elif var == 'z':
                        c = coeff
                
                # Extract constant term from left side (numbers without variables)
                const_pattern = r'([+-]\s*\d+\.?\d*)\s*(?![xyz])'
                const_matches = re.findall(const_pattern, left_side)
                
                if const_matches:
                    for const_str in const_matches:
                        const_str = const_str.replace(' ', '')
                        d = d - float(const_str)  # Move left side constant to right side
            
            normal = np.array([a, b, c])
            print(f"     Plane normal vector: ({a:.4f}, {b:.4f}, {c:.4f})")
            print(f"     Plane constant: {d:.4f}")
            
            # Line: P = start_point + t * direction
            # Plane: normal · P = d
            # Intersection: normal · (start_point + t * direction) = d
            # Therefore: t = (d - normal · start_point) / (normal · direction)
            
            denominator = np.dot(normal, direction)
            
            if abs(denominator) < 1e-10:
                print("ERROR: Line and plane are parallel")
                return None
            
            numerator = d - np.dot(normal, start_point)
            t = numerator / denominator
            
            # Calculate intersection point
            intersection = start_point + t * direction
            
            print(f"     Parameter t: {t:.4f}")
            print(f"     Intersection point: ({intersection[0]:.4f}, {intersection[1]:.4f}, {intersection[2]:.4f})")
            
            return intersection
            
        except Exception as e:
            print(f"ERROR: Line-plane intersection calculation failed: {str(e)}")
            return None
    
    def calculate_distance_to_next_same_plane(self, sequence_info, c_direction):
        """
        Calculate distance from first reference plane to next first reference plane 
        when moving in c_user direction.
        
        Args:
            sequence_info: Sequence analysis result
            c_direction: c_user direction unit vector
            
        Returns:
            float: Distance or None
        """
        try:
            print(f"\nCalculating distance: First reference plane → Next first reference plane:")
            
            # 1. Find next first reference plane position from sequence
            pattern = sequence_info['sequence_pattern']['pattern']
            pattern_length = sequence_info['sequence_pattern']['pattern_length']
            
            # First reference plane is usually first in pattern (A or B)
            # For BA pattern: B(0) → A(1) → B(2), so next B is at index 2
            next_same_plane_index = pattern_length
            print(f"   Pattern: {pattern}")
            print(f"   First reference plane index: 0")
            print(f"   Next first reference plane index: {next_same_plane_index}")
            
            # 2. Get equations for first reference plane and next first reference plane
            first_plane_eq = self.get_plane_equation_by_sequence_index(0)
            next_plane_eq = self.get_plane_equation_by_sequence_index(next_same_plane_index)
            
            if not first_plane_eq or not next_plane_eq:
                print("ERROR: Cannot retrieve plane equations")
                return None
            
            print(f"   First reference plane equation: {first_plane_eq}")
            print(f"   Next first reference plane equation: {next_plane_eq}")
            
            # 3. Find starting point on first reference plane
            start_point = self.get_point_on_plane(first_plane_eq)
            if start_point is None:
                print("ERROR: Cannot find starting point")
                return None
            
            # 4. Draw line in c_user direction to find intersection with next plane
            intersection = self.calculate_line_plane_intersection(start_point, c_direction, next_plane_eq)
            if intersection is None:
                print("ERROR: Cannot find intersection point")
                return None
            
            # 5. Calculate distance
            distance = np.linalg.norm(intersection - start_point)
            print(f"   Calculated distance: {distance:.4f} Å")
            
            return distance
            
        except Exception as e:
            print(f"ERROR: Distance calculation failed: {str(e)}")
            return None
    
    def get_plane_equation_by_sequence_index(self, sequence_index):
        """
        Get plane equation corresponding to sequence index.
        
        Args:
            sequence_index (int): Index in sequence
            
        Returns:
            str: Plane equation or None
        """
        try:
            import csv
            import os
            
            # Read all plane information from CSV file
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                return None
            
            # Collect plane information (only those with plane_type)
            planes_info = {}
            
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    plane_id = row.get('plane_id', '')
                    plane_equation = row.get('plane_equation', '')
                    plane_type = row.get('plane_type', '')
                    
                    if plane_id and plane_equation and plane_type and plane_type != 'Unknown':
                        # Extract D value
                        from .plane_analyzer import PlaneAnalyzer
                        analyzer = PlaneAnalyzer()
                        d_value = analyzer.extract_d_value_from_equation(plane_equation)
                        
                        if d_value is not None:
                            planes_info[plane_id] = {
                                'plane_type': plane_type,
                                'd_value': d_value,
                                'plane_equation': plane_equation
                            }
            
            # Sort by D value
            sorted_planes = sorted(planes_info.items(), key=lambda x: x[1]['d_value'])
            
            # Return plane corresponding to sequence index
            if 0 <= sequence_index < len(sorted_planes):
                plane_id, plane_info = sorted_planes[sequence_index]
                return plane_info['plane_equation']
            
            return None
            
        except Exception as e:
            print(f"ERROR: Failed to retrieve plane equation: {str(e)}")
            return None 