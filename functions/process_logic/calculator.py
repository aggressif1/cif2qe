import math
import numpy as np

class CellCalculator:
    """Class for performing unit cell related calculations"""
    
    @staticmethod
    def calculate_volume(cell_params):
        """
        Calculate the volume of a unit cell.
        
        Args:
            cell_params (dict): Unit cell parameters
            
        Returns:
            float: Volume of unit cell (Å³)
        """
        a, b, c = cell_params['a'], cell_params['b'], cell_params['c']
        alpha, beta, gamma = cell_params['alpha'], cell_params['beta'], cell_params['gamma']
        
        # Check if all parameters are available
        if None in [a, b, c, alpha, beta, gamma]:
            return None
        
        # Convert angles to radians
        alpha_rad = math.radians(alpha)
        beta_rad = math.radians(beta)
        gamma_rad = math.radians(gamma)
        
        # Calculate volume
        volume = a * b * c * math.sqrt(
            1 - math.cos(alpha_rad)**2 - math.cos(beta_rad)**2 - math.cos(gamma_rad)**2 +
            2 * math.cos(alpha_rad) * math.cos(beta_rad) * math.cos(gamma_rad)
        )
        
        return volume
    
    @staticmethod
    def calculate_cell_vectors(cell_params):
        """
        Calculate unit cell vectors.
        
        Args:
            cell_params (dict): Unit cell parameters
            
        Returns:
            dict: Cell vectors {'a': [x, y, z], 'b': [x, y, z], 'c': [x, y, z]}
        """
        a, b, c = cell_params['a'], cell_params['b'], cell_params['c']
        alpha, beta, gamma = cell_params['alpha'], cell_params['beta'], cell_params['gamma']
        
        # Check if all parameters are available
        if None in [a, b, c, alpha, beta, gamma]:
            return None
        
        # Convert angles to radians
        alpha_rad = math.radians(alpha)
        beta_rad = math.radians(beta)
        gamma_rad = math.radians(gamma)
        
        # a vector is aligned with x-axis
        ax = a
        ay = 0.0
        az = 0.0
        
        # b vector is in xy plane
        bx = b * math.cos(gamma_rad)
        by = b * math.sin(gamma_rad)
        bz = 0.0
        
        # Calculate c vector
        cx = c * math.cos(beta_rad)
        cy = c * (math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) / math.sin(gamma_rad)
        cz = math.sqrt(c*c - cx*cx - cy*cy)
        
        return {
            'a': [ax, ay, az],
            'b': [bx, by, bz],
            'c': [cx, cy, cz]
        }
    
    @staticmethod
    def determine_crystal_system(cell_params):
        """
        Determine crystal system.
        
        Args:
            cell_params (dict): Unit cell parameters
            
        Returns:
            str: Crystal system name
        """
        a, b, c = cell_params['a'], cell_params['b'], cell_params['c']
        alpha, beta, gamma = cell_params['alpha'], cell_params['beta'], cell_params['gamma']
        
        def is_equal(x, y, tolerance=1e-3):
            """Compare two values for equality (with tolerance)"""
            return abs(x - y) < tolerance
        
        def is_angle_equal(x, y, tolerance=0.1):
            """Compare two angles for equality (with tolerance)"""
            return abs(x - y) < tolerance
        
        # Check if all angles are 90 degrees
        all_90 = (is_angle_equal(alpha, 90.0) and 
                 is_angle_equal(beta, 90.0) and 
                 is_angle_equal(gamma, 90.0))
        
        # Check if all angles are equal
        all_angles_equal = (is_angle_equal(alpha, beta) and 
                          is_angle_equal(beta, gamma))
        
        # Check if all lengths are equal
        all_lengths_equal = (is_equal(a, b) and is_equal(b, c))
        
        # Determine crystal system
        if all_lengths_equal and all_angles_equal:
            if is_angle_equal(alpha, 90.0):
                return "Cubic"
            else:
                return "Rhombohedral"
        
        elif all_90:
            if all_lengths_equal:
                return "Cubic"
            elif is_equal(a, b) or is_equal(b, c) or is_equal(a, c):
                return "Tetragonal"
            else:
                return "Orthorhombic"
        
        elif is_equal(a, b) and is_angle_equal(alpha, beta):
            return "Hexagonal" if is_angle_equal(gamma, 120.0) else "Tetragonal"
        
        elif is_equal(a, b):
            return "Tetragonal"
        
        else:
            return "Monoclinic" if is_angle_equal(alpha, 90.0) and is_angle_equal(gamma, 90.0) else "Triclinic" 