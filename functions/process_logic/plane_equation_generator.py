import numpy as np
from .plane_calculator import PlaneCalculator

class PlaneEquationGenerator:
    """Class for generating plane equations from Miller indices"""
    
    def __init__(self):
        """Initialize PlaneEquationGenerator"""
        self.plane_calculator = PlaneCalculator()
    
    def generate_plane_equation(self, miller_indices, cell_vectors):
        """
        Generate plane equation from Miller indices and cell vectors.
        
        Args:
            miller_indices (dict): Miller index information {'h': int, 'k': int, 'l': int}
            cell_vectors (dict): Cell vector information {'a': [x,y,z], 'b': [x,y,z], 'c': [x,y,z]}
            
        Returns:
            dict: Plane equation information or None (if failed)
        """
        print("\nGenerating plane equation...")
        
        try:
            # Extract Miller indices
            h = miller_indices['h']
            k = miller_indices['k']
            l = miller_indices['l']
            
            print(f"Miller indices: ({h} {k} {l})")
            print("Cell vector information:")
            for key, vector in cell_vectors.items():
                print(f"   {key} = [{vector[0]:8.4f}, {vector[1]:8.4f}, {vector[2]:8.4f}] Å")
            
            # Calculate plane equation using existing PlaneCalculator
            a, b, c, d = PlaneCalculator.calculate_plane_equation(h, k, l, cell_vectors)
            
            # Format plane equation
            equation_str = PlaneCalculator.format_plane_equation((a, b, c, d))
            
            print(f"SUCCESS: Calculated plane equation: {equation_str}")
            
            # Construct result data
            plane_equation = {
                'miller_indices': miller_indices,
                'coefficients': {'a': a, 'b': b, 'c': c, 'd': d},
                'equation_string': equation_str,
                'original_d': d  # Store original D value
            }
            
            return plane_equation
            
        except Exception as e:
            print(f"ERROR: Plane equation generation failed: {str(e)}")
            return None
    
    def create_parallel_planes_with_d_values(self, base_plane_equation, d_min, d_max, d_step):
        """
        Generate parallel planes by changing D values from base plane equation.
        
        Args:
            base_plane_equation (dict): Base plane equation information
            d_min (float): Minimum D value
            d_max (float): Maximum D value
            d_step (float): D value increment
            
        Returns:
            list: List of parallel planes
        """
        try:
            planes = []
            coeffs = base_plane_equation['coefficients']
            a, b, c = coeffs['a'], coeffs['b'], coeffs['c']
            
            # Generate parallel planes by varying D values
            d_current = d_min
            plane_id = 1
            
            while d_current <= d_max:
                # Generate new plane equation
                new_equation_str = PlaneCalculator.format_plane_equation((a, b, c, d_current))
                
                plane_info = {
                    'plane_id': plane_id,
                    'coefficients': {'a': a, 'b': b, 'c': c, 'd': d_current},
                    'equation_string': new_equation_str,
                    'd_value': d_current
                }
                
                planes.append(plane_info)
                d_current += d_step
                plane_id += 1
            
            return planes
            
        except Exception as e:
            print(f"ERROR: Parallel plane generation failed: {str(e)}")
            return []
    
    def format_plane_info(self, plane_equation):
        """
        Format and display plane equation information.
        
        Args:
            plane_equation (dict): Plane equation information
        """
        if not plane_equation:
            print("ERROR: Invalid plane equation information")
            return
        
        miller = plane_equation['miller_indices']
        print(f"\nMiller indices: ({miller['h']} {miller['k']} {miller['l']})")
        print(f"Plane equation: {plane_equation['equation_string']}")
        
        coeffs = plane_equation['coefficients']
        print(f"Coefficient information:")
        print(f"   A = {coeffs['a']:8.4f}")
        print(f"   B = {coeffs['b']:8.4f}")  
        print(f"   C = {coeffs['c']:8.4f}")
        print(f"   D = {coeffs['d']:8.4f}") 