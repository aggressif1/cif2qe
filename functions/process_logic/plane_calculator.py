import numpy as np

class PlaneCalculator:
    """Class for calculating plane equations from Miller indices"""
    
    @staticmethod
    def calculate_reciprocal_vectors(cell_vectors):
        """
        Calculate reciprocal lattice vectors from real space lattice vectors.
        
        Args:
            cell_vectors (dict): Real space lattice vectors
            
        Returns:
            tuple: (a*, b*, c*) reciprocal lattice vectors
        """
        # Real space vectors
        a = np.array(cell_vectors['a'])
        b = np.array(cell_vectors['b'])
        c = np.array(cell_vectors['c'])
        
        # Unit cell volume
        volume = np.abs(np.dot(a, np.cross(b, c)))
        
        # Calculate reciprocal lattice vectors
        a_star = 2 * np.pi * np.cross(b, c) / volume
        b_star = 2 * np.pi * np.cross(c, a) / volume
        c_star = 2 * np.pi * np.cross(a, b) / volume
        
        return a_star, b_star, c_star
    
    @staticmethod
    def calculate_plane_equation(h, k, l, cell_vectors):
        """
        Calculate plane equation coefficients from Miller indices and cell vectors.
        Plane equation: ax + by + cz + d = 0
        
        Args:
            h, k, l (int): Miller indices
            cell_vectors (dict): Cell vector information {'a': [x,y,z], 'b': [x,y,z], 'c': [x,y,z]}
            
        Returns:
            tuple: (a, b, c, d) plane equation coefficients (float)
        """
        # Calculate reciprocal lattice vectors
        a_star, b_star, c_star = PlaneCalculator.calculate_reciprocal_vectors(cell_vectors)
        
        # Calculate plane normal vector (Miller indices in reciprocal space)
        # G = ha* + kb* + lc*
        normal = h * a_star + k * b_star + l * c_star
        
        # Normalize normal vector
        norm = np.linalg.norm(normal)
        if norm > 0:
            normal = normal / norm
            
        # Plane equation coefficients (components of normalized normal vector)
        a, b, c = normal.tolist()  # Convert numpy array to Python list
        
        # Calculate d coefficient (so plane passes through origin)
        d = 0.0
        
        return float(a), float(b), float(c), float(d)
    
    @staticmethod
    def format_plane_equation(plane_eq):
        """
        Format plane equation coefficients as Ax + By + Cz + D = 0 string.
        
        Args:
            plane_eq (tuple): (a, b, c, d) plane equation coefficients
            
        Returns:
            str: Formatted plane equation string
        """
        if isinstance(plane_eq, (tuple, list)) and len(plane_eq) >= 4:
            # Extract coefficients from tuple or list
            a, b, c, d = plane_eq[0], plane_eq[1], plane_eq[2], plane_eq[3]
        else:
            return str(plane_eq)  # Return as is if unexpected format
            
        # Convert numpy.float64 type to float
        try:
            a = float(a)
            b = float(b) 
            c = float(c)
            d = float(d)
        except (ValueError, TypeError):
            return str(plane_eq)  # Return original if conversion fails
        
        # Collect only valid terms
        terms = []
        
        # x term
        if abs(a) >= 1e-10:  # Only if not zero
            if a == 1.0:
                terms.append("x")
            elif a == -1.0:
                terms.append("-x")
            else:
                terms.append(f"{a:.4f}x")
        
        # y term
        if abs(b) >= 1e-10:  # Only if not zero
            if b == 1.0:
                if terms:  # If there are previous terms
                    terms.append(" + y")
                else:
                    terms.append("y")
            elif b == -1.0:
                terms.append(" - y")
            elif b > 0:
                if terms:  # If there are previous terms
                    terms.append(f" + {b:.4f}y")
                else:
                    terms.append(f"{b:.4f}y")
            else:  # b < 0
                terms.append(f" - {abs(b):.4f}y")
        
        # z term
        if abs(c) >= 1e-10:  # Only if not zero
            if c == 1.0:
                if terms:  # If there are previous terms
                    terms.append(" + z")
                else:
                    terms.append("z")
            elif c == -1.0:
                terms.append(" - z")
            elif c > 0:
                if terms:  # If there are previous terms
                    terms.append(f" + {c:.4f}z")
                else:
                    terms.append(f"{c:.4f}z")
            else:  # c < 0
                terms.append(f" - {abs(c):.4f}z")
        
        # Constant term
        if abs(d) >= 1e-10:  # Only if not zero
            if d > 0:
                if terms:  # If there are previous terms
                    terms.append(f" + {d:.4f}")
                else:
                    terms.append(f"{d:.4f}")
            else:  # d < 0
                terms.append(f" - {abs(d):.4f}")
        
        # If no terms, return "0"
        if not terms:
            return "0 = 0"
        
        # Connect all terms to complete plane equation
        equation = "".join(terms) + " = 0"
        
        return equation
    
    def find_atoms_on_plane(self, atoms, plane_eq, tolerance=0.0005):
        """
        Find atoms belonging to a given plane equation.
        
        Args:
            atoms (list): List of atomic information. Each atom has format {'element': str, 'cart_x': float, 'cart_y': float, 'cart_z': float}
            plane_eq (tuple): (a, b, c, d) plane equation coefficients
            tolerance (float): Tolerance (default: 0.01Å)
            
        Returns:
            list: List of atoms belonging to the plane
        """
        a, b, c, d = plane_eq
        atoms_on_plane = []
        
        for atom in atoms:
            x = atom['cart_x']
            y = atom['cart_y']
            z = atom['cart_z']
            
            # Calculate distance by substituting coordinates into plane equation
            distance = abs(a*x + b*y + c*z + d)
            
            # Add atom to list if distance is within tolerance
            if distance <= tolerance:
                atoms_on_plane.append(atom)
        
        return atoms_on_plane 