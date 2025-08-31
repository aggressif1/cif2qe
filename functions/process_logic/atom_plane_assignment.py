import numpy as np
from .plane_calculator import PlaneCalculator

class AtomPlaneAssignment:
    """Class for assigning atoms to planes"""
    
    def __init__(self):
        """Initialize AtomPlaneAssignment"""
        self.plane_calculator = PlaneCalculator()
    
    def assign_atoms_to_planes(self, supercell_atoms, plane_equation, tolerance=0.0005):
        """
        Assign supercell atoms to planes.
        Uses optimized method that generates planes based on atoms.
        
        Args:
            supercell_atoms (list): Supercell atom information
            plane_equation (dict): Basic plane equation information
            tolerance (float): Tolerance for plane assignment (default: 0.0005Å)
            
        Returns:
            dict: Plane assignment result or None (if failed)
        """
        print("\nStarting atom-based plane assignment...")
        
        try:
            # Perform atom-based plane assignment
            plane_assignments = self._assign_atoms_by_seed_method(
                supercell_atoms, plane_equation, tolerance
            )
            
            if not plane_assignments:
                print("ERROR: Plane assignment failed")
                return None
            
            # Calculate statistics
            total_planes = len(plane_assignments)
            total_atoms = sum(assignment['atom_count'] for assignment in plane_assignments.values())
            
            print(f"SUCCESS: Plane assignment completed:")
            print(f"   - Generated planes: {total_planes}")
            print(f"   - Assigned atoms: {total_atoms}")
            
            # Calculate projection range (for compatibility)
            projection_range = self._calculate_projection_range(supercell_atoms, plane_equation)
            
            # Construct result data (maintaining existing interface compatibility)
            result = {
                'plane_assignments': plane_assignments,  # Structure expected by existing code
                'projection_range': projection_range,    # Existing code compatibility
                'plane_spacing': 'adaptive',             # Using adaptive spacing
                'tolerance': tolerance,
                'total_planes': total_planes,
                'total_atoms': total_atoms,
                'method': 'seed_based_optimization'
            }
            
            return result
            
        except Exception as e:
            print(f"ERROR: Atom plane assignment failed: {str(e)}")
            return None
    
    def _assign_atoms_by_seed_method(self, supercell_atoms, plane_equation, tolerance):
        """
        Optimized method to generate and assign planes using atoms as seeds
        
        Args:
            supercell_atoms (list): Supercell atom information
            plane_equation (dict): Basic plane equation information
            tolerance (float): Tolerance
            
        Returns:
            dict: Atom assignment result per plane
        """
        print("Starting seed-based plane generation...")
        
        # Initialize
        unassigned_atoms = supercell_atoms.copy()
        plane_assignments = {}
        plane_id = 1
        
        # Repeat until all atoms are assigned
        while unassigned_atoms:
            # Select first unassigned atom as seed
            seed_atom = unassigned_atoms[0]
            
            # Calculate D value from seed atom coordinates to generate plane
            d_value = self._calculate_d_value_from_atom(seed_atom, plane_equation)
            current_plane = self._create_plane_with_d_value(plane_equation, d_value, plane_id)
            
            # Find all atoms belonging to this plane
            atoms_on_plane = self._find_atoms_on_current_plane(
                unassigned_atoms, current_plane, tolerance
            )
            
            # Assign atoms to plane
            if atoms_on_plane:
                plane_assignments[plane_id] = {
                    'plane_info': current_plane,
                    'atoms': atoms_on_plane,
                    'atom_count': len(atoms_on_plane)
                }
                
                # Remove assigned atoms from unassigned list
                for atom in atoms_on_plane:
                    if atom in unassigned_atoms:
                        unassigned_atoms.remove(atom)
                
                print(f"   SUCCESS: Plane {plane_id}: {len(atoms_on_plane)} atoms assigned (D={d_value:.4f})")
                plane_id += 1
            else:
                # Exception case: even seed atom doesn't belong to plane
                print(f"   WARNING: Seed atom doesn't belong to plane: {seed_atom.get('atom_id', 'Unknown')}")
                unassigned_atoms.remove(seed_atom)
        
        print(f"SUCCESS: Total {plane_id-1} planes generated")
        return plane_assignments
    
    def _calculate_d_value_from_atom(self, atom, plane_equation):
        """
        Calculate D value of plane equation from atom coordinates.
        
        Args:
            atom (dict): Atom information
            plane_equation (dict): Basic plane equation information
            
        Returns:
            float: Calculated D value
        """
        coeffs = plane_equation['coefficients']
        a, b, c = coeffs['a'], coeffs['b'], coeffs['c']
        
        x = atom['cart_x']
        y = atom['cart_y']
        z = atom['cart_z']
        
        # From ax + by + cz + d = 0, d = -(ax + by + cz)
        d_value = -(a*x + b*y + c*z)
        return d_value
    
    def _create_plane_with_d_value(self, base_plane_equation, d_value, plane_id):
        """
        Generate new plane information with specific D value.
        
        Args:
            base_plane_equation (dict): Basic plane equation information
            d_value (float): New D value
            plane_id (int): Plane ID
            
        Returns:
            dict: New plane information
        """
        coeffs = base_plane_equation['coefficients']
        a, b, c = coeffs['a'], coeffs['b'], coeffs['c']
        
        # Generate plane equation string
        equation_string = self._format_plane_equation(a, b, c, d_value)
        
        return {
            'plane_id': plane_id,
            'coefficients': {
                'a': a,
                'b': b,
                'c': c,
                'd': d_value
            },
            'equation_string': equation_string,
            'd_value': d_value
        }
    
    def _find_atoms_on_current_plane(self, atoms, plane_info, tolerance):
        """
        Find all atoms belonging to current plane.
        
        Args:
            atoms (list): Atoms to examine
            plane_info (dict): Plane information
            tolerance (float): Tolerance
            
        Returns:
            list: Atoms belonging to plane
        """
        coeffs = plane_info['coefficients']
        a, b, c, d = coeffs['a'], coeffs['b'], coeffs['c'], coeffs['d']
        
        atoms_on_plane = []
        for atom in atoms:
            x = atom['cart_x']
            y = atom['cart_y']
            z = atom['cart_z']
            
            # Calculate distance from plane
            distance = abs(a*x + b*y + c*z + d)
            
            if distance <= tolerance:
                atoms_on_plane.append(atom)
        
        return atoms_on_plane
    
    def _format_plane_equation(self, a, b, c, d):
        """
        Format plane equation coefficients as string.
        
        Args:
            a, b, c, d (float): Plane equation coefficients
            
        Returns:
            str: Formatted plane equation string
        """
        # Utilize existing formatting method of PlaneCalculator
        return PlaneCalculator.format_plane_equation((a, b, c, d))
    
    def _calculate_projection_range(self, atoms, plane_equation):
        """
        Calculate projection range of all atoms (for existing code compatibility).
        
        Args:
            atoms (list): List of atom information
            plane_equation (dict): Plane equation information
            
        Returns:
            dict: Projection range information
        """
        coeffs = plane_equation['coefficients']
        a, b, c = coeffs['a'], coeffs['b'], coeffs['c']
        
        projection_values = []
        for atom in atoms:
            x, y, z = atom['cart_x'], atom['cart_y'], atom['cart_z']
            projection = -(a*x + b*y + c*z)
            projection_values.append(projection)
        
        return {
            'min': min(projection_values),
            'max': max(projection_values)
        }
    
    # Existing methods (maintaining compatibility)
    def calculate_atom_projection_value(self, atom, plane_equation):
        """
        Calculate projection value of individual atom with respect to plane.
        
        Args:
            atom (dict): Atom information
            plane_equation (dict): Plane equation information
            
        Returns:
            float: Projection value
        """
        coeffs = plane_equation['coefficients']
        a, b, c = coeffs['a'], coeffs['b'], coeffs['c']
        
        x = atom['cart_x']
        y = atom['cart_y']
        z = atom['cart_z']
        
        # Calculate projection value
        projection = -(a*x + b*y + c*z)
        
        return projection 