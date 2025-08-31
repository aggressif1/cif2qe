class ReferencePlaneCalculator:
    """Class for setting and managing the first reference plane"""
    
    def __init__(self):
        self.reference_plane = None
    
    def set_reference_plane(self, plane_counts):
        """
        Set the first reference plane based on atom count information for each plane.
        
        Args:
            plane_counts (dict): Atom count information for each plane
                                {plane_id: {'count': int, 'projection_value': str, 'equation': str}}
        
        Returns:
            dict: First reference plane information
        """
        if not plane_counts:
            return None
        
        # Calculate total number of planes
        total_planes = len(plane_counts)
        
        # Set middle plane as first reference plane
        middle_plane_number = (total_planes + 1) // 2  # Exact middle for odd number, front middle for even number
        reference_plane_id = f"Plane {middle_plane_number}"
        
        if reference_plane_id in plane_counts:
            ref_info = plane_counts[reference_plane_id]
            
            # Generate first reference plane information
            self.reference_plane = {
                'plane_id': reference_plane_id,
                'plane_number': middle_plane_number,
                'atom_count': ref_info['count'],
                'projection_value': ref_info['projection_value'],
                'equation': ref_info['equation'],
                'total_planes': total_planes
            }
            
            return self.reference_plane
        
        return None
    
    def display_reference_plane_info(self):
        """Display first reference plane information on screen."""
        if not self.reference_plane:
            print("First reference plane has not been set.")
            return
        
        ref = self.reference_plane
        print(f"\n[ First Reference Plane Settings ]")
        print(f"Total number of planes: {ref['total_planes']}")
        print(f"First reference plane: {ref['plane_id']}")
        print(f"  - Number of atoms: {ref['atom_count']}")
        print(f"  - Projection value: {ref['projection_value']}")
        print(f"  - Plane equation: {ref['equation']}")
    
    def get_reference_plane(self):
        """Return first reference plane information."""
        return self.reference_plane
    
    def is_reference_plane_set(self):
        """Check if first reference plane has been set."""
        return self.reference_plane is not None
    
    def get_reference_plane_number(self):
        """Return the number of the first reference plane."""
        if self.reference_plane:
            return self.reference_plane['plane_number']
        return None
    
    def get_relative_position_from_reference(self, plane_number):
        """
        Calculate how many planes away a specific plane is from the first reference plane.
        
        Args:
            plane_number (int): Target plane number
        
        Returns:
            int: Relative position from first reference plane (negative for below, positive for above)
        """
        if not self.reference_plane:
            return None
        
        return plane_number - self.reference_plane['plane_number'] 