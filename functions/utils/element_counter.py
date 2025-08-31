"""
Element counting utility
"""

class ElementCounter:
    """Utility class responsible for counting elements by type"""
    
    @staticmethod
    def count_elements(atoms):
        """
        Calculate the count of each element type from atoms.
        
        Args:
            atoms (list): List of atom information
            
        Returns:
            dict: Count by element type
        """
        element_counts = {}
        for atom in atoms:
            element = atom.get('element', atom.get('atom_type', 'Unknown'))
            element_counts[element] = element_counts.get(element, 0) + 1
        return element_counts 