"""
Phase3 Miller index based statistics generation related functions
"""
from ..utils.element_counter import ElementCounter

class MillerStatisticsGenerator:
    """Phase3 Miller index based statistics generation class"""
    
    def __init__(self):
        """Initialize Phase3MillerStatisticsGenerator"""
        pass
    
    def generate_comparison_statistics(self, unit_cell_group, expanded_supercell):
        """
        Generate comparison statistics.
        
        Args:
            unit_cell_group (dict): Unit cell group
            expanded_supercell (dict): Expanded supercell
            
        Returns:
            dict: Comparison statistics
        """
        unit_atoms = unit_cell_group.get('atoms', [])
        super_atoms = expanded_supercell.get('atoms', [])
        supercell_size = expanded_supercell.get('size', {})
        
        unit_elements = ElementCounter.count_elements(unit_atoms)
        super_elements = ElementCounter.count_elements(super_atoms)
        
        return {
            'unit_cell_total': len(unit_atoms),
            'supercell_total': len(super_atoms),
            'unit_cell_elements': unit_elements,
            'supercell_elements': super_elements,
            'supercell_size': supercell_size,
            'expansion_factor': supercell_size.get('nx', 1) * supercell_size.get('ny', 1) * supercell_size.get('nz', 1)
        }
    
 