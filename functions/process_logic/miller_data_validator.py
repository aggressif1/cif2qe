"""
Phase3 Miller index based data validation related functions
"""

class MillerDataValidator:
    """Phase3 Miller index based data validation class"""
    
    def __init__(self):
        """Initialize Phase3MillerDataValidator"""
        pass
    
    def validate_input_data(self, step3_data):
        """
        Validate Step3 data
        
        Args:
            step3_data (dict): Step3 data
            
        Returns:
            bool: Validation success status
        """
        if not step3_data:
            print("ERROR: No Step3 data available.")
            return False
        
        return True
    
    def validate_miller_data_availability(self, unit_cell_group, expanded_supercell):
        """
        Validate Miller index based data availability
        
        Args:
            unit_cell_group (dict): Unit cell group
            expanded_supercell (dict): Expanded supercell
            
        Returns:
            tuple: (bool, str) - (validation result, message)
        """
        if not unit_cell_group or not expanded_supercell:
            return False, "Miller index based unit cell or supercell has not been generated."
        
        return True, "Miller index based data is properly prepared." 