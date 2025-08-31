"""
Miller indices display utility
"""

class MillerDisplayHelper:
    """Utility class responsible for Miller indices display"""
    
    @staticmethod
    def display_miller_indices(data, data_paths=None, default_value=0, prefix="   - "):
        """
        Extract and display Miller indices from various data structures
        
        Args:
            data (dict): Data containing Miller indices
            data_paths (list): List of Miller indices access paths (in priority order)
            default_value: Default value (0 or 'N/A')
            prefix (str): Output prefix
            
        Returns:
            bool: Whether Miller indices were found and displayed
        """
        miller_indices = None
        
        # Default paths (in priority order)
        if data_paths is None:
            data_paths = [
                'miller_indices',
                ('phase3_data', 'miller_indices'),
                ('unit_cell_group', 'miller_indices'),
                ('expanded_supercell', 'miller_indices')
            ]
        
        # Search for Miller indices by path
        for path in data_paths:
            if isinstance(path, str):
                # Simple key access
                miller_indices = data.get(path)
            elif isinstance(path, tuple):
                # Nested key access
                temp_data = data
                for key in path:
                    if isinstance(temp_data, dict) and key in temp_data:
                        temp_data = temp_data[key]
                    else:
                        temp_data = None
                        break
                miller_indices = temp_data
            
            # Stop if Miller indices are found
            if miller_indices:
                break
        
        # Display Miller indices
        if miller_indices:
            h = miller_indices.get('h', default_value)
            k = miller_indices.get('k', default_value)
            l = miller_indices.get('l', default_value)
            print(f"{prefix}Miller indices: ({h} {k} {l})")
            return True
        
        return False
    
    @staticmethod
    def extract_miller_indices(data, data_paths=None, default_value=0):
        """
        Extract Miller indices only (without display)
        
        Args:
            data (dict): Data containing Miller indices
            data_paths (list): List of Miller indices access paths
            default_value: Default value
            
        Returns:
            tuple: (h, k, l) or None
        """
        miller_indices = None
        
        # Default paths
        if data_paths is None:
            data_paths = [
                'miller_indices',
                ('phase3_data', 'miller_indices'),
                ('unit_cell_group', 'miller_indices'),
                ('expanded_supercell', 'miller_indices')
            ]
        
        # Search for Miller indices by path
        for path in data_paths:
            if isinstance(path, str):
                miller_indices = data.get(path)
            elif isinstance(path, tuple):
                temp_data = data
                for key in path:
                    if isinstance(temp_data, dict) and key in temp_data:
                        temp_data = temp_data[key]
                    else:
                        temp_data = None
                        break
                miller_indices = temp_data
            
            if miller_indices:
                break
        
        if miller_indices:
            h = miller_indices.get('h', default_value)
            k = miller_indices.get('k', default_value)
            l = miller_indices.get('l', default_value)
            return (h, k, l)
        
        return None 