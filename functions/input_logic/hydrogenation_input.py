"""
Hydrogenation configuration input related functions
"""

class HydrogenationInput:
    """Class for handling hydrogenation configuration input"""
    
    def __init__(self):
        """Initialize HydrogenationInput"""
        pass
    
    def get_hydrogenation_configuration(self):
        """
        Get hydrogenation configuration from user input.
        
        Returns:
            dict: hydrogenation configuration information
        """
        try:
            print("\nHydrogenation Process Configuration")
            print("=" * 50)
            
            while True:  # Loop for configuration re-entry
                # Select hydrogenation direction (based on user-defined cell vectors)
                direction_config = self._get_direction_configuration()
                if not direction_config:
                    return None
                
                # Input hydrogen count (only if direction setting is successful)
                hydrogen_count = self._get_hydrogen_count()
                if hydrogen_count is None:
                    return None
                
                # Configuration confirmation
                config = {
                    'direction_config': direction_config,
                    'hydrogen_count': hydrogen_count
                }
                
                self._display_configuration_summary(config)
                
                # User confirmation
                confirm = input("\nProceed with hydrogenation using the above settings? (y/n): ").strip().lower()
                if confirm in ['y', 'yes']:
                    return config
                elif confirm in ['n', 'no']:
                    print("Re-entering configuration.\n")
                    continue  # Start over from the beginning
                else:
                    print("ERROR: Please enter y or n.")
                    continue  # Only re-ask confirmation question
            
        except Exception as e:
            print(f"ERROR: Error during hydrogenation configuration input: {str(e)}")
            return None
    
    def _get_direction_configuration(self):
        """Configure hydrogenation directions based on user-defined cell vectors"""
        print("\nHydrogenation Direction Configuration (Based on User-defined Cell Vectors):")
        print("Configure hydrogenation ranges for positive and negative directions of each cell vector.")
        print("Minimum range: +1.0 Å (hydrogen placement starts at least 1.0 Å away from surface)")
        
        direction_config = {}
        
        # Each user-defined cell vector direction (extended cell vectors with vacuum layer)
        cell_vectors = ['a_user_extended_withvac', 'b_user_extended_withvac', 'c_user_extended_withvac']
        
        for vector_name in cell_vectors:
            print(f"\n{vector_name} direction hydrogenation configuration:")
            
            # Positive direction configuration
            positive_config = self._get_direction_range(vector_name, "positive direction (+)")
            if positive_config is None:
                return None
            
            # Negative direction configuration
            negative_config = self._get_direction_range(vector_name, "negative direction (-)")
            if negative_config is None:
                return None
            
            direction_config[vector_name] = {
                'positive': positive_config,
                'negative': negative_config
            }
        
        # Check if all directions are disabled
        all_disabled = True
        for vector_name, directions in direction_config.items():
            if directions['positive']['enabled'] or directions['negative']['enabled']:
                all_disabled = False
                break
        
        if all_disabled:
            print("\nWARNING: Hydrogenation is disabled in all directions.")
            print("There is no space to add hydrogen atoms.")
            
            while True:
                choice = input("\nWould you like to reconfigure? (y: reconfigure, n: exit program): ").strip().lower()
                if choice in ['y', 'yes']:
                    return self._get_direction_configuration()  # Recursive call to reconfigure
                elif choice in ['n', 'no']:
                    print("Ending hydrogenation process.")
                    return None
                else:
                    print("ERROR: Please enter y or n.")
        
        return direction_config
    
    def _get_direction_range(self, vector_name, direction_desc):
        """Configure hydrogenation range for specific direction"""
        print(f"   {direction_desc} hydrogenation configuration:")
        
        # Check if hydrogenation should be enabled
        while True:
            enable_str = input(f"     Enable {direction_desc} hydrogenation? (y/n): ").strip().lower()
            if enable_str in ['y', 'yes']:
                enabled = True
                break
            elif enable_str in ['n', 'no']:
                enabled = False
                print(f"     {direction_desc} hydrogenation disabled")
                return {'enabled': False}
            else:
                print("     ERROR: Please enter y or n.")
        
        if not enabled:
            return {'enabled': False}
        
        # Range configuration
        while True:
            try:
                min_range_str = input(f"     Minimum range (Å, default: 1.0): ").strip()
                if not min_range_str:
                    min_range = 1.0
                else:
                    min_range = float(min_range_str)
                
                if min_range < 1.0:
                    print("     ERROR: Minimum range must be at least 1.0 Å.")
                    continue
                
                max_range_str = input(f"     Maximum range (Å, default: 10.0): ").strip()
                if not max_range_str:
                    max_range = 10.0
                else:
                    max_range = float(max_range_str)
                
                if max_range <= min_range:
                    print("     ERROR: Maximum range must be greater than minimum range.")
                    continue
                
                print(f"     {direction_desc}: {min_range:.1f} ~ {max_range:.1f} Å")
                
                return {
                    'enabled': True,
                    'min_range': min_range,
                    'max_range': max_range
                }
                
            except ValueError:
                print("     ERROR: Please enter a valid number.")
            except KeyboardInterrupt:
                print("\nERROR: Cancelled by user.")
                return None
    
    def _get_hydrogen_count(self):
        """Input hydrogen count"""
        print("\nNumber of Hydrogen Atoms to Add:")
        
        while True:
            try:
                count_str = input("   Number of hydrogen atoms (1-100, default: 10): ").strip()
                
                if not count_str:
                    count = 10  # Default value
                else:
                    count = int(count_str)
                
                if count < 1:
                    print("ERROR: Number of hydrogen atoms must be at least 1.")
                    continue
                elif count > 100:
                    print("ERROR: Number of hydrogen atoms is limited to 100 or less.")
                    continue
                
                print(f"   Number of hydrogen atoms: {count}")
                return count
                
            except ValueError:
                print("ERROR: Please enter a valid integer.")
            except KeyboardInterrupt:
                print("\nERROR: Cancelled by user.")
                return None
    
    def _display_configuration_summary(self, config):
        """Display configuration summary (in compressed table format)"""
        print("\nHydrogenation Configuration Summary:")
        print("-" * 40)
        
        direction_config = config['direction_config']
        hydrogen_count = config['hydrogen_count']
        
        # Table header
        print("Direction               Positive    Negative")
        print("-" * 40)
        
        # Display each direction configuration in one line
        for vector_name, directions in direction_config.items():
            pos_config = directions['positive']
            neg_config = directions['negative']
            
            # Positive direction information
            if pos_config['enabled']:
                pos_info = f"{pos_config['min_range']:.0f}~{pos_config['max_range']:.0f}Å"
            else:
                pos_info = "Disabled"
            
            # Negative direction information
            if neg_config['enabled']:
                neg_info = f"{neg_config['min_range']:.0f}~{neg_config['max_range']:.0f}Å"
            else:
                neg_info = "Disabled"
            
            # Display abbreviated direction name
            short_name = vector_name.replace('_user_extended_withvac', '')
            print(f"{short_name:20s}    {pos_info:8s}    {neg_info:8s}")
        
        print(f"\nNumber of hydrogen atoms: {hydrogen_count}") 