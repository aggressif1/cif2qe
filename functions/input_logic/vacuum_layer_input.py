"""
Vacuum layer configuration input related functions
"""

class VacuumLayerInput:
    """Class responsible for vacuum layer configuration input"""
    
    def __init__(self):
        """Initialize VacuumLayerInput"""
        pass
    
    def get_vacuum_configuration(self):
        """
        Get vacuum layer configuration from user input.
        
        Returns:
            dict: Vacuum layer configuration information
                - coordinate_system: 'cartesian' or 'user_defined'
                - vacuum_layers: Vacuum layer thickness for each direction
        """
        try:
            print("\n🌌 Vacuum Layer Addition Configuration")
            print("=" * 50)
            
            # Select coordinate system
            coordinate_system = self._get_coordinate_system()
            if not coordinate_system:
                return None
            
            # Input vacuum layer thickness
            vacuum_layers = self._get_vacuum_layer_thickness(coordinate_system)
            if not vacuum_layers:
                return None
            
            # Configuration verification
            config = {
                'coordinate_system': coordinate_system,
                'vacuum_layers': vacuum_layers
            }
            
            self._display_configuration_summary(config)
            
            # User confirmation
            confirm = input("\nDo you want to add vacuum layers with the above settings? (y/n): ").strip().lower()
            if confirm not in ['y', 'yes']:
                print("❌ Vacuum layer addition cancelled.")
                return None
            
            return config
            
        except Exception as e:
            print(f"❌ Error during vacuum layer configuration input: {str(e)}")
            return None
    
    def _get_coordinate_system(self):
        """Select coordinate system"""
        print("\n📐 Select vacuum layer addition direction:")
        print("1. Cartesian coordinate system (+x, +y, +z directions)")
        print("2. User-defined cell vectors (+a_user, +b_user, +c_user directions)")
        
        while True:
            try:
                choice = input("Selection (1 or 2): ").strip()
                
                if choice == '1':
                    print("✅ Cartesian coordinate system (+x, +y, +z) selected")
                    return 'cartesian'
                elif choice == '2':
                    print("✅ User-defined cell vectors (+a_user, +b_user, +c_user) selected")
                    return 'user_defined'
                else:
                    print("❌ Invalid selection. Please enter 1 or 2.")
                    
            except KeyboardInterrupt:
                print("\n❌ Cancelled by user.")
                return None
    
    def _get_vacuum_layer_thickness(self, coordinate_system):
        """Input vacuum layer thickness"""
        if coordinate_system == 'cartesian':
            directions = ['x', 'y', 'z']
            print("\n📏 Enter vacuum layer thickness for each Cartesian coordinate direction (in Å):")
        else:
            directions = ['a_user', 'b_user', 'c_user']
            print("\n📏 Enter vacuum layer thickness for each user-defined cell vector direction (in Å):")
        
        vacuum_layers = {}
        
        for direction in directions:
            while True:
                try:
                    thickness_str = input(f"   +{direction} direction vacuum layer thickness (≥0, default: 0): ").strip()
                    
                    if not thickness_str:
                        thickness = 0.0
                    else:
                        thickness = float(thickness_str)
                    
                    if thickness < 0:
                        print("❌ Vacuum layer thickness must be 0 or greater.")
                        continue
                    
                    vacuum_layers[direction] = thickness
                    print(f"   ✅ +{direction} direction: {thickness:.2f} Å")
                    break
                    
                except ValueError:
                    print("❌ Please enter a valid number.")
                except KeyboardInterrupt:
                    print("\n❌ Cancelled by user.")
                    return None
        
        # Check if all directions are 0
        total_vacuum = sum(vacuum_layers.values())
        if total_vacuum == 0:
            print("⚠️  Vacuum layer thickness is 0 in all directions.")
            confirm = input("Do you want to continue without adding vacuum layers? (y/n): ").strip().lower()
            if confirm not in ['y', 'yes']:
                return None
        
        return vacuum_layers
    
    def _display_configuration_summary(self, config):
        """Display configuration summary"""
        print("\n📋 Vacuum Layer Addition Configuration Summary:")
        print("-" * 30)
        
        coord_system = config['coordinate_system']
        if coord_system == 'cartesian':
            print("Coordinate system: Cartesian coordinate system (+x, +y, +z)")
        else:
            print("Coordinate system: User-defined cell vectors (+a_user, +b_user, +c_user)")
        
        print("Vacuum layer thickness:")
        for direction, thickness in config['vacuum_layers'].items():
            print(f"   +{direction}: {thickness:.2f} Å")
        
        total_vacuum = sum(config['vacuum_layers'].values())
        print(f"Total vacuum layer volume increase: {total_vacuum:.2f} Å (linear sum)") 