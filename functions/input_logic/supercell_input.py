class SupercellInput:
    """Class responsible for handling supercell input processing"""
    
    def get_supercell_size(self):
        """
        Returns the supercell size. Currently fixed at 6x6x6.
        
        Returns:
            tuple: (nx, ny, nz) supercell size
        """
        return (6, 6, 6)
    
    @staticmethod
    def get_vacuum_size():
        """
        Get vacuum region size from user input.
        
        Returns:
            float: Vacuum region size (Å)
        """
        while True:
            try:
                print("\nEnter the vacuum region size (Å) (0 = no vacuum):")
                vacuum = float(input().strip())
                
                # Check for negative values
                if vacuum < 0:
                    print("\nError: Vacuum size cannot be negative.")
                    continue
                
                return vacuum
                
            except ValueError:
                print("\nError: Please enter a valid number.") 