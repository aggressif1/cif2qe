class MillerIndexInput:
    """Class for handling Miller index input"""
    
    def __init__(self):
        """Initialize MillerIndexInput"""
        pass
    
    def get_miller_indices(self, crystal_system=None):
        """
        Get Miller indices from user input.
        
        Args:
            crystal_system (str, optional): crystal system information
        
        Returns:
            dict: Miller index information {'h': int, 'k': int, 'l': int, 'm': int (hexagonal only)} or None (if failed)
        """
        print("\nMiller Index Input")
        print("=" * 50)
        
        # Check if hexagonal crystal system
        is_hexagonal = crystal_system and crystal_system.lower() == "hexagonal"
        
        try:
            if is_hexagonal:
                print("Hexagonal crystal system detected!")
                print("Please enter Miller-Bravais indices (h k l m)")
                print("Note: m must equal -(h+k+l)")
                print("Example: 1 0 -1 0 or 1 1 -2 0")
                input_prompt = "Miller-Bravais indices (h k l m): "
                expected_count = 4
            else:
                print("Enter Miller indices (e.g., 1 0 0 or 1 1 1)")
                print("Note: 0 0 0 is not valid.")
                input_prompt = "Miller indices (h k l): "
                expected_count = 3
            
            while True:
                user_input = input(input_prompt).strip()
                
                if user_input.lower() == 'exit':
                    print("Miller index input cancelled.")
                    return None
                
                # Parse input
                try:
                    indices = user_input.split()
                    if len(indices) != expected_count:
                        print(f"ERROR: Exactly {expected_count} values must be entered. Please try again.")
                        continue
                    
                    indices = list(map(int, indices))
                    
                    if is_hexagonal:
                        h, k, l, m = indices
                        
                        # Miller-Bravais validation: m = -(h+k+l)
                        expected_m = -(h + k + l)
                        if m != expected_m:
                            print(f"ERROR: Invalid Miller-Bravais indices. m should be {expected_m}, but got {m}")
                            print("Miller-Bravais rule: m = -(h+k+l)")
                            continue
                        
                        # Validity check
                        if h == 0 and k == 0 and l == 0:
                            print("ERROR: (0 0 0 m) is not valid. Please try again.")
                            continue
                        
                        # Input confirmation
                        print(f"\nEntered Miller-Bravais indices: ({h} {k} {l} {m})")
                        print(f"Equivalent Miller indices: ({h} {k} {l})")
                        confirm = input("Proceed with these indices? (y/n): ").strip().lower()
                        
                        if confirm in ['y', 'yes', '']:
                            miller_indices = {'h': h, 'k': k, 'l': l, 'm': m}
                            print(f"Miller-Bravais indices set: ({h} {k} {l} {m})")
                            return miller_indices
                        else:
                            print("Please enter again.")
                            continue
                    else:
                        h, k, l = indices
                        
                        # Validity check
                        if h == 0 and k == 0 and l == 0:
                            print("ERROR: (0 0 0) is not valid. Please try again.")
                            continue
                        
                        # Input confirmation
                        print(f"\nEntered Miller indices: ({h} {k} {l})")
                        confirm = input("Proceed with these indices? (y/n): ").strip().lower()
                        
                        if confirm in ['y', 'yes', '']:
                            miller_indices = {'h': h, 'k': k, 'l': l}
                            print(f"Miller indices set: ({h} {k} {l})")
                            return miller_indices
                        else:
                            print("Please enter again.")
                            continue
                        
                except ValueError:
                    print("ERROR: Only integers are allowed. Please try again.")
                    continue
                except Exception as e:
                    print(f"ERROR: Input processing error: {str(e)}")
                    continue
        
        except KeyboardInterrupt:
            print("\n\nWARNING: Miller index input interrupted.")
            return None
        except Exception as e:
            print(f"ERROR: Error during Miller index input: {str(e)}")
            return None
    
    def validate_miller_indices(self, h, k, l, m=None):
        """
        Validate Miller indices.
        
        Args:
            h, k, l (int): Miller indices
            m (int, optional): Miller-Bravais index (hexagonal only)
            
        Returns:
            bool: validity status
        """
        # All zeros case is invalid
        if h == 0 and k == 0 and l == 0:
            return False
        
        # Check if integers
        indices_to_check = [h, k, l]
        if m is not None:
            indices_to_check.append(m)
            # Miller-Bravais rule check: m = -(h+k+l)
            if m != -(h + k + l):
                return False
        
        if not all(isinstance(x, int) for x in indices_to_check):
            return False
        
        return True
    
    def format_miller_indices(self, miller_indices):
        """
        Format Miller indices as string.
        
        Args:
            miller_indices (dict): Miller index information
            
        Returns:
            str: formatted Miller index string
        """
        if not miller_indices:
            return "Invalid Miller indices"
        
        h = miller_indices.get('h', 0)
        k = miller_indices.get('k', 0)
        l = miller_indices.get('l', 0)
        
        # Miller-Bravais indices for hexagonal system
        if 'm' in miller_indices:
            m = miller_indices.get('m', 0)
            return f"({h} {k} {l} {m})"
        else:
            return f"({h} {k} {l})" 