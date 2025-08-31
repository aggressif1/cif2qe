import os
from ..process_logic.bond_calculator import BondCalculator
from .display_manager_base import BaseDisplayManager

class ConsoleWriter(BaseDisplayManager):
    """Class responsible for console output"""
    
    @staticmethod
    def write_cell_info(cell_info):
        """
        Output unit cell information.
        
        Args:
            cell_info (dict): Unit cell information
                {
                    'params': dict,  # Unit cell parameters
                    'vectors': dict, # Cell vectors
                    'volume': float, # Volume
                    'system': str    # Crystal system
                }
        """
        print("\n============================================================")
        print("Unit Cell Information")
        print("============================================================")
        
        # Output lattice constants and angles
        print("\n[ Lattice Constants ]")
        print(f"a = {cell_info['params']['a']:.4f} Å")
        print(f"b = {cell_info['params']['b']:.4f} Å")
        print(f"c = {cell_info['params']['c']:.4f} Å")
        
        print("\n[ Lattice Angles ]")
        print(f"α = {cell_info['params']['alpha']:.4f}°")
        print(f"β = {cell_info['params']['beta']:.4f}°")
        print(f"γ = {cell_info['params']['gamma']:.4f}°")
        
        # Output cell vectors
        print("\n[ Cell Vectors ]")
        for vector_name in ['a', 'b', 'c']:
            vector = cell_info['vectors'][vector_name]
            print(f"{vector_name} = [{vector[0]:8.4f}, {vector[1]:8.4f}, {vector[2]:8.4f}] Å")
        
        # Output volume and crystal system
        print(f"\n[ Volume ] {cell_info['volume']:.4f} Å³")
        print(f"[ Crystal System ] {cell_info['system']}")
    
    @staticmethod
    def show_distance_analysis_and_get_bond_range(bond_calc, atoms, cell_vectors=None):
        """
        Show actual interatomic distances and get bond distance range from user input.
        
        Args:
            bond_calc: Already created BondCalculator instance
            atoms (list): List of atomic information
            cell_vectors (dict): Cell vector information
            
        Returns:
            tuple: (minimum bond distance, maximum bond distance) (Å)
        """
        print("\n[ Bond Analysis Configuration ]")
        
        # Create extended atom list to calculate actual interatomic distances
        if cell_vectors:
            from ..process_logic.periodicity_calculator import PeriodicityCalculator
            periodicity_calc = PeriodicityCalculator(cell_vectors)
            extended_atoms = periodicity_calc.create_periodic_images(atoms, range_a=1, range_b=1, range_c=1)  # Nearest images only
        else:
            extended_atoms = atoms
        
        bond_types = {}
        
        # Calculate distances between all atom pairs and classify by bond type
        for i, atom1 in enumerate(atoms):
            for atom2 in extended_atoms:
                if not (atom2.get('is_original', False) and atom2.get('original_index', -1) == i):
                    distance = bond_calc.calculate_distance(atom1, atom2)
                    if distance > 0.1:  # Exclude too close distances (same atom)
                        # Normalize element pairs (sort alphabetically)
                        elements = sorted([atom1['element'], atom2['element']])
                        bond_type = f"{elements[0]}-{elements[1]}"
                        
                        if bond_type not in bond_types:
                            bond_types[bond_type] = []
                        bond_types[bond_type].append(distance)
        
        # Sort distances for each bond type
        for bond_type in bond_types:
            bond_types[bond_type].sort()
        
        # Display bond type information
        print("[ Bond Types Found in CIF File ]")
        print("Bond Type    |  Count  |  Min Distance (Å)")
        print("-" * 35)
        
        bond_options = {}
        option_num = 1
        for bond_type, distances in bond_types.items():
            min_dist = min(distances)
            bond_options[option_num] = bond_type
            
            print(f"{option_num:2d}. {bond_type:8s}  |  {len(distances):4d}   |  {min_dist:7.3f}")
            option_num += 1
        
        # Bond type selection
        if not bond_options:
            print("\nBond information not found.")
            return 0.5, 3.0  # Return default values
        
        print(f"\nSelect bond type to analyze (1-{len(bond_options)}):")
        print("0. All bonds (user-defined range)")
        
        bond_type = None
        while True:
            try:
                choice = input("Selection: ").strip()
                choice_num = int(choice)
                
                if choice_num == 0:
                    bond_type = None
                    break
                elif 1 <= choice_num <= len(bond_options):
                    bond_type = bond_options[choice_num]
                    selected_distances = bond_types[bond_type]
                    print(f"\nSelected bond: {bond_type}")
                    print(f"Distance range: {min(selected_distances):.3f}Å ~ {max(selected_distances):.3f}Å")
                    break
                else:
                    print(f"Please enter a number between 1 and {len(bond_options)} or 0.")
            except ValueError:
                print("Please enter a valid number.")
        
        # Set distance range based on selected bond
        if bond_type:
            selected_distances = bond_types[bond_type]
            default_min = min(selected_distances)
            default_max = max(selected_distances)
            
            print(f"\nSet distance range for {bond_type} bond:")
        else:
            default_min = 0.5
            default_max = 3.0
            print("\nSet bond distance range:")
        
        # Input minimum distance
        while True:
            try:
                min_dist = input(f"Minimum bond distance (Å, default: {default_min:.3f}): ").strip()
                if not min_dist:
                    min_distance = default_min
                    break
                min_distance = float(min_dist)
                if min_distance >= 0:
                    break
                else:
                    print("Please enter a value greater than or equal to 0.")
            except ValueError:
                print("Please enter a valid number.")
        
        # Input maximum distance
        while True:
            try:
                max_dist = input(f"Maximum bond distance (Å, default: {default_max:.3f}): ").strip()
                if not max_dist:
                    max_distance = default_max
                    break
                max_distance = float(max_dist)
                if max_distance > min_distance:
                    break
                else:
                    print(f"Please enter a value greater than the minimum bond distance ({min_distance:.3f}Å).")
            except ValueError:
                print("Please enter a valid number.")
        
        if bond_type:
            print(f"Set {bond_type} bond distance range: {min_distance:.3f}Å ~ {max_distance:.3f}Å")
        else:
            print(f"Set bond distance range: {min_distance:.3f}Å ~ {max_distance:.3f}Å")
        
        return min_distance, max_distance, bond_type
    
    @staticmethod
    def write_atomic_positions(atoms, cell_vectors=None):
        """
        Output atomic position information.
        
        Args:
            atoms (list): List of atomic information
            cell_vectors (dict): Cell vector information (for periodic boundary conditions)
            
        Returns:
            dict: Bond information (optional)
        """
        if not atoms:
            print("\nNo atomic position information available.")
            return None
        
        print("\n============================================================")
        print("Atomic Position Information")
        print("============================================================")
        
        # Calculate element counts
        element_counts = {}
        for atom in atoms:
            element = atom['element']
            element_counts[element] = element_counts.get(element, 0) + 1
        
        # Output element counts
        print("\nElement counts:")
        print("-" * 20)
        for element, count in element_counts.items():
            print(f"{element}: {count}")
        
        # Calculate bond information (considering periodic boundary conditions)
        bond_calc = BondCalculator()
        
        # Get bond distance range from user (with actual distance display)
        min_bond_dist, max_bond_dist, bond_type = ConsoleWriter.show_distance_analysis_and_get_bond_range(bond_calc, atoms, cell_vectors)
        
        # Calculate bond information
        bonds = bond_calc.find_bonds(atoms, cell_vectors, min_bond_distance=min_bond_dist, max_bond_distance=max_bond_dist)
        
        # Output atomic positions and bond information
        print("\nAtomic positions and bond information:")
        print("-" * 75)
        print("No.   Element Fractional Coordinates         Cartesian Coordinates (Å)")
        print("-" * 75)
        
        for i, atom in enumerate(atoms, 1):
            # If fractional coordinates exist
            frac_coords = ""
            if all(key in atom for key in ['x', 'y', 'z']):
                frac_coords = f"({atom['x']:7.4f}, {atom['y']:7.4f}, {atom['z']:7.4f})"
            
            # If Cartesian coordinates exist
            cart_coords = ""
            if all(key in atom for key in ['cart_x', 'cart_y', 'cart_z']):
                cart_coords = f"({atom['cart_x']:8.4f}, {atom['cart_y']:8.4f}, {atom['cart_z']:8.4f})"
            
            print(f"{i:4d}  {atom['element']:<6s} {frac_coords:28s} {cart_coords}")
            
            # Output bond information
            bond_info_lines = bond_calc.format_bond_info(bonds, i-1)  # Index starts from 0
            for bond_line in bond_info_lines:
                print(f"      {bond_line}")
            print()  # Add blank line between each atom's information
        
        # Return bond information (for use in Step2)
        bond_info = {
            'bond_type': bond_type,
            'min_distance': min_bond_dist,
            'max_distance': max_bond_dist,
            'bonds': bonds
        }
        
        return bond_info
    
    @staticmethod
    def write_error(message):
        """
        Output error messages.

        Args:
            message (str): Error message
        """
        BaseDisplayManager.show_error(message)
    
    @staticmethod
    def write_section_header(title):
        """
        Output section headers.

        Args:
            title (str): Section title
        """
        BaseDisplayManager.show_section_header(title)
    
    @staticmethod
    def display_cif_files(cif_files):
        """
        Display list of found CIF files in hierarchical tree structure.

        Args:
            cif_files (dict): Dictionary of CIF file paths by folder
        """
        BaseDisplayManager.display_file_list(cif_files, "Found CIF Files")
        
 