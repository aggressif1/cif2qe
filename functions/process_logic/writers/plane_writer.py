import os
import csv
import numpy as np
from .supercell_writer import SupercellWriter

class PlaneWriter(SupercellWriter):
    """Class for outputting plane equation results"""
    
    def __init__(self):
        super().__init__()
        self.plane_count = 0  # Counter to track plane numbers
        self.plane_atoms_file = None  # CSV file writer
        self.writer = None
    
    def initialize_plane_csv(self, output_dir):
        """
        Initialize CSV file to store plane analysis results.
        
        Args:
            output_dir (str): output directory path
        """
        filepath = os.path.join(output_dir, 'plane_analysis.csv')
        self.plane_atoms_file = open(filepath, 'w', newline='', encoding='utf-8')
        self.writer = csv.writer(self.plane_atoms_file)
        self.writer.writerow(['Plane_ID', 'D_Value', 'Element', 'Cart_X', 'Cart_Y', 'Cart_Z'])
        print(f"\nPlane analysis results will be saved to plane_analysis.csv file.")
    
    def close_plane_csv(self):
        """Close CSV file."""
        if self.plane_atoms_file:
            self.plane_atoms_file.close()
            print("\nPlane analysis completed.")
    
    def display_plane_equation(self, h, k, l, coefficients):
        """
        Display plane equation.
        
        Args:
            h, k, l (int): Miller indices
            coefficients (tuple): (a, b, c, d) coefficients of plane equation
        """
        print("\n[ Plane Equation Results ]")
        print("[ Miller Indices ] ({} {} {})".format(h, k, l))
        print("[ Plane Equation ]")
        
        a, b, c, d = coefficients
        print(f"{a:8.4f}x + {b:8.4f}y + {c:8.4f}z + {d:8.4f} = 0")
    
    def save_atoms_on_plane(self, d_value, atoms):
        """
        Save atoms belonging to a plane with specific D value to CSV file.
        
        Args:
            d_value (float): D value of plane equation
            atoms (list): list of atoms belonging to the plane
        """
        if not atoms or not self.writer:
            return
            
        self.plane_count += 1
        plane_id = f"Plane {self.plane_count}"
        
        # Save atom information to CSV file
        for atom in atoms:
            self.writer.writerow([
                plane_id,
                f"{d_value:.4f}",
                atom['element'],
                f"{atom['cart_x']:.6f}",
                f"{atom['cart_y']:.6f}",
                f"{atom['cart_z']:.6f}"
            ])
        
        # Display progress (every 10 planes)
        if self.plane_count % 10 == 0:
            print(f"Processing... {self.plane_count} planes processed")

    @staticmethod
    def display_atoms_on_plane(d_value, atoms):
        """
        Display atoms belonging to a plane with specific D value.
        
        Args:
            d_value (float): D value of plane equation
            atoms (list): list of atoms belonging to the plane. Each atom has format {'element': str, 'cart_x': float, 'cart_y': float, 'cart_z': float}
        """
        print(f"\n[ Atoms on D = {d_value:.4f} Plane ]")
        print(f"Total {len(atoms)} atoms found.")
        
        # Classify by atom type
        atoms_by_element = {}
        for atom in atoms:
            element = atom['element']
            if element not in atoms_by_element:
                atoms_by_element[element] = []
            atoms_by_element[element].append(atom)
        
        # Output by atom type
        for element, atom_list in atoms_by_element.items():
            print(f"\n{element} atoms ({len(atom_list)}):")
            for atom in atom_list:
                print(f"  Position: [{atom['cart_x']:8.4f}, {atom['cart_y']:8.4f}, {atom['cart_z']:8.4f}]") 

    def display_plane_analysis_results(self, miller_indices, plane_equation, plane_assignments):
        """
        Display Step3 plane analysis results.
        
        Args:
            miller_indices (dict): Miller index information
            plane_equation (dict): basic plane equation information
            plane_assignments (dict): plane assignment results
        """
        print("\n" + "=" * 80)
        print("Plane Analysis Results")
        print("=" * 80)
        
        # Output Miller indices and basic plane equation
        h, k, l = miller_indices['h'], miller_indices['k'], miller_indices['l']
        print(f"\nMiller indices: ({h} {k} {l})")
        print(f"Basic plane equation: {plane_equation['equation_string']}")
        
        # Output plane assignment results
        assignments = plane_assignments['plane_assignments']
        print(f"\nNumber of generated planes: {len(assignments)}")
        
        if not assignments:
            print("ERROR: No planes with assigned atoms.")
            return
        
        # Output plane information (sorted by projection value - closest to origin first)
        # Sort by projection value and assign new Plane IDs
        sorted_assignments = []
        for plane_id, assignment in assignments.items():
            d_value = assignment['plane_info']['d_value']
            sorted_assignments.append((abs(d_value), plane_id, assignment))  # Sort by absolute value
        
        sorted_assignments.sort(key=lambda x: x[0])  # Ascending order by absolute value
        
        # Create new plane_assignments dictionary (change keys to display IDs)
        new_assignments = {}
        
        # Plane information table header (expanded width)
        print("\n" + "=" * 120)
        print(f"{'Plane ID':<10} {'Plane Equation':<60} {'Atom Count':<8} {'Projection Value':<12}")
        print("=" * 120)
        
        total_atoms_assigned = 0
        for new_plane_id, (_, original_plane_id, assignment) in enumerate(sorted_assignments, 1):
            plane_info = assignment['plane_info']
            atom_count = assignment['atom_count']
            total_atoms_assigned += atom_count
            
            # Plane equation string (display full without truncation)
            equation_str = plane_info['equation_string']
            
            # Projection value (D value)
            d_value = plane_info['d_value']
            
            print(f"Plane {new_plane_id:<3} {equation_str:<60} {atom_count:<8} {d_value:<12.4f}")
            
            # Update assignment data with new plane_id
            new_plane_key = f"Plane {new_plane_id}"
            new_assignments[new_plane_key] = assignment.copy()
            # Update plane_id in plane_info as well
            new_assignments[new_plane_key]['plane_info'] = assignment['plane_info'].copy()
            new_assignments[new_plane_key]['plane_info']['plane_id'] = new_plane_id
        
        print("=" * 120)
        print(f"Total assigned atoms: {total_atoms_assigned}")
        
        # Update plane_assignments data (unified with display IDs)
        plane_assignments['plane_assignments'] = new_assignments
        
        # Projection value range information
        proj_range = plane_assignments['projection_range']
        print(f"\nProjection value range: {proj_range['min']:.4f} ~ {proj_range['max']:.4f}")
        # Output plane spacing (check type as it might be string)
        plane_spacing = plane_assignments['plane_spacing']
        if isinstance(plane_spacing, (int, float)):
            print(f"Plane spacing: {plane_spacing:.4f}")
        else:
            print(f"Plane spacing: {plane_spacing}")
        print(f"Tolerance: {plane_assignments['tolerance']:.4f}")
        
        # Output element statistics (using new assignments)
        self._display_element_statistics(new_assignments)
    
    def _display_element_statistics(self, assignments):
        """
        Display element distribution statistics across planes.
        
        Args:
            assignments (dict): plane assignment results
        """
        print(f"\nElement Distribution Across Planes:")
        
        # Collect element statistics
        element_stats = {}
        
        for plane_id, assignment in assignments.items():
            atoms = assignment['atoms']
            
            for atom in atoms:
                element = atom.get('element', atom.get('atom_type', 'Unknown'))
                
                if element not in element_stats:
                    element_stats[element] = {'total': 0, 'planes': set()}
                
                element_stats[element]['total'] += 1
                element_stats[element]['planes'].add(plane_id)
        
        # Output statistics
        print(f"{'Element':<8} {'Total Atoms':<12} {'Plane Count':<12}")
        print("-" * 32)
        
        for element in sorted(element_stats.keys()):
            total_atoms = element_stats[element]['total']
            plane_count = len(element_stats[element]['planes'])
            print(f"{element:<8} {total_atoms:<12} {plane_count:<12}")
    
    def display_detailed_plane_info(self, plane_id, assignment, show_atoms=False):
        """
        Display detailed information for specific plane.
        
        Args:
            plane_id (int): plane ID
            assignment (dict): plane assignment information
            show_atoms (bool): whether to show detailed atom information
        """
        plane_info = assignment['plane_info']
        atoms = assignment['atoms']
        
        print(f"\n{'='*60}")
        print(f"Plane {plane_id} Detailed Information")
        print(f"{'='*60}")
        
        print(f"Plane equation: {plane_info['equation_string']}")
        print(f"D value: {plane_info['d_value']:.4f}")
        print(f"Number of atoms: {len(atoms)}")
        
        if show_atoms and atoms:
            print(f"\nAtom list:")
            print(f"{'No.':<6} {'Element':<6} {'X Coord':<12} {'Y Coord':<12} {'Z Coord':<12}")
            print("-" * 54)
            
            for i, atom in enumerate(atoms, 1):
                element = atom.get('element', atom.get('atom_type', 'Unknown'))
                x, y, z = atom['cart_x'], atom['cart_y'], atom['cart_z']
                print(f"{i:<6} {element:<6} {x:<12.4f} {y:<12.4f} {z:<12.4f}") 