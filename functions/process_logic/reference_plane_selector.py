import csv
import os

class ReferencePlaneSelector:
    """Class for selecting reference planes and updating CSV files"""
    
    def __init__(self):
        """Initialize ReferencePlaneSelector"""
        self.reference_plane_info = None
    
    def select_reference_plane(self, plane_assignments):
        """
        Select the most central plane among planes with the highest atom count as the first reference plane.
        
        Args:
            plane_assignments (dict): Plane assignment results
            
        Returns:
            dict: Reference plane information or None
        """
        print("\nSelecting first reference plane...")
        
        try:
            assignments = plane_assignments['plane_assignments']
            
            if not assignments:
                print("ERROR: No assigned planes found.")
                return None
            
            # Step 1: Find planes with highest atom count
            max_atom_count = max(assignment['atom_count'] for assignment in assignments.values())
            max_atom_planes = [
                (plane_id, assignment) for plane_id, assignment in assignments.items()
                if assignment['atom_count'] == max_atom_count
            ]
            
            print(f"Maximum atom count: {max_atom_count}")
            print(f"Number of planes with maximum atom count: {len(max_atom_planes)}")
            
            # Step 2: Select middle Plane ID among planes with highest atom count
            plane_ids = [plane_id for plane_id, _ in max_atom_planes]
            
            # Extract and sort by number if plane_id is in "Plane X" format
            def extract_plane_number(plane_id):
                if isinstance(plane_id, str) and plane_id.startswith("Plane "):
                    return int(plane_id.split()[1])
                return plane_id
            
            plane_ids.sort(key=extract_plane_number)  # Sort by Plane ID order
            
            # Select middle value
            middle_index = len(plane_ids) // 2
            selected_plane_id = plane_ids[middle_index]
            selected_assignment = assignments[selected_plane_id]
            
            # Store reference plane information (simplified structure, excluding atom information)
            self.reference_plane_info = {
                'plane_id': selected_plane_id,
                'atom_count': selected_assignment['atom_count'],
                'plane_equation': selected_assignment['plane_info']['equation_string'],
                'd_value': float(selected_assignment['plane_info']['d_value'])  # Convert numpy type
            }
            
            print(f"SUCCESS: First reference plane selection completed:")
            print(f"   - Plane ID: {self.reference_plane_info['plane_id']}")
            print(f"   - D-value: {self.reference_plane_info['d_value']:.4f}")
            print(f"   - Atom count: {self.reference_plane_info['atom_count']}")
            
            # Output only brief summary information (excluding massive atom information)
            summary_info = {
                'plane_id': self.reference_plane_info['plane_id'],
                'd_value': self.reference_plane_info['d_value'],
                'atom_count': self.reference_plane_info['atom_count']
            }
            
            # Output in structured format
            print(f"   - Plane ID: {summary_info['plane_id']}")
            print(f"   - D-value: {summary_info['d_value']:.4f}")
            print(f"   - Atom count: {summary_info['atom_count']}")
            
            # Save first reference plane information to temporary file
            try:
                from ..utils.temporary_file_handler import TemporaryFileHandler
                temp_handler = TemporaryFileHandler()
                
                # Load existing Phase3 data (if any)
                existing_data = temp_handler.load_phase_results("phase3") or {}
                
                # Add first reference plane information
                existing_data['first_reference_plane'] = {
                    'plane_id': selected_plane_id,
                    'atom_count': selected_assignment['atom_count'],
                    'plane_equation': selected_assignment['plane_info']['equation_string'],
                    'd_value': float(selected_assignment['plane_info']['d_value'])
                }
                
                # Save to temporary file
                temp_handler.save_phase_results("phase3", existing_data)
                print(f"First reference plane information saved to temp file: {selected_plane_id}")
                
            except Exception as e:
                print(f"WARNING: Failed to save first reference plane information to temp file: {str(e)}")
            
            return summary_info
            
        except Exception as e:
            print(f"ERROR: Reference plane selection failed: {str(e)}")
            return None
    
    def update_csv_with_reference_plane(self, csv_file_path, plane_assignments):
        """
        Update CSV file to add reference plane information.
        
        Args:
            csv_file_path (str): CSV file path
            plane_assignments (dict): Plane assignment results
        """
        print(f"\nUpdating CSV file: {csv_file_path}")
        
        try:
            if not self.reference_plane_info:
                print("ERROR: No reference plane information available.")
                return False
            
            # Create coordinate set of atoms belonging to reference plane
            reference_plane_id = self.reference_plane_info['plane_id']
            assignments = plane_assignments['plane_assignments']
            
            if reference_plane_id not in assignments:
                print(f"ERROR: Cannot find assignment information for reference plane {reference_plane_id}.")
                return False
                
            reference_atoms = assignments[reference_plane_id]['atoms']
            reference_coords_set = set()
            
            for atom in reference_atoms:
                coord_key = (
                    round(atom['cart_x'], 6),
                    round(atom['cart_y'], 6), 
                    round(atom['cart_z'], 6)
                )
                reference_coords_set.add(coord_key)
            
            # Read CSV file
            if not os.path.exists(csv_file_path):
                print(f"ERROR: CSV file does not exist: {csv_file_path}")
                return False
            
            rows = []
            with open(csv_file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)
                rows.append(header)
                
                # Find required column indices in header
                try:
                    cart_x_idx = header.index('cart_x')
                    cart_y_idx = header.index('cart_y')
                    cart_z_idx = header.index('cart_z')
                    plane_eq_idx = header.index('plane_equation')
                    ref_plane_idx = header.index('reference_plane')
                    plane_type_idx = header.index('plane_type')
                except ValueError as e:
                    print(f"ERROR: Cannot find required columns: {e}")
                    return False
                
                # Process each row
                updated_count = 0
                for row in reader:
                    try:
                        # Check atom coordinates
                        cart_x = round(float(row[cart_x_idx]), 6)
                        cart_y = round(float(row[cart_y_idx]), 6)
                        cart_z = round(float(row[cart_z_idx]), 6)
                        coord_key = (cart_x, cart_y, cart_z)
                        
                        # Check if belongs to reference plane
                        if coord_key in reference_coords_set:
                            # Update only reference plane atoms (keep existing plane_equation)
                            row[ref_plane_idx] = "First reference plane"  # Set reference plane as "First reference plane"
                            row[plane_type_idx] = "Type A"
                            updated_count += 1
                        # Keep existing values for other atoms (do not overwrite)
                        
                        rows.append(row)
                        
                    except (ValueError, IndexError) as e:
                        print(f"WARNING: Error processing row: {e}")
                        rows.append(row)  # Keep original row
            
            # Save updated CSV file
            with open(csv_file_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerows(rows)
            
            print(f"SUCCESS: CSV file update completed:")
            print(f"   - First reference plane atoms: {updated_count}")
            print(f"   - Other atoms: {len(rows) - 1 - updated_count}")  # Excluding header
            
            return True
            
        except Exception as e:
            print(f"ERROR: CSV file update failed: {str(e)}")
            return False
    
    def get_reference_plane_info(self):
        """
        Get reference plane information.
        
        Returns:
            dict: Reference plane information
        """
        return self.reference_plane_info
    
    def display_reference_plane_summary(self):
        """
        Display reference plane summary.
        """
        if not self.reference_plane_info:
            print("No reference plane information available.")
            return
        
        print("\nReference plane summary:")
        print(f"   Plane ID: {self.reference_plane_info['plane_id']}")
        print(f"   D-value: {self.reference_plane_info['d_value']:.4f}")
        print(f"   Atom count: {self.reference_plane_info['atom_count']}")
        print(f"   Plane equation: {self.reference_plane_info['plane_equation']}")
    
    def update_csv_with_all_plane_info(self, csv_file_path, plane_assignments):
        """
        Update CSV file with all plane information.
        
        Args:
            csv_file_path (str): CSV file path
            plane_assignments (dict): Plane assignment results
        """
        print(f"\nUpdating CSV with all plane information: {csv_file_path}")
        
        try:
            assignments = plane_assignments['plane_assignments']
            
            if not assignments:
                print("ERROR: No plane assignment information available.")
                return False
            
            # Create atom coordinate mapping to plane information
            atom_to_plane_map = {}
            
            for plane_id, assignment in assignments.items():
                plane_equation = assignment['plane_info']['equation_string']
                atoms = assignment['atoms']
                
                for atom in atoms:
                    coord_key = (
                        round(atom['cart_x'], 6),
                        round(atom['cart_y'], 6),
                        round(atom['cart_z'], 6)
                    )
                    atom_to_plane_map[coord_key] = {
                        'plane_id': plane_id,
                        'plane_equation': plane_equation,
                        'reference_plane': 'Unknown',  # Set as Unknown initially
                        'plane_type': 'Unknown'
                    }
            
            # Read and update CSV file
            if not os.path.exists(csv_file_path):
                print(f"ERROR: CSV file does not exist: {csv_file_path}")
                return False
            
            rows = []
            with open(csv_file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)
                rows.append(header)
                
                # Find required column indices
                try:
                    cart_x_idx = header.index('cart_x')
                    cart_y_idx = header.index('cart_y')
                    cart_z_idx = header.index('cart_z')
                    plane_id_idx = header.index('plane_id')
                    plane_eq_idx = header.index('plane_equation')
                    ref_plane_idx = header.index('reference_plane')
                    plane_type_idx = header.index('plane_type')
                except ValueError as e:
                    print(f"ERROR: Cannot find required columns: {e}")
                    return False
                
                # Process each row
                updated_count = 0
                for row in reader:
                    try:
                        # Check atom coordinates
                        cart_x = round(float(row[cart_x_idx]), 6)
                        cart_y = round(float(row[cart_y_idx]), 6)
                        cart_z = round(float(row[cart_z_idx]), 6)
                        coord_key = (cart_x, cart_y, cart_z)
                        
                        # Update plane information if found
                        if coord_key in atom_to_plane_map:
                            plane_info = atom_to_plane_map[coord_key]
                            row[plane_id_idx] = plane_info['plane_id']
                            row[plane_eq_idx] = plane_info['plane_equation']
                            row[ref_plane_idx] = plane_info['reference_plane']
                            row[plane_type_idx] = plane_info['plane_type']
                            updated_count += 1
                        
                        rows.append(row)
                        
                    except (ValueError, IndexError) as e:
                        print(f"WARNING: Error processing row: {e}")
                        rows.append(row)  # Keep original row
            
            # Save updated CSV file
            with open(csv_file_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerows(rows)
            
            print(f"SUCCESS: All plane information update completed:")
            print(f"   - Updated atoms: {updated_count}")
            print(f"   - Total rows: {len(rows) - 1}")  # Excluding header
            
            return True
            
        except Exception as e:
            print(f"ERROR: All plane information update failed: {str(e)}")
            return False 