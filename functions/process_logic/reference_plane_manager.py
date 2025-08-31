"""
Reference plane management functions separated from Phase3
"""
import csv
import os

class ReferencePlaneManager:
    """Phase3 dedicated reference plane management class"""
    
    def __init__(self):
        """Initialize ReferencePlaneManager"""
        pass
    
    def perform_iterative_reference_plane_analysis(self, plane_assignments, initial_reference_plane, supercell_atoms, cell_vectors, miller_indices, comparison_analyzer):
        """
        Perform iterative reference plane analysis.
        
        Args:
            plane_assignments: Plane assignment results
            initial_reference_plane: Initial reference plane
            supercell_atoms: Supercell atoms
            cell_vectors: Cell vectors
            miller_indices: Miller indices
            comparison_analyzer: Plane comparison analyzer (PlaneComparisonAnalyzer instance)
            
        Returns:
            dict: Analysis results
        """
        print("\nStarting iterative reference plane analysis")
        print("=" * 60)
        
        try:
            current_reference_plane = initial_reference_plane
            reference_count = 1
            max_iterations = 50  # Maximum iteration count
            
            for iteration in range(max_iterations):
                print(f"\nIteration {iteration + 1}: Reference plane {current_reference_plane['plane_id']}")
                
                # Perform plane comparison analysis with current reference plane
                from .comparison_analyzer import Phase3ComparisonAnalyzer
                phase3_comparison_analyzer = ComparisonAnalyzer()
                phase3_comparison_analyzer.test_plane_comparison_analysis(
                    plane_assignments, current_reference_plane, comparison_analyzer, supercell_atoms, miller_indices, reference_count
                )
                
                # Find next reference plane
                next_reference_plane = self.find_next_reference_plane(current_reference_plane)
                
                if not next_reference_plane:
                    print(f"SUCCESS: No more reference planes available. Analysis completed.")
                    break
                
                # Update CSV
                reference_count += 1
                self.update_csv_for_new_reference_plane(next_reference_plane, reference_count)
                
                current_reference_plane = next_reference_plane
                
                print(f"Moving to next reference plane: {current_reference_plane['plane_id']}")
            
            # Display final summary
            self.display_all_plane_classification_summary()
            
            return {
                'final_reference_plane': current_reference_plane,
                'total_iterations': iteration + 1,
                'reference_count': reference_count
            }
            
        except Exception as e:
            print(f"ERROR: Iterative reference plane analysis failed: {str(e)}")
            return None
    
    def find_next_reference_plane(self, current_reference_plane):
        """
        Find the next reference plane.
        
        Args:
            current_reference_plane: Current reference plane information
            
        Returns:
            dict: Next reference plane information or None
        """
        try:
            print(f"\nSearching for next reference plane...")
            
            # Sort Unknown planes from CSV by atom count in descending order
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("ERROR: CSV file not found")
                return None
            
            # Collect Unknown planes and atom counts (only planes with 4 or more atoms)
            unknown_planes = {}
            min_atoms_for_reference = 4  # Minimum atom count for reference plane
            
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    plane_id = row.get('plane_id', '')
                    reference_plane = row.get('reference_plane', '')
                    
                    if reference_plane == 'Unknown' and plane_id:
                        unknown_planes[plane_id] = unknown_planes.get(plane_id, 0) + 1
            
            # Filter planes that meet minimum atom count condition
            valid_unknown_planes = {plane_id: count for plane_id, count in unknown_planes.items() 
                                  if count >= min_atoms_for_reference}
            
            if not valid_unknown_planes:
                if unknown_planes:
                    print(f"WARNING: All remaining Unknown planes have fewer than {min_atoms_for_reference} atoms.")
                    print("SUCCESS: Ending reference plane setting task and proceeding to next step.")
                else:
                    print("SUCCESS: All planes have been classified")
                return None
            
            # Select Unknown plane with most atoms
            next_plane_id = max(valid_unknown_planes.items(), key=lambda x: x[1])[0]
            atom_count = valid_unknown_planes[next_plane_id]
            
            print(f"Next reference plane candidate: {next_plane_id} (atom count: {atom_count})")
            
            # Construct plane information
            next_reference_plane = {
                'plane_id': next_plane_id,
                'atom_count': atom_count,
                'plane_type': 'B'  # Second reference plane onwards are B type
            }
            
            return next_reference_plane
            
        except Exception as e:
            print(f"ERROR: Next reference plane search failed: {str(e)}")
            return None
    
    def update_csv_for_new_reference_plane(self, new_reference_plane, reference_count):
        """
        Update CSV for new reference plane.
        
        Args:
            new_reference_plane: New reference plane information
            reference_count: Reference plane number
        """
        try:
            print(f"\nUpdating CSV: Setting reference plane {reference_count}")
            
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("ERROR: CSV file not found")
                return
            
            plane_id = new_reference_plane['plane_id']
            reference_name = f"Reference plane {reference_count}"
            
            # Read CSV file
            rows = []
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                fieldnames = reader.fieldnames
                
                for row in reader:
                    if row.get('plane_id') == plane_id:
                        row['reference_plane'] = reference_name
                        # Set plane_type alphabetically according to reference plane number
                        plane_type_letter = chr(ord('A') + (reference_count - 1) % 26)
                        row['plane_type'] = f"Type {plane_type_letter}"
                    rows.append(row)
            
            # Write CSV file
            with open(csv_path, 'w', encoding='utf-8', newline='') as file:
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)
            
            print(f"SUCCESS: CSV update completed: {plane_id} → {reference_name}")
            
        except Exception as e:
            print(f"ERROR: CSV update failed: {str(e)}")
    
    def display_all_plane_classification_summary(self):
        """
        Display summary of all plane classification results.
        """
        try:
            print("\n" + "=" * 80)
            print("Overall plane classification results summary")
            print("=" * 80)
            
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("ERROR: CSV file not found")
                return
            
            # Collect plane statistics
            plane_stats = {}
            reference_planes = {}
            
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    plane_id = row.get('plane_id', '')
                    reference_plane = row.get('reference_plane', '')
                    plane_type = row.get('plane_type', '')
                    
                    if plane_id:
                        if plane_id not in plane_stats:
                            plane_stats[plane_id] = {
                                'atom_count': 0,
                                'reference_plane': reference_plane,
                                'plane_type': plane_type
                            }
                        plane_stats[plane_id]['atom_count'] += 1
                        
                        # Track reference planes
                        if reference_plane and reference_plane != 'Unknown':
                            reference_planes[reference_plane] = reference_planes.get(reference_plane, 0) + 1
            
            # Display reference plane summary
            print(f"\nReference plane classification:")
            if reference_planes:
                for ref_name, count in sorted(reference_planes.items()):
                    print(f"   {ref_name}: {count} atoms")
            else:
                print("   No reference planes found")
            
            # Display plane type summary
            type_summary = {}
            unknown_count = 0
            
            for plane_id, stats in plane_stats.items():
                plane_type = stats['plane_type']
                atom_count = stats['atom_count']
                
                if stats['reference_plane'] == 'Unknown':
                    unknown_count += atom_count
                else:
                    if plane_type not in type_summary:
                        type_summary[plane_type] = {'planes': 0, 'atoms': 0}
                    type_summary[plane_type]['planes'] += 1
                    type_summary[plane_type]['atoms'] += atom_count
            
            print(f"\nPlane type classification:")
            if type_summary:
                for plane_type, stats in sorted(type_summary.items()):
                    print(f"   {plane_type}: {stats['planes']} planes, {stats['atoms']} atoms")
            
            if unknown_count > 0:
                unknown_planes = sum(1 for stats in plane_stats.values() if stats['reference_plane'] == 'Unknown')
                print(f"   Unknown: {unknown_planes} planes, {unknown_count} atoms")
            
            print(f"\nTotal planes: {len(plane_stats)}")
            print(f"Total atoms: {sum(stats['atom_count'] for stats in plane_stats.values())}")
            
        except Exception as e:
            print(f"ERROR: Failed to display plane classification summary: {str(e)}")
    
    def update_plane_types_in_csv(self, sequence_pattern):
        """
        Update plane types in CSV based on sequence pattern.
        
        Args:
            sequence_pattern (dict): Sequence pattern information
        """
        try:
            print(f"\nUpdating plane types based on sequence pattern...")
            
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("ERROR: CSV file not found")
                return
            
            pattern = sequence_pattern.get('pattern', [])
            if not pattern:
                print("WARNING: No sequence pattern available")
                return
            
            print(f"   Sequence pattern: {' '.join(pattern)}")
            
            # Read and update CSV
            rows = []
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                fieldnames = reader.fieldnames
                
                for row in reader:
                    reference_plane = row.get('reference_plane', '')
                    
                    # Update plane type based on reference plane
                    if reference_plane and reference_plane != 'Unknown':
                        plane_type = self.determine_plane_type_from_reference(reference_plane)
                        if plane_type:
                            row['plane_type'] = plane_type
                    
                    rows.append(row)
            
            # Write updated CSV
            with open(csv_path, 'w', encoding='utf-8', newline='') as file:
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)
            
            print(f"SUCCESS: Plane type update completed")
            
        except Exception as e:
            print(f"ERROR: Plane type update failed: {str(e)}")
    
    def determine_plane_type_from_reference(self, reference_id):
        """
        Determine plane type from reference ID.
        
        Args:
            reference_id (str): Reference plane ID
            
        Returns:
            str: Plane type or None
        """
        if 'Reference plane 1' in reference_id or 'First reference plane' in reference_id:
            return 'Type A'
        elif 'Reference plane 2' in reference_id or 'Second reference plane' in reference_id:
            return 'Type B'
        elif 'Reference plane 3' in reference_id or 'Third reference plane' in reference_id:
            return 'Type C'
        elif 'Reference plane 4' in reference_id or 'Fourth reference plane' in reference_id:
            return 'Type D'
        elif 'Reference plane 5' in reference_id or 'Fifth reference plane' in reference_id:
            return 'Type E'
        else:
            # For additional reference planes, assign types cyclically
            import re
            match = re.search(r'Reference plane (\d+)', reference_id)
            if match:
                ref_num = int(match.group(1))
                type_letter = chr(ord('A') + (ref_num - 1) % 26)
                return f'Type {type_letter}'
        
        return None 