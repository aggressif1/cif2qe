"""
Plane analysis related functions
"""
import csv
import os

class PlaneAnalyzer:
    """Plane analysis class"""

    def __init__(self):
        """Initialize PlaneAnalyzer"""
        pass

    def get_unclassified_planes(self, current_reference_plane_id):
        """
        Find planes with 'Unknown' reference_plane from CSV file.
        Exclude all reference planes and already classified planes.

        Args:
            current_reference_plane_id (str): Current reference plane ID (for verification purposes)

        Returns:
            list: List of Unknown plane information
        """
        try:
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print(f"ERROR: CSV file not found: {csv_path}")
                return []

            # Read CSV file
            plane_counts = {}  # plane_id -> atom_count
            unknown_planes = set()
            all_reference_planes = set()  # Collect all reference planes
            classified_planes = set()  # Collect already classified planes

            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)

                for row in reader:
                    plane_id = row.get('plane_id', '')
                    reference_plane = row.get('reference_plane', '')

                    # Count atoms
                    if plane_id:
                        plane_counts[plane_id] = plane_counts.get(plane_id, 0) + 1

                    # Collect all reference planes (reference plane format)
                    if reference_plane and 'reference_plane' in reference_plane and plane_id:
                        all_reference_planes.add(plane_id)

                    # Collect already classified planes (non-reference planes that are not Unknown)
                    if reference_plane and reference_plane != 'Unknown' and 'reference_plane' not in reference_plane and plane_id:
                        classified_planes.add(plane_id)

                    # Collect Unknown planes
                    if reference_plane == 'Unknown' and plane_id:
                        unknown_planes.add(plane_id)

            # Include only Unknown planes excluding reference planes and already classified planes, with at least 4 atoms
            result_planes = []
            for plane_id in unknown_planes:
                # Exclude reference planes and already classified planes
                if plane_id not in all_reference_planes and plane_id not in classified_planes:
                    atom_count = plane_counts.get(plane_id, 0)
                    if atom_count >= 4:
                        result_planes.append({
                            'plane_id': plane_id,
                            'atom_count': atom_count
                        })

            # Sort by plane ID
            result_planes.sort(key=lambda x: x['plane_id'])

            print(f"Number of analyzable Unknown planes: {len(result_planes)}")
            print(f"Excluded reference planes: {len(all_reference_planes)}")
            print(f"Excluded classified planes: {len(classified_planes)}")

            for plane in result_planes[:5]:  # Preview first 5 only
                print(f"   - {plane['plane_id']}: {plane['atom_count']} atoms")
            if len(result_planes) > 5:
                print(f"   ... and {len(result_planes) - 5} more planes")

            return result_planes

        except Exception as e:
            print(f"ERROR: Unknown plane search failed: {str(e)}")
            return []

    def extract_d_value_from_equation(self, plane_equation):
        """
        Extract D value from plane equation.

        Args:
            plane_equation (str): Plane equation string

        Returns:
            float: D value or None
        """
        try:
            import re

            # Extract d from "ax + by + cz + d = 0" or "ax + by + cz - d = 0" format
            # Find constant term before equals sign
            if '=' in plane_equation:
                left_side = plane_equation.split('=')[0].strip()

                # Find constant term without variables (usually the last term)
                # Use regex to find all terms and extract last constant term
                pattern = r'([+-]?\s*\d+\.?\d*)\s*(?![xyz])'
                matches = re.findall(pattern, left_side)

                if matches:
                    # Use last constant term as D value
                    d_str = matches[-1].replace(' ', '')
                    d_value = float(d_str)
                    return -d_value  # Move d to right side: ax + by + cz + d = 0 becomes -d

            return None

        except Exception as e:
            print(f"WARNING: D value extraction failed: {str(e)}")
            return None