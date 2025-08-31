import re
import numpy as np

class CIFParser:
    """Class for parsing CIF file content"""
    
    @staticmethod
    def parse_cell_parameters(content):
        """
        Extract unit cell parameters from CIF file content.
        
        Args:
            content (str): CIF file content
            
        Returns:
            dict: Unit cell parameters {'a': float, 'b': float, 'c': float,
                                'alpha': float, 'beta': float, 'gamma': float}
        """
        unit_cell = {
            'a': None, 'b': None, 'c': None,
            'alpha': None, 'beta': None, 'gamma': None
        }
        
        if not content:
            return unit_cell
        
        lines = content.split('\n')
        
        for line in lines:
            line = line.strip()
            
            # Extract lattice constants
            if line.startswith('_cell_length_a'):
                match = re.search(r'_cell_length_a\s+([0-9.-]+)', line)
                if match:
                    unit_cell['a'] = float(match.group(1))
            
            elif line.startswith('_cell_length_b'):
                match = re.search(r'_cell_length_b\s+([0-9.-]+)', line)
                if match:
                    unit_cell['b'] = float(match.group(1))
            
            elif line.startswith('_cell_length_c'):
                match = re.search(r'_cell_length_c\s+([0-9.-]+)', line)
                if match:
                    unit_cell['c'] = float(match.group(1))
            
            # Extract lattice angles
            elif line.startswith('_cell_angle_alpha'):
                match = re.search(r'_cell_angle_alpha\s+([0-9.-]+)', line)
                if match:
                    unit_cell['alpha'] = float(match.group(1))
            
            elif line.startswith('_cell_angle_beta'):
                match = re.search(r'_cell_angle_beta\s+([0-9.-]+)', line)
                if match:
                    unit_cell['beta'] = float(match.group(1))
            
            elif line.startswith('_cell_angle_gamma'):
                match = re.search(r'_cell_angle_gamma\s+([0-9.-]+)', line)
                if match:
                    unit_cell['gamma'] = float(match.group(1))
        
        return unit_cell
    
    @staticmethod
    def parse_atomic_positions(content):
        """
        Extract atomic position information from CIF file content.
        
        Args:
            content (str): CIF file content
            
        Returns:
            list: List of atomic information [{'label': str, 'element': str, 
                                  'x': float, 'y': float, 'z': float}, ...]
        """
        atoms = []
        
        if not content:
            return atoms
        
        lines = content.split('\n')
        
        # Find atomic position information in loop_ section
        in_atom_loop = False
        atom_site_headers = {}
        data_started = False
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Check for loop_ start
            if line == 'loop_':
                # Check next lines to determine if this is _atom_site_ section
                j = i + 1
                temp_headers = {}
                while j < len(lines) and lines[j].strip().startswith('_atom_site_'):
                    header_line = lines[j].strip()
                    header_match = re.match(r'_atom_site_(\w+)', header_line)
                    if header_match:
                        header_name = header_match.group(1)
                        temp_headers[header_name] = len(temp_headers)
                    j += 1
                
                # Process this loop_ if _atom_site_ headers are found
                if temp_headers:
                    in_atom_loop = True
                    atom_site_headers = temp_headers
                    data_started = False
                continue
            
            # Process atom loop_ section with headers set
            elif in_atom_loop and atom_site_headers:
                # Skip header lines
                if line.startswith('_atom_site_'):
                    continue
                
                # Skip empty lines or comments
                if not line or line.startswith('#'):
                    continue
                
                # End when new loop_ or other section starts
                if line == 'loop_' or (line.startswith('_') and not line.startswith('_atom_site_')):
                    in_atom_loop = False
                    data_started = False
                    continue
                
                # Process data lines
                columns = line.split()
                if len(columns) >= len(atom_site_headers):
                    data_started = True
                    atom_data = {}
                    
                    # Extract data according to header positions
                    for header, index in atom_site_headers.items():
                        if index < len(columns):
                            value = columns[index]
                            
                            # Remove uncertainty values (values in parentheses)
                            value = re.sub(r'\([^)]*\)', '', value)
                            
                            if header in ['fract_x', 'fract_y', 'fract_z']:
                                try:
                                    atom_data[header] = float(value)
                                except ValueError:
                                    atom_data[header] = 0.0
                            else:
                                atom_data[header] = value
                    
                    # Add atom if all required information is present
                    if all(key in atom_data for key in ['fract_x', 'fract_y', 'fract_z']):
                        atom = {
                            'label': atom_data.get('label', 'Unknown'),
                            'element': atom_data.get('type_symbol', 
                                      atom_data.get('label', 'Unknown')),
                            'x': atom_data['fract_x'],
                            'y': atom_data['fract_y'],
                            'z': atom_data['fract_z']
                        }
                        
                        # Clean element symbol (remove numbers and special characters)
                        atom['element'] = re.sub(r'[0-9+\-_]', '', atom['element'])
                        
                        atoms.append(atom)
        
        return atoms
    
    @staticmethod
    def parse_symmetry_operations(content):
        """
        Extract symmetry operations from CIF file.
        
        Args:
            content (str): CIF file content
            
        Returns:
            list: List of symmetry operation strings
        """
        symmetry_ops = []
        
        if not content:
            return symmetry_ops
        
        lines = content.split('\n')
        in_symmetry_section = False
        
        for line in lines:
            line = line.strip()
            
            # Start of symmetry operation section
            if line.startswith('_symmetry_equiv_pos_as_xyz') or 'symmetry_equiv_pos_site_id' in line:
                in_symmetry_section = True
                continue
            
            # Process symmetry operation data
            elif in_symmetry_section and line and not line.startswith('_') and not line.startswith('#'):
                # Lines starting with numbers (symmetry operations)
                parts = line.split()
                if len(parts) >= 2:
                    # First is number, second and onwards are symmetry operations
                    symmetry_op = ' '.join(parts[1:]).strip("'\"")
                    symmetry_ops.append(symmetry_op)
            
            # End symmetry operation section when new section starts
            elif in_symmetry_section and line.startswith('_') and not 'symmetry' in line:
                break
        
        # Add identity operation if no symmetry operations found
        if not symmetry_ops:
            symmetry_ops.append('x, y, z')
        
        return symmetry_ops
    
    @staticmethod
    def apply_symmetry_operations(base_atoms, symmetry_ops):
        """
        Generate all atoms by applying symmetry operations to base atoms.
        
        Args:
            base_atoms (list): Base atomic positions
            symmetry_ops (list): List of symmetry operations
            
        Returns:
            list: All atoms with symmetry operations applied
        """
        all_atoms = []
        
        for base_atom in base_atoms:
            for sym_op in symmetry_ops:
                new_atom = CIFParser._apply_single_symmetry_operation(base_atom, sym_op)
                if new_atom:
                    # Check for duplicate atoms (remove atoms at same position)
                    is_duplicate = False
                    for existing_atom in all_atoms:
                        if (abs(existing_atom['x'] - new_atom['x']) < 1e-6 and
                            abs(existing_atom['y'] - new_atom['y']) < 1e-6 and
                            abs(existing_atom['z'] - new_atom['z']) < 1e-6 and
                            existing_atom['element'] == new_atom['element']):
                            is_duplicate = True
                            break
                    
                    if not is_duplicate:
                        all_atoms.append(new_atom)
        
        return all_atoms
    
    @staticmethod
    def _apply_single_symmetry_operation(atom, sym_op):
        """
        Apply symmetry operation to a single atom.
        
        Args:
            atom (dict): Atom information
            sym_op (str): Symmetry operation string (e.g., 'x, y, z' or '-x, y+1/2, -z+1/2')
            
        Returns:
            dict: New atom with symmetry operation applied
        """
        try:
            # Fractional coordinates of atom
            x, y, z = atom['x'], atom['y'], atom['z']
            
            # Parse symmetry operation
            operations = [op.strip() for op in sym_op.split(',')]
            if len(operations) != 3:
                return None
            
            # Apply symmetry operation to each coordinate
            new_coords = []
            for op in operations:
                # Variable substitution
                op = op.replace('x', str(x))
                op = op.replace('y', str(y))
                op = op.replace('z', str(z))
                
                # Handle fractional expressions
                op = op.replace('1/2', '0.5')
                op = op.replace('1/3', str(1/3))
                op = op.replace('2/3', str(2/3))
                op = op.replace('1/4', '0.25')
                op = op.replace('3/4', '0.75')
                op = op.replace('1/6', str(1/6))
                op = op.replace('5/6', str(5/6))
                
                # Calculate result
                result = eval(op)
                
                # Normalize fractional coordinates to [0, 1) range
                result = result % 1.0
                
                new_coords.append(result)
            
            # Create new atom
            new_atom = {
                'label': atom['label'],
                'element': atom['element'],
                'x': new_coords[0],
                'y': new_coords[1],
                'z': new_coords[2]
            }
            
            return new_atom
            
        except Exception as e:
            print(f"Warning: Failed to apply symmetry operation '{sym_op}' to atom {atom['label']}: {e}")
            return None 