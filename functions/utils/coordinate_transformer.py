#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coordinate transformation utility
Provides mutual conversion between fractional coordinates and absolute coordinates.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Union

class CoordinateTransformer:
    """Coordinate transformation class"""
    
    def __init__(self):
        """Initialize CoordinateTransformer"""
        pass
    
    def create_cell_matrix(self, a_vector: List[float], b_vector: List[float], c_vector: List[float]) -> np.ndarray:
        """
        Create transformation matrix from cell vectors.
        
        Args:
            a_vector (List[float]): a cell vector [x, y, z]
            b_vector (List[float]): b cell vector [x, y, z]
            c_vector (List[float]): c cell vector [x, y, z]
            
        Returns:
            np.ndarray: 3x3 cell matrix (composed of column vectors)
        """
        try:
            # Create cell matrix: 3x3 matrix with each cell vector as columns
            cell_matrix = np.column_stack([
                np.array(a_vector),
                np.array(b_vector),
                np.array(c_vector)
            ])
            
            return cell_matrix
            
        except Exception as e:
            print(f"ERROR: Error during cell matrix creation: {str(e)}")
            return None
    
    def fractional_to_absolute_single(self, frac_coords: List[float], cell_matrix: np.ndarray) -> Optional[np.ndarray]:
        """
        Convert fractional coordinates of a single atom to absolute coordinates.
        
        Args:
            frac_coords (List[float]): fractional coordinates [u, v, w]
            cell_matrix (np.ndarray): cell transformation matrix
            
        Returns:
            Optional[np.ndarray]: absolute coordinates [x, y, z] or None (if failed)
        """
        try:
            frac_pos = np.array(frac_coords)
            abs_pos = cell_matrix @ frac_pos
            return abs_pos
            
        except Exception as e:
            print(f"ERROR: Error during single coordinate transformation: {str(e)}")
            return None
    
    def absolute_to_fractional_single(self, abs_coords: List[float], cell_matrix: np.ndarray) -> Optional[np.ndarray]:
        """
        Convert absolute coordinates of a single atom to fractional coordinates.
        
        Args:
            abs_coords (List[float]): absolute coordinates [x, y, z]
            cell_matrix (np.ndarray): cell transformation matrix
            
        Returns:
            Optional[np.ndarray]: fractional coordinates [u, v, w] or None (if failed)
        """
        try:
            abs_pos = np.array(abs_coords)
            # Transformation using inverse matrix: fractional coordinates = cell_matrix^(-1) × absolute coordinates
            inv_cell_matrix = np.linalg.inv(cell_matrix)
            frac_pos = inv_cell_matrix @ abs_pos
            return frac_pos
            
        except Exception as e:
            print(f"ERROR: Error during single coordinate inverse transformation: {str(e)}")
            return None
    
    def convert_fractional_to_absolute_coordinates(self, atoms: List[Dict], a_vector: List[float], 
                                                 b_vector: List[float], c_vector: List[float],
                                                 frac_coords_key: str = 'frac_coords_extended',
                                                 abs_coords_keys: Tuple[str, str, str] = ('abs_x', 'abs_y', 'abs_z'),
                                                 verbose: bool = True) -> Optional[List[Dict]]:
        """
        Convert fractional coordinates of multiple atoms to absolute coordinates.
        
        Args:
            atoms (List[Dict]): list of atoms
            a_vector (List[float]): a cell vector
            b_vector (List[float]): b cell vector
            c_vector (List[float]): c cell vector
            frac_coords_key (str): fractional coordinates key name
            abs_coords_keys (Tuple[str, str, str]): absolute coordinates key names (x, y, z)
            verbose (bool): whether to show detailed output
            
        Returns:
            Optional[List[Dict]]: list of atoms with added absolute coordinates or None (if failed)
        """
        try:
            if verbose:
                print(f"   Number of atoms to convert: {len(atoms)}")
            
            # Create cell matrix
            cell_matrix = self.create_cell_matrix(a_vector, b_vector, c_vector)
            if cell_matrix is None:
                return None
            
            if verbose:
                print(f"   Used cell matrix:")
                print(f"     a_vector: [{a_vector[0]:.6f}, {a_vector[1]:.6f}, {a_vector[2]:.6f}]")
                print(f"     b_vector: [{b_vector[0]:.6f}, {b_vector[1]:.6f}, {b_vector[2]:.6f}]")
                print(f"     c_vector: [{c_vector[0]:.6f}, {c_vector[1]:.6f}, {c_vector[2]:.6f}]")
            
            absolute_atoms = []
            conversion_count = 0
            
            for atom in atoms:
                # Extract fractional coordinates
                if frac_coords_key in atom:
                    frac_coords = atom[frac_coords_key]
                else:
                    if verbose:
                        print(f"   WARNING: Atom {atom.get('element', 'X')} does not have {frac_coords_key} information.")
                    continue
                
                # Convert fractional coordinates to absolute coordinates
                abs_pos = self.fractional_to_absolute_single(frac_coords, cell_matrix)
                if abs_pos is None:
                    continue
                
                # Create new atom dictionary (existing info + absolute coordinates)
                new_atom = atom.copy()
                new_atom[abs_coords_keys[0]] = abs_pos[0]  # abs_x
                new_atom[abs_coords_keys[1]] = abs_pos[1]  # abs_y
                new_atom[abs_coords_keys[2]] = abs_pos[2]  # abs_z
                new_atom['original_frac_coords'] = frac_coords  # preserve original fractional coordinates
                
                absolute_atoms.append(new_atom)
                conversion_count += 1
                
                # Output conversion results for first 3 atoms
                if verbose and conversion_count <= 3:
                    element = atom.get('element', 'X')
                    print(f"   Atom {element}: fractional coordinates ({frac_coords[0]:.6f}, {frac_coords[1]:.6f}, {frac_coords[2]:.6f}) → absolute coordinates ({abs_pos[0]:.6f}, {abs_pos[1]:.6f}, {abs_pos[2]:.6f})")
            
            if verbose:
                print(f"   Conversion completed: {conversion_count}/{len(atoms)} atoms")
            
            if conversion_count == 0:
                if verbose:
                    print("   ERROR: No atoms were converted. Please check fractional coordinate information.")
                return None
            
            return absolute_atoms
            
        except Exception as e:
            print(f"   ERROR: Error during fractional → absolute coordinate conversion: {str(e)}")
            return None
    
    def convert_absolute_to_fractional_coordinates(self, atoms: List[Dict], a_vector: List[float], 
                                                 b_vector: List[float], c_vector: List[float],
                                                 abs_coords_keys: Tuple[str, str, str] = ('x', 'y', 'z'),
                                                 frac_coords_key: str = 'frac_coords',
                                                 verbose: bool = True) -> Optional[List[Dict]]:
        """
        Convert absolute coordinates of multiple atoms to fractional coordinates.
        
        Args:
            atoms (List[Dict]): list of atoms
            a_vector (List[float]): a cell vector
            b_vector (List[float]): b cell vector
            c_vector (List[float]): c cell vector
            abs_coords_keys (Tuple[str, str, str]): absolute coordinates key names (x, y, z)
            frac_coords_key (str): fractional coordinates key name
            verbose (bool): whether to show detailed output
            
        Returns:
            Optional[List[Dict]]: list of atoms with added fractional coordinates or None (if failed)
        """
        try:
            if verbose:
                print(f"   Number of atoms to convert: {len(atoms)}")
            
            # Create cell matrix
            cell_matrix = self.create_cell_matrix(a_vector, b_vector, c_vector)
            if cell_matrix is None:
                return None
            
            if verbose:
                print(f"   Used cell matrix:")
                print(f"     a_vector: [{a_vector[0]:.6f}, {a_vector[1]:.6f}, {a_vector[2]:.6f}]")
                print(f"     b_vector: [{b_vector[0]:.6f}, {b_vector[1]:.6f}, {b_vector[2]:.6f}]")
                print(f"     c_vector: [{c_vector[0]:.6f}, {c_vector[1]:.6f}, {c_vector[2]:.6f}]")
            
            fractional_atoms = []
            conversion_count = 0
            
            for atom in atoms:
                # Extract absolute coordinates
                if all(key in atom for key in abs_coords_keys):
                    abs_coords = [atom[abs_coords_keys[0]], atom[abs_coords_keys[1]], atom[abs_coords_keys[2]]]
                else:
                    if verbose:
                        print(f"   WARNING: Atom {atom.get('element', 'X')} does not have absolute coordinate information.")
                    continue
                
                # Convert absolute coordinates to fractional coordinates
                frac_pos = self.absolute_to_fractional_single(abs_coords, cell_matrix)
                if frac_pos is None:
                    continue
                
                # Create new atom dictionary (existing info + fractional coordinates)
                new_atom = atom.copy()
                new_atom[frac_coords_key] = frac_pos.tolist()
                new_atom['original_abs_coords'] = abs_coords  # preserve original absolute coordinates
                
                fractional_atoms.append(new_atom)
                conversion_count += 1
                
                # Output conversion results for first 3 atoms
                if verbose and conversion_count <= 3:
                    element = atom.get('element', 'X')
                    print(f"   Atom {element}: absolute coordinates ({abs_coords[0]:.6f}, {abs_coords[1]:.6f}, {abs_coords[2]:.6f}) → fractional coordinates ({frac_pos[0]:.6f}, {frac_pos[1]:.6f}, {frac_pos[2]:.6f})")
            
            if verbose:
                print(f"   Conversion completed: {conversion_count}/{len(atoms)} atoms")
            
            if conversion_count == 0:
                if verbose:
                    print("   ERROR: No atoms were converted. Please check absolute coordinate information.")
                return None
            
            return fractional_atoms
            
        except Exception as e:
            print(f"   ERROR: Error during absolute → fractional coordinate conversion: {str(e)}")
            return None
    
    def get_cell_parameters(self, a_vector: List[float], b_vector: List[float], c_vector: List[float]) -> Dict[str, float]:
        """
        Calculate lattice parameters from cell vectors.
        
        Args:
            a_vector (List[float]): a cell vector
            b_vector (List[float]): b cell vector
            c_vector (List[float]): c cell vector
            
        Returns:
            Dict[str, float]: lattice parameters dictionary {a, b, c, alpha, beta, gamma}
        """
        try:
            a_vec = np.array(a_vector)
            b_vec = np.array(b_vector)
            c_vec = np.array(c_vector)
            
            # Lattice parameter lengths
            a = np.linalg.norm(a_vec)
            b = np.linalg.norm(b_vec)
            c = np.linalg.norm(c_vec)
            
            # Lattice parameter angles (degrees)
            alpha = np.degrees(np.arccos(np.dot(b_vec, c_vec) / (b * c)))
            beta = np.degrees(np.arccos(np.dot(a_vec, c_vec) / (a * c)))
            gamma = np.degrees(np.arccos(np.dot(a_vec, b_vec) / (a * b)))
            
            return {
                'a': a,
                'b': b,
                'c': c,
                'alpha': alpha,
                'beta': beta,
                'gamma': gamma
            }
            
        except Exception as e:
            print(f"ERROR: Error during lattice parameter calculation: {str(e)}")
            return {} 