"""
Phase5 collision detection related functions
"""
import numpy as np

class CollisionDetector:
    """Phase5 collision detection class"""
    
    def __init__(self):
        """Initialize Phase5CollisionDetector"""
        self.collision_threshold = 2.0  # Consider collision if distance is below 2.0 Å
    
    def detect_collisions(self, structure):
        """
        Detect collisions between hydrogen atoms and other atoms.
        
        Args:
            structure (dict): Structure data with added hydrogen
            
        Returns:
            list: List of collision information
        """
        try:
            atoms = structure.get('atoms', [])
            if not atoms:
                return []
            
            # Separate hydrogen atoms and other atoms
            hydrogen_atoms = []
            other_atoms = []
            
            for i, atom in enumerate(atoms):
                if atom.get('element') == 'H':
                    hydrogen_atoms.append((i, atom))
                else:
                    other_atoms.append((i, atom))
            
            collisions = []
            
            # Check collisions for each hydrogen atom
            for h_idx, h_atom in hydrogen_atoms:
                h_pos = np.array([h_atom.get('x', 0), h_atom.get('y', 0), h_atom.get('z', 0)])
                
                for o_idx, o_atom in other_atoms:
                    o_pos = np.array([o_atom.get('x', 0), o_atom.get('y', 0), o_atom.get('z', 0)])
                    
                    distance = np.linalg.norm(h_pos - o_pos)
                    
                    if distance <= self.collision_threshold:
                        collision_info = {
                            'hydrogen_index': h_idx,
                            'other_index': o_idx,
                            'hydrogen_atom': h_atom,
                            'other_atom': o_atom,
                            'distance': distance,
                            'hydrogen_pos': h_pos.tolist(),
                            'other_pos': o_pos.tolist()
                        }
                        collisions.append(collision_info)
            
            return collisions
            
        except Exception as e:
            print(f"ERROR: Error during collision detection: {str(e)}")
            return [] 