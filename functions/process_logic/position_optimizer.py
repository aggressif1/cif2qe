"""
Phase5 position optimization related functions
"""
import numpy as np
import copy

class PositionOptimizer:
    """Phase5 position optimization class"""
    
    def __init__(self):
        """Initialize Phase5PositionOptimizer"""
        self.adjustment_distance = 1.0  # Adjust by 1.0 Å
    
    def optimize_positions(self, structure, collisions):
        """
        Optimize positions of colliding hydrogen atoms.
        
        Args:
            structure (dict): Structure data with added hydrogen
            collisions (list): List of collision information
            
        Returns:
            dict: Structure data with optimized positions
        """
        try:
            optimized_structure = copy.deepcopy(structure)
            atoms = optimized_structure.get('atoms', [])
            
            print(f"   Resolving collisions: {len(collisions)} collisions")
            
            for collision in collisions:
                h_idx = collision['hydrogen_index']
                o_idx = collision['other_index']
                distance = collision['distance']
                
                h_pos = np.array(collision['hydrogen_pos'])
                o_pos = np.array(collision['other_pos'])
                
                # Calculate collision direction vector (from other atom to hydrogen)
                direction = h_pos - o_pos
                direction_norm = np.linalg.norm(direction)
                
                if direction_norm > 0:
                    # Normalize to unit vector
                    direction_unit = direction / direction_norm
                    
                    # Move hydrogen away from other atom by 1Å
                    new_h_pos = o_pos + direction_unit * (distance + self.adjustment_distance)
                    
                    # Update atom position
                    atoms[h_idx]['x'] = new_h_pos[0]
                    atoms[h_idx]['y'] = new_h_pos[1]
                    atoms[h_idx]['z'] = new_h_pos[2]
                    
                    print(f"     H{h_idx} position adjusted: {distance:.3f}Å → {distance + self.adjustment_distance:.3f}Å")
            
            print("   SUCCESS: Position optimization completed")
            return optimized_structure
            
        except Exception as e:
            print(f"ERROR: Error during position optimization: {str(e)}")
            return structure  # Return original on error 