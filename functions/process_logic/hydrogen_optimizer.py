#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase5 hydrogen position optimizer
Responsible for hydrogen position optimization and collision detection.
"""

class HydrogenOptimizer:
    """Phase5 hydrogen position optimizer"""
    
    def __init__(self):
        """Initialize Phase5HydrogenOptimizer"""
        from .collision_detector import CollisionDetector
        from .position_optimizer import PositionOptimizer

        self.collision_detector = CollisionDetector()
        self.position_optimizer = PositionOptimizer()
    
    def validate_input_data(self, phase4_data):
        """Validate input data"""
        if not phase4_data:
            print("ERROR: No Phase4 data available.")
            return False
        
        # Check structure type
        structure_type = phase4_data.get('structure_type')
        if structure_type != 'slab':
            print("ERROR: Phase5 can only be executed with slab structures.")
            print(f"   Current structure type: {structure_type}")
            return False
        
        # Check slab structure data
        if not phase4_data.get('slab_structure'):
            print("ERROR: No slab structure data available.")
            return False
        
        return True
    
    def optimize_hydrogen_positions(self, hydrogenated_structure):
        """Optimize hydrogen positions"""
        try:
            # Collision detection
            print("   Checking hydrogen-atom collisions...")
            collisions = self.collision_detector.detect_collisions(hydrogenated_structure)
            
            if collisions:
                print(f"   WARNING: {len(collisions)} collisions detected")
                
                # Position optimization
                print("   Optimizing positions to resolve collisions...")
                optimized_structure = self.position_optimizer.optimize_positions(
                    hydrogenated_structure, collisions
                )
                
                # Re-check
                remaining_collisions = self.collision_detector.detect_collisions(optimized_structure)
                if remaining_collisions:
                    print(f"   WARNING: {len(remaining_collisions)} collisions still remain")
                else:
                    print("   SUCCESS: All collisions resolved")
                
                return optimized_structure
            else:
                print("   SUCCESS: No collisions - optimization unnecessary")
                return hydrogenated_structure
                
        except Exception as e:
            print(f"ERROR: Error during hydrogen position optimization: {str(e)}")
            return hydrogenated_structure  # Return original on error 