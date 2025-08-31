"""
Sequence analysis related functions separated from Phase3
"""
import csv
import os

class SequenceAnalyzer:
    """Sequence analysis class"""
    
    def __init__(self):
        """Initialize SequenceAnalyzer"""
        pass
    
    def analyze_d_value_sequence(self):
        """
        Analyze repeating sequences by sorting plane_type according to D values.
        
        Returns:
            dict: Sequence analysis result information
        """
        print("\n🔍 Step 3: D-value based sequence analysis")
        print("-" * 50)
        
        try:
            # Read plane information from CSV file
            csv_path = os.path.join('output', 'supercell_coordinates.csv')
            if not os.path.exists(csv_path):
                print("❌ CSV file not found")
                return None
            
            # Collect information by plane
            planes_info = {}
            
            with open(csv_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    plane_id = row.get('plane_id', '')
                    plane_equation = row.get('plane_equation', '')
                    plane_type = row.get('plane_type', '')
                    
                    if plane_id and plane_equation and plane_type and plane_type != 'Unknown':
                        # Extract D value
                        from .plane_analyzer import PlaneAnalyzer
                        analyzer = PlaneAnalyzer()
                        d_value = analyzer.extract_d_value_from_equation(plane_equation)
                        
                        if d_value is not None:
                            if plane_id not in planes_info:
                                planes_info[plane_id] = {
                                    'plane_type': plane_type,
                                    'd_value': d_value,
                                    'plane_equation': plane_equation,
                                    'atom_count': 0
                                }
                            planes_info[plane_id]['atom_count'] += 1
            
            if not planes_info:
                print("❌ No plane information to analyze")
                return None
            
            print(f"📊 Number of planes to analyze: {len(planes_info)}")
            
            # Sort by D value
            sorted_planes = sorted(planes_info.items(), key=lambda x: x[1]['d_value'])
            
            # Display all planes in table format
            self._display_sorted_planes_table(sorted_planes)
            
            # Find repeating sequence
            sequence_pattern = self.find_repeating_sequence(sorted_planes)
            
            if sequence_pattern:
                return {
                    'sorted_planes': sorted_planes,
                    'sequence_pattern': sequence_pattern,
                    'total_planes': len(planes_info)
                }
            else:
                print("❌ Unable to find repeating sequence")
                return None
                
        except Exception as e:
            print(f"❌ Sequence analysis failed: {str(e)}")
            return None
    
    def find_repeating_sequence(self, sorted_planes):
        """
        Find repeating sequences in sorted planes.
        
        Args:
            sorted_planes (list): List of plane information sorted by D values
            
        Returns:
            dict: Repeating sequence information or None
        """
        if len(sorted_planes) < 4:
            return None
        
        # Extract plane_type sequence
        type_sequence = [info['plane_type'] for _, info in sorted_planes]
        
        # Efficiently analyze and output plane_type sequence
        print(f"\n🔍 plane_type sequence analysis:")
        print(f"   Total number of planes: {len(type_sequence)}")
        
        # Analyze unique types and frequencies in sequence
        unique_types = list(set(type_sequence))
        type_counts = {t: type_sequence.count(t) for t in unique_types}
        
        print(f"   Number of unique types: {len(unique_types)}")
        print(f"   Frequency by type: {', '.join([f'{t}({count})' for t, count in sorted(type_counts.items())])}")
        
        # Display entire sequence concisely (considering repetition patterns)
        if len(type_sequence) <= 10:
            # Short sequences display entirely
            print(f"   Complete sequence: {' → '.join(type_sequence)}")
        else:
            # Long sequences display start/end parts and pattern hints
            start_part = ' → '.join(type_sequence[:5])
            end_part = ' → '.join(type_sequence[-3:])
            print(f"   Sequence start: {start_part} → ...")
            print(f"   Sequence end: ... → {end_part}")
            
            # Provide simple pattern hints
            if len(unique_types) == 2:
                # Check if two types appear alternately
                alternating = all(type_sequence[i] != type_sequence[i+1] for i in range(len(type_sequence)-1))
                if alternating:
                    print(f"   Pattern hint: {unique_types[0]} and {unique_types[1]} appear alternately")
        
        # Detect repeating patterns
        pattern_info = self.detect_repeating_pattern(type_sequence)
        
        return pattern_info
    
    def detect_repeating_pattern(self, type_sequence):
        """
        Detects repeating patterns in type sequence.
        
        Args:
            type_sequence (list): plane_type sequence
            
        Returns:
            dict: Pattern information or default pattern
        """
        sequence_length = len(type_sequence)
        
        # Analyze excluding Unclassified types
        classified_sequence = [t for t in type_sequence if t != 'Unclassified']
        classified_length = len(classified_sequence)
        
        if classified_length < 6:  # Need at least 6 for pattern analysis
            print(f"⚠️ Too few classified planes ({classified_length})")
            return self._create_default_pattern(type_sequence)
        
        print(f"🔍 Pattern analysis using only classified planes: {classified_length}")
        print(f"   Classified sequence: {' → '.join(classified_sequence[:10])}{'...' if classified_length > 10 else ''}")
        
        # Step 1: Find exact match pattern
        exact_pattern = self._find_exact_pattern(classified_sequence, sequence_length)
        if exact_pattern:
            return exact_pattern
        
        # Step 2: Find flexible pattern (allow some missing)
        flexible_pattern = self._find_flexible_pattern(classified_sequence, sequence_length)
        if flexible_pattern:
            return flexible_pattern
        
        # If pattern not found, create default pattern
        return self._create_default_pattern(type_sequence)
    
    def _create_default_pattern(self, type_sequence):
        """
        Creates default pattern.
        
        Args:
            type_sequence (list): plane_type sequence
            
        Returns:
            dict: Default pattern information or None
        """
        # Check if 'Type A' and 'Type B' types exist
        has_type_a = 'Type A' in type_sequence
        has_type_b = 'Type B' in type_sequence
        
        if has_type_a and has_type_b:
            # Default pattern: composed of Type A and Type B types
            default_pattern = ['Type A', 'Type B']
            pattern_str = " → ".join(default_pattern)
            print(f"\n⚠️ Could not find clear repeating pattern.")
            print(f"   📋 Using default pattern: [{pattern_str}]")
        elif has_type_a:
            # If only Type A exists
            default_pattern = ['Type A']
            print(f"\n⚠️ Could not find clear repeating pattern.")
            print(f"   📋 Using default pattern: [Type A]")
        elif has_type_b:
            # If only Type B exists
            default_pattern = ['Type B']
            print(f"\n⚠️ Could not find clear repeating pattern.")
            print(f"   📋 Using default pattern: [Type B]")
        else:
            # If neither exists
            print("\n❌ Could not find valid repeating pattern")
            return None
        
        # Return default pattern information
        return {
            'pattern': default_pattern,
            'pattern_length': len(default_pattern),
            'repetitions': 1,  # No actual repetition but set to 1 for processing
            'coverage': 0.5,   # Arbitrary coverage value
            'is_default': True # Indicates this is a default pattern
        }
    
    def _display_sorted_planes_table(self, sorted_planes):
        """
        Outputs planes sorted by D-value in table format.
        
        Args:
            sorted_planes (list): List of plane information sorted by D-value
        """
        print(f"\n📋 Planes sorted by D-value (total {len(sorted_planes)}):")
        print("=" * 135)
        print(f"{'No.':<4} {'Plane ID':<10} {'Plane Equation':<50} {'Plane Type':<15} {'Atoms':<8} {'D-value':<10}")
        print("-" * 135)
        
        for idx, (plane_id, info) in enumerate(sorted_planes, 1):
            plane_equation = info['plane_equation']
            plane_type = info['plane_type']
            d_value = info['d_value']
            atom_count = info['atom_count']
            
            # Truncate if plane equation is too long
            if len(plane_equation) > 47:
                plane_equation = plane_equation[:44] + "..."
            
            print(f"{idx:<4} {plane_id:<10} {plane_equation:<50} {plane_type:<15} {atom_count:<8} {d_value:<10.4f}")
        
        print("=" * 135)
    
    def _find_exact_pattern(self, classified_sequence, total_sequence_length):
        """
        Finds exact match pattern.
        
        Args:
            classified_sequence (list): 분류된 평면 시퀀스
            total_sequence_length (int): 전체 시퀀스 길이
            
        Returns:
            dict: 패턴 정보 또는 None
        """
        classified_length = len(classified_sequence)
        max_pattern_length = min(20, classified_length // 2)
        
        for pattern_length in range(2, max_pattern_length + 1):
            # Extract first pattern
            pattern = classified_sequence[:pattern_length]
            
            # Check if pattern repeats
            repetitions = 0
            valid_pattern = True
            
            for start_idx in range(0, classified_length - pattern_length + 1, pattern_length):
                current_segment = classified_sequence[start_idx:start_idx + pattern_length]
                
                if current_segment == pattern:
                    repetitions += 1
                elif len(current_segment) == pattern_length:
                    # Fail if complete length segment differs from pattern
                    valid_pattern = False
                    break
                # Ignore last incomplete segment
            
            # Need at least 3 repetitions to be considered valid pattern
            if valid_pattern and repetitions >= 3:
                coverage = (repetitions * pattern_length) / classified_length
                pattern_str = " → ".join(pattern)
                
                print(f"\n✅ 완전 일치 패턴 발견!")
                print(f"   📋 패턴: [{pattern_str}]")
                print(f"   📏 패턴 길이: {pattern_length}")
                print(f"   🔄 반복 횟수: {repetitions}회")
                print(f"   📊 분류된 평면 커버리지: {coverage:.1%}")
                print(f"   📊 전체 평면 커버리지: {(repetitions * pattern_length) / total_sequence_length:.1%}")
                
                return {
                    'pattern': pattern,
                    'pattern_length': pattern_length,
                    'repetitions': repetitions,
                    'coverage': coverage,
                    'classified_coverage': coverage,
                    'total_coverage': (repetitions * pattern_length) / total_sequence_length,
                    'pattern_type': 'exact'
                }
        
        return None
    
    def _find_flexible_pattern(self, classified_sequence, total_sequence_length):
        """
        Finds flexible pattern (allows some missing or order changes).
        
        Args:
            classified_sequence (list): 분류된 평면 시퀀스
            total_sequence_length (int): 전체 시퀀스 길이
            
        Returns:
            dict: 패턴 정보 또는 None
        """
        print(f"\n🔍 유연한 패턴 탐지 시작...")
        
        classified_length = len(classified_sequence)
        
        # Analyze frequency of unique types
        from collections import Counter
        type_counts = Counter(classified_sequence)
        unique_types = list(type_counts.keys())
        
        print(f"   고유 타입 수: {len(unique_types)}")
        print(f"   타입별 빈도: {dict(type_counts)}")
        
        # Generate pattern based on first occurrence order of unique types in sequence
        if len(unique_types) >= 8:  # Need at least 8 types for meaningful pattern
            # Find first occurrence position of each type in sequence and sort by order
            first_occurrence = {}
            for i, type_name in enumerate(classified_sequence):
                if type_name not in first_occurrence:
                    first_occurrence[type_name] = i
            
            # Sort types by first occurrence order
            sequence_ordered_types = sorted(unique_types, key=lambda x: first_occurrence[x])
            
            print(f"   시퀀스 순서: {' → '.join(sequence_ordered_types)}")
            
            # Try various pattern lengths (8~18)
            for pattern_length in range(8, min(19, len(sequence_ordered_types) + 1)):
                candidate_pattern = sequence_ordered_types[:pattern_length]
                
                # Check how well this pattern appears in sequence
                pattern_score = self._calculate_pattern_score(classified_sequence, candidate_pattern)
                
                if pattern_score['coverage'] >= 0.6:  # At least 60% coverage
                    pattern_str = " → ".join(candidate_pattern)
                    
                    print(f"\n✅ 유연한 패턴 발견!")
                    print(f"   📋 패턴: [{pattern_str}]")
                    print(f"   📏 패턴 길이: {pattern_length}")
                    print(f"   🔄 예상 반복 횟수: {pattern_score['estimated_repetitions']:.1f}회")
                    print(f"   📊 패턴 커버리지: {pattern_score['coverage']:.1%}")
                    print(f"   📊 전체 평면 커버리지: {pattern_score['coverage'] * classified_length / total_sequence_length:.1%}")
                    
                    return {
                        'pattern': candidate_pattern,
                        'pattern_length': pattern_length,
                        'repetitions': pattern_score['estimated_repetitions'],
                        'coverage': pattern_score['coverage'],
                        'classified_coverage': pattern_score['coverage'],
                        'total_coverage': pattern_score['coverage'] * classified_length / total_sequence_length,
                        'pattern_type': 'flexible'
                    }
        
        return None
    
    def _calculate_pattern_score(self, sequence, pattern):
        """
        Calculates score for how well a given pattern appears in sequence.
        
        Args:
            sequence (list): 분석할 시퀀스
            pattern (list): 패턴 후보
            
        Returns:
            dict: 패턴 점수 정보
        """
        pattern_set = set(pattern)
        
        # Ratio of types included in pattern in sequence
        pattern_elements_in_sequence = [elem for elem in sequence if elem in pattern_set]
        coverage = len(pattern_elements_in_sequence) / len(sequence) if sequence else 0
        
        # Expected repetition count (relative to pattern length)
        estimated_repetitions = len(pattern_elements_in_sequence) / len(pattern) if pattern else 0
        
        return {
            'coverage': coverage,
            'estimated_repetitions': estimated_repetitions,
            'pattern_elements_count': len(pattern_elements_in_sequence)
        }