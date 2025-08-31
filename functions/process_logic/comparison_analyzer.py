"""
Functions related to plane comparison analysis
"""

class ComparisonAnalyzer:
    """Plane comparison analysis class"""

    def __init__(self):
        """Initialize ComparisonAnalyzer"""
        pass
    
    def analyze_failure_reason(self, analysis_result, step1_result=None, step2_result=None, step3_result=None):
        """
        Analyzes plane analysis failure reasons in detail.
        
        Args:
            analysis_result: Atom analysis result
            step1_result: Position distance verification result
            step2_result: Element type verification result
            step3_result: Bond analysis result
            
        Returns:
            dict: Detailed failure reason information
        """
        failure_details = {
            'main_reason': None,
            'sub_reason': None,
            'details': {}
        }

        if not analysis_result:
            failure_details['main_reason'] = 'atom_analysis_failed'
            failure_details['sub_reason'] = 'Atom analysis failed'
            return failure_details

        if step1_result and not step1_result.get('all_passed', False):
            failure_details['main_reason'] = 'position_check_failed'
            failure_details['sub_reason'] = 'Position distance verification failed'
            failure_details['details'] = {
                'max_distance': step1_result.get('max_distance', 0.0),
                'threshold': step1_result.get('threshold', 0.0),
                'failed_atoms': step1_result.get('failed_atoms', [])
            }
            return failure_details

        if step2_result and not step2_result.get('all_passed', False):
            failure_details['main_reason'] = 'element_mismatch'
            failure_details['sub_reason'] = 'Element type mismatch'
            failure_details['details'] = {
                'expected_elements': step2_result.get('expected_elements', []),
                'found_elements': step2_result.get('found_elements', [])
            }
            return failure_details

        if step3_result and not step3_result.get('all_passed', False):
            failure_details['main_reason'] = 'bond_analysis_failed'
            failure_details['sub_reason'] = 'Bond analysis mismatch'
            failure_details['details'] = {
                'bond_differences': step3_result.get('bond_differences', [])
            }
            return failure_details

        return failure_details

    def test_plane_comparison_analysis(self, plane_assignments, reference_plane, comparison_analyzer, all_supercell_atoms, miller_indices, reference_count=1):
        """
        Tests plane comparison analysis (2 steps: center point, atom selection and bond analysis)
        
        Args:
            plane_assignments (dict): Plane assignment results
            reference_plane (dict): Reference plane information
            comparison_analyzer: Plane comparison analyzer
            all_supercell_atoms: All supercell atom information
            miller_indices: Miller indices
        """
        try:
            assignments = plane_assignments['plane_assignments']
            reference_plane_id = reference_plane['plane_id']
            
            print(f"📋 Reference plane #{reference_count}: {reference_plane_id}")
            print(f"📊 Total number of planes: {len(assignments)}")
            
            # Analyze reference plane
            print(f"\n🎯 Analysis of reference plane #{reference_count} ({reference_plane_id}):")
            if reference_plane_id in assignments:
                ref_assignment = assignments[reference_plane_id]
                
                print(f"\n📊 Reference plane basic info: Atom count: {ref_assignment.get('atom_count', 0)} | Plane equation: {ref_assignment.get('plane_info', {}).get('equation_string', 'Unknown')}")
                
                # Reference plane uses all atoms (not just 4 selected)
                print("\n🔍 Using all atoms in reference plane...")
                reference_plane_atoms = ref_assignment.get('atoms', [])
                
                if not reference_plane_atoms:
                    print("❌ No atoms in reference plane.")
                    return
                
                print(f"✅ Reference plane atom count: {len(reference_plane_atoms)}")
                
                # Add atom_id to reference plane atoms (if not present)
                for i, atom in enumerate(reference_plane_atoms):
                    if 'atom_id' not in atom:
                        atom['atom_id'] = atom.get('label', f"{atom.get('element', 'Unknown')}{i+1}")
                
            else:
                print(f"❌ Cannot find reference plane {reference_plane_id} in assignments.")
                return

            # Find "Unknown" planes from CSV file
            from .plane_analyzer import PlaneAnalyzer
            plane_analyzer = PlaneAnalyzer()
            unclassified_planes = plane_analyzer.get_unclassified_planes(reference_plane_id)
            
            # Match with plane information from plane_assignments
            test_planes = []
            for plane_info in unclassified_planes:
                plane_id = plane_info['plane_id']
                if plane_id in assignments:
                    assignment = assignments[plane_id]
                    test_planes.append((plane_id, assignment))
            
            if not test_planes:
                print("No 'Unknown' planes to compare.")
                return
            
            print(f"\n🔍 Analyzing comparison target planes (total {len(test_planes)}):")
            
            # List to store analysis results
            comparison_results = []
            
            for i, (plane_id, assignment) in enumerate(test_planes):
                print(f"\nComparison plane: {plane_id}")
                
                # Initialize result storage for each plane
                plane_result = {
                    'plane_id': plane_id,
                    'atom_count': assignment.get('atom_count', 0),
                    'analysis_success': False,
                    'step1_passed': False,
                    'step2_passed': False,
                    'step3_passed': False,
                    'is_equivalent': False,
                    'error_message': None,
                    'failure_details': None,  # Field to store detailed failure information
                    'verification_result': None,
                    'bond_analysis': None,
                    'cloned_group': None,
                    'moved_group': None
                }
                
                # Step 1: Center point and atom selection
                analysis_result = comparison_analyzer.analyze_comparison_plane_atoms(assignment)
                if not analysis_result:
                    plane_result['failure_details'] = self.analyze_failure_reason(analysis_result)
                    comparison_results.append(plane_result)
                    print(f"{plane_id} atom analysis failed")
                    continue
                
                plane_result['analysis_success'] = True
                comparison_analyzer.display_analysis_summary(analysis_result)
                
                # Step 1.5: Verify plane position of selected atoms
                print(f"\n{plane_id} atom position verification...")
                plane_info = assignment.get('plane_info', {})
                selected_atoms = analysis_result['selected_atoms']
                verification_result = comparison_analyzer.verify_atoms_on_plane(
                    selected_atoms, plane_info
                )
                
                if verification_result:
                    plane_result['step1_passed'] = True
                    plane_result['verification_result'] = verification_result
                else:
                    plane_result['failure_details'] = self.analyze_failure_reason(
                        analysis_result, step1_result=verification_result
                    )
                    comparison_results.append(plane_result)
                    print(f"{plane_id} atom position verification failed")
                    continue
                
                # Step 2: Atom cloning
                print(f"\n{plane_id} atom cloning...")
                cloned_group = comparison_analyzer.clone_selected_atoms(selected_atoms)
                if not cloned_group:
                    plane_result['error_message'] = "Atom cloning failed"
                    comparison_results.append(plane_result)
                    print(f"{plane_id} atom cloning failed")
                    continue
                
                plane_result['cloned_group'] = cloned_group
                
                # Step 3: Move cloned atoms to reference plane
                print(f"\n{plane_id} moving cloned atoms to reference plane...")
                # Convert Miller indices to tuple
                miller_tuple = (miller_indices['h'], miller_indices['k'], miller_indices['l'])
                moved_group = comparison_analyzer.move_cloned_atoms_to_reference_plane(
                    cloned_group,
                    reference_plane,
                    miller_tuple  # Pass in converted format
                )
                
                if not moved_group:
                    plane_result['error_message'] = "Atom movement failed"
                    comparison_results.append(plane_result)
                    print(f"{plane_id} atom movement failed")
                    continue
                
                plane_result['moved_group'] = moved_group
                
                # Step 4: Equivalence test Step 1 - Position distance check
                print(f"{plane_id} equivalence test step 1: Position distance check")
                step1_result = comparison_analyzer.perform_equivalence_test_step1(
                    moved_group,
                    reference_plane_atoms
                )
                
                if step1_result and step1_result.get('all_passed', False):
                    plane_result['step1_passed'] = True
                    print(f"{plane_id} Step 1 passed")
                else:
                    print(f"{plane_id} Step 1 failed")
                    plane_result['failure_details'] = self.analyze_failure_reason(
                        analysis_result, step1_result=step1_result
                    )
                    comparison_results.append(plane_result)
                    continue
                
                # Step 5: Equivalence test Step 2 - Check element types
                print(f"�� {plane_id} 동등성 검사 2단계: 원소 종류 확인")
                step2_result = comparison_analyzer.perform_equivalence_test_step2(
                    moved_group,
                    reference_plane_atoms,
                    step1_result
                )
                
                if step2_result and step2_result.get('all_passed', False):
                    plane_result['step2_passed'] = True
                    print(f"✅ {plane_id} Step 2 통과")
                else:
                    print(f"❌ {plane_id} Step 2 실패")
                    plane_result['failure_details'] = self.analyze_failure_reason(
                        analysis_result, step1_result=step1_result, step2_result=step2_result
                    )
                    comparison_results.append(plane_result)
                    continue
                
                # Step 6: Equivalence test Step 3 - Bonding analysis comparison
                if all_supercell_atoms:
                    print(f"🔍 {plane_id} 동등성 검사 3단계: 결합 분석 비교")
                    step3_result = comparison_analyzer.perform_equivalence_test_step3(
                        moved_group['moved_atoms'],
                        reference_plane_atoms,
                        reference_plane,
                        all_supercell_atoms,
                        step2_result  # Pass atom pair matching information established in Step 2
                    )
                    
                    if step3_result and step3_result.get('all_passed', False):
                        plane_result['step3_passed'] = True
                        print(f"✅ {plane_id} Step 3 통과")
                    else:
                        print(f"❌ {plane_id} Step 3 실패")
                        plane_result['failure_details'] = self.analyze_failure_reason(
                            analysis_result, step1_result=step1_result, 
                            step2_result=step2_result, step3_result=step3_result
                        )
                        comparison_results.append(plane_result)
                        continue
                
                # If all steps pass, determine as equivalent plane
                if plane_result['step1_passed'] and plane_result['step2_passed'] and plane_result['step3_passed']:
                    plane_result['is_equivalent'] = True
                    print(f"🎉 {plane_id}는 기준평면과 동등한 평면입니다!")
                    
                    # Update CSV file
                    import os
                    csv_path = os.path.join('output', 'supercell_coordinates.csv')
                    reference_plane_name = f"제{reference_count}기준평면"
                    
                    # Determine type of plane equivalent to reference plane (same type as reference plane)
                    # Set plane_type in alphabetical order according to reference plane number
                    plane_type_letter = chr(ord('A') + (reference_count - 1) % 26)
                    plane_type = f"Type {plane_type_letter}"
                    
                    success = comparison_analyzer.update_csv_for_equivalent_plane(
                        csv_path, plane_id, reference_plane_name, plane_type
                    )
                    
                    if success:
                        print(f"✅ {plane_id} CSV 업데이트 완료: {reference_plane_name}과 동등, 타입: {plane_type}")
                    else:
                        print(f"⚠️ {plane_id} CSV 업데이트 실패")
                
                comparison_results.append(plane_result)
                print(f"✅ {plane_id} 분석 완료")
            
            # Final result summary
            print("\n📊 전체 분석 결과 요약:")
            print(f"   - 총 비교 평면 수: {len(test_planes)} | 분석 성공: {sum(1 for r in comparison_results if r['analysis_success'])}개 | 위치 검증 통과: {sum(1 for r in comparison_results if r['step1_passed'])}개 | 원소 종류 통과: {sum(1 for r in comparison_results if r['step2_passed'])}개 | 결합 분석 통과: {sum(1 for r in comparison_results if r['step3_passed'])}개 | 동등한 평면 수: {sum(1 for r in comparison_results if r['is_equivalent'])}개")
            
            # Display final analysis result summary
            self.display_comparison_summary(comparison_results, reference_plane_id, reference_count)
            
            return comparison_results
            
        except Exception as e:
            print(f"❌ 평면 비교 분석 실패: {str(e)}")
            return None
    
    def display_comparison_summary(self, comparison_results, reference_plane_id, reference_count=1):
        """
        Displays plane comparison analysis result summary in table format.
        
        Args:
            comparison_results: 비교 분석 결과 리스트
            reference_plane_id: 기준평면 ID
            reference_count: 기준평면 번호 (기본값: 1)
        """
        try:
            if not comparison_results:
                print("❌ 표시할 분석 결과가 없습니다.")
                return

            print("\n" + "=" * 100)
            print("📊 평면 비교 분석 결과 요약")
            print("=" * 100)
            
            total_planes = len(comparison_results)
            successful_analyses = sum(1 for result in comparison_results if result['analysis_success'])
            equivalent_planes = sum(1 for result in comparison_results if result['is_equivalent'])
            
            print(f"🎯 제{reference_count}기준평면: {reference_plane_id}")
            print(f"📋 총 비교 평면 수: {total_planes} | ✅ 분석 성공: {successful_analyses}/{total_planes} ({successful_analyses/total_planes*100:.1f}%) | 🎉 동등 평면: {equivalent_planes}/{total_planes} ({equivalent_planes/total_planes*100:.1f}%)")
            
            # Table header
            print("\n" + "-" * 100)
            print(f"{'평면ID':<10} {'원자수':<8} {'분석':<6} {'Step1':<6} {'Step2':<6} {'Step3':<6} {'동등성':<8} {'실패원인':<30}")
            print("-" * 100)
            
            # Display equivalent planes first
            equivalent_results = [r for r in comparison_results if r['is_equivalent']]
            non_equivalent_results = [r for r in comparison_results if not r['is_equivalent']]
            
            # Equivalent planes
            for result in sorted(equivalent_results, key=lambda x: int(x['plane_id'].replace('Plane ', ''))):
                plane_id = result['plane_id']
                atom_count = result['atom_count']
                analysis = "✅" if result['analysis_success'] else "❌"
                step1 = "✅" if result['step1_passed'] else "❌"
                step2 = "✅" if result['step2_passed'] else "❌"
                step3 = "✅" if result['step3_passed'] else "❌"
                equivalent = "🎉 동등" if result['is_equivalent'] else "❌ 다름"
                failure_reason = "-"
                
                print(f"{plane_id:<10} {atom_count:<8} {analysis:<6} {step1:<6} {step2:<6} {step3:<6} {equivalent:<8} {failure_reason:<30}")
            
            # Separator line
            if equivalent_results and non_equivalent_results:
                print("-" * 100)
            
            # Non-equivalent planes
            for result in sorted(non_equivalent_results, key=lambda x: int(x['plane_id'].replace('Plane ', ''))):
                plane_id = result['plane_id']
                atom_count = result['atom_count']
                analysis = "✅" if result['analysis_success'] else "❌"
                step1 = "✅" if result['step1_passed'] else "❌"
                step2 = "✅" if result['step2_passed'] else "❌"
                step3 = "✅" if result['step3_passed'] else "❌"
                equivalent = "🎉 동등" if result['is_equivalent'] else "❌ 다름"
                
                # Simplify failure reason
                failure_reason = "-"
                if result['failure_details']:
                    main_reason = result['failure_details']['main_reason']
                    if main_reason == 'position_check_failed':
                        failure_reason = "위치거리 검증실패"
                    elif main_reason == 'element_mismatch':
                        failure_reason = "원소종류 불일치"
                    elif main_reason == 'bond_analysis_failed':
                        failure_reason = "결합분석 불일치"
                    elif main_reason == 'atom_analysis_failed':
                        failure_reason = "원자분석 실패"
                    else:
                        failure_reason = "기타 오류"
                elif result['error_message']:
                    failure_reason = result['error_message'][:20] + "..." if len(result['error_message']) > 20 else result['error_message']
                
                print(f"{plane_id:<10} {atom_count:<8} {analysis:<6} {step1:<6} {step2:<6} {step3:<6} {equivalent:<8} {failure_reason:<30}")
            
            print("-" * 100)
            
            # Statistical summary
            step1_passed = sum(1 for r in comparison_results if r['step1_passed'])
            step2_passed = sum(1 for r in comparison_results if r['step2_passed'])
            step3_passed = sum(1 for r in comparison_results if r['step3_passed'])
            
            print(f"📊 단계별 통과율: Step1({step1_passed}/{total_planes}, {step1_passed/total_planes*100:.1f}%) | Step2({step2_passed}/{total_planes}, {step2_passed/total_planes*100:.1f}%) | Step3({step3_passed}/{total_planes}, {step3_passed/total_planes*100:.1f}%)")
            print("=" * 100)
                        
        except Exception as e:
            print(f"⚠️ 요약 정보 표시 실패: {str(e)}") 