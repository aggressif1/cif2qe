"""
Temporary file storage handler
Creating and managing temporary files for data transfer between Phases
"""

import os
import json
from typing import Dict, Any
from datetime import datetime


class TemporaryFileHandler:
    """Class responsible for temporary file storage and management"""
    
    def __init__(self, base_dir: str = None):
        """
        Initialize temporary file handler
        
        Args:
            base_dir (str): base directory (if None, use project root/modules/temporary)
        """
        if base_dir is None:
            # Find project root
            current_dir = os.path.dirname(os.path.abspath(__file__))
            project_root = os.path.dirname(os.path.dirname(current_dir))
            self.temp_dir = os.path.join(project_root, "modules", "temporary")
        else:
            self.temp_dir = base_dir
        
        # Create directory
        os.makedirs(self.temp_dir, exist_ok=True)
    
    def save_phase_results(self, phase_name: str, data: Dict[str, Any]) -> bool:
        """
        Save Phase results to temporary file.
        
        Args:
            phase_name (str): Phase name (e.g., "phase1", "phase2")
            data (Dict[str, Any]): data to save
            
        Returns:
            bool: save success status
        """
        try:
            # Save entire data as JSON file
            json_path = os.path.join(self.temp_dir, f"{phase_name}_results.json")
            json_data = self._prepare_json_data(data)
            
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(json_data, f, indent=2, ensure_ascii=False)
            
            # Save text summary file
            summary_path = os.path.join(self.temp_dir, f"{phase_name}_summary.txt")
            summary_content = self._generate_summary(phase_name, data)
            
            with open(summary_path, 'w', encoding='utf-8') as f:
                f.write(summary_content)
            
            print(f"{phase_name.upper()} temporary file save completed:")
            print(f"   - JSON: {json_path}")
            print(f"   - Summary: {summary_path}")
            
            return True
            
        except Exception as e:
            print(f"WARNING: Error saving {phase_name.upper()} temporary file: {str(e)}")
            return False
    
    def load_phase_results(self, phase_name: str) -> Dict[str, Any]:
        """
        Load Phase results from temporary file.
        
        Args:
            phase_name (str): Phase name (e.g., "phase1", "phase2")
            
        Returns:
            Dict[str, Any]: loaded data, empty dictionary if failed
        """
        try:
            json_path = os.path.join(self.temp_dir, f"{phase_name}_results.json")
            
            if not os.path.exists(json_path):
                print(f"WARNING: {phase_name.upper()} result file not found: {json_path}")
                return {}
            
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            print(f"{phase_name.upper()} result load completed: {json_path}")
            return data
            
        except Exception as e:
            print(f"WARNING: Error loading {phase_name.upper()} result: {str(e)}")
            return {}
    
    def cleanup_phase_files(self, phase_name: str) -> bool:
        """
        Delete temporary files for specific Phase.
        
        Args:
            phase_name (str): Phase name
            
        Returns:
            bool: deletion success status
        """
        try:
            json_path = os.path.join(self.temp_dir, f"{phase_name}_results.json")
            summary_path = os.path.join(self.temp_dir, f"{phase_name}_summary.txt")
            
            removed_files = []
            for file_path in [json_path, summary_path]:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    removed_files.append(os.path.basename(file_path))
            
            if removed_files:
                print(f"{phase_name.upper()} temporary file deletion completed: {', '.join(removed_files)}")
            
            return True
            
        except Exception as e:
            print(f"WARNING: Error deleting {phase_name.upper()} temporary files: {str(e)}")
            return False
    
    def cleanup_all_files(self) -> bool:
        """
        Delete all temporary files.
        
        Returns:
            bool: deletion success status
        """
        try:
            if not os.path.exists(self.temp_dir):
                return True
            
            removed_files = []
            for filename in os.listdir(self.temp_dir):
                file_path = os.path.join(self.temp_dir, filename)
                if os.path.isfile(file_path):
                    os.remove(file_path)
                    removed_files.append(filename)
            
            if removed_files:
                print(f"All temporary files deletion completed: {len(removed_files)} files")
            
            return True
            
        except Exception as e:
            print(f"WARNING: Error deleting temporary files: {str(e)}")
            return False
    
    def _prepare_json_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Convert data to JSON serializable format."""
        json_data = {}
        
        for key, value in data.items():
            if key == 'atoms' and isinstance(value, list):
                # Convert each atom dictionary in atoms list
                json_data[key] = []
                for atom in value:
                    if isinstance(atom, dict):
                        atom_dict = {}
                        for atom_key, atom_value in atom.items():
                            if isinstance(atom_value, (int, float, str, bool, type(None))):
                                atom_dict[atom_key] = atom_value
                            else:
                                atom_dict[atom_key] = str(atom_value)
                        json_data[key].append(atom_dict)
                    else:
                        json_data[key].append(str(atom))
            elif isinstance(value, dict):
                # Dictionary type conversion
                json_data[key] = {}
                for sub_key, sub_value in value.items():
                    if isinstance(sub_value, (int, float, str, bool, type(None))):
                        json_data[key][sub_key] = sub_value
                    else:
                        json_data[key][sub_key] = str(sub_value)
            elif isinstance(value, list):
                # List type conversion
                json_data[key] = []
                for item in value:
                    if isinstance(item, (int, float, str, bool, type(None))):
                        json_data[key].append(item)
                    elif isinstance(item, dict):
                        json_data[key].append(self._prepare_json_data(item))
                    else:
                        json_data[key].append(str(item))
            elif isinstance(value, (int, float, str, bool, type(None))):
                json_data[key] = value
            else:
                json_data[key] = str(value)
        
        return json_data
    
    def _generate_summary(self, phase_name: str, data: Dict[str, Any]) -> str:
        """Generate summary text for each Phase."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        summary = f"=== {phase_name.upper()} RESULT SUMMARY ===\n"
        summary += f"Generated: {timestamp}\n\n"
        
        # Phase1 specific summary
        if phase_name.lower() == "phase1" and all(key in data for key in ['filename', 'crystal_system', 'volume', 'atoms', 'cell_params']):
            summary += f"Filename: {data['filename']}\n"
            summary += f"Crystal system: {data['crystal_system']}\n"
            summary += f"Unit cell volume: {data['volume']:.4f} Ų\n"
            summary += f"Number of atoms: {len(data['atoms'])}\n"
            summary += f"Unit cell parameters:\n"
            params = data['cell_params']
            summary += f"  a = {params['a']:.4f} Å\n"
            summary += f"  b = {params['b']:.4f} Å\n"
            summary += f"  c = {params['c']:.4f} Å\n"
            summary += f"  α = {params['alpha']:.2f}°\n"
            summary += f"  β = {params['beta']:.2f}°\n"
            summary += f"  γ = {params['gamma']:.2f}°\n"
        else:
            # General summary
            summary += f"Number of data items: {len(data)}\n"
            summary += f"Main keys: {', '.join(list(data.keys())[:10])}\n"
            if len(data) > 10:
                summary += f"... and {len(data) - 10} additional items\n"
        
        return summary 