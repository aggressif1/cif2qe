import os
import glob
from pathlib import Path

class FileReader:
    """Class responsible for file reading operations"""
    
    @staticmethod
    def read_file(file_path):
        """
        Read and return the contents of a file.
        
        Args:
            file_path (str): Path to the file to read
            
        Returns:
            str: Contents of the file
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return f.read()
        except Exception as e:
            print(f"Error reading file {file_path}: {str(e)}")
            return None
    
    @staticmethod
    def find_cif_files(input_dir="INPUT_cif"):
        """
        Find CIF files in the specified directory.
        
        Args:
            input_dir (str): Directory path to search for CIF files
            
        Returns:
            dict: Dictionary in the format {folder_name: [file_paths]}
        """
        cif_files = {}
        
        try:
            # Create input directory if it doesn't exist
            if not os.path.exists(input_dir):
                os.makedirs(input_dir)
            
            # Search all subdirectories
            for root, dirs, files in os.walk(input_dir):
                folder_name = os.path.basename(root)
                
                # Find CIF files
                cif_paths = []
                for file in files:
                    if file.lower().endswith('.cif'):
                        full_path = os.path.join(root, file)
                        cif_paths.append(full_path)
                
                if cif_paths:
                    cif_files[folder_name] = cif_paths
            
            return cif_files
            
        except Exception as e:
            print(f"Error finding CIF files: {str(e)}")
            return {} 