import os
import glob

class CIFReader:
    """Class for reading CIF files"""
    
    @staticmethod
    def find_cif_files(base_dir='.'):
        """
        Find CIF files in directories.
        
        Args:
            base_dir (str): Base directory to start searching from
            
        Returns:
            dict: Dictionary in the format {folder_name: [CIF file paths]}
        """
        cif_files = {}
        
        # Search for .cif files in current directory and all subdirectories
        for root, _, files in os.walk(base_dir):
            cif_in_dir = [f for f in files if f.endswith('.cif')]
            if cif_in_dir:
                folder_name = os.path.basename(root) or 'Current Directory'
                cif_files[folder_name] = [os.path.join(root, f) for f in cif_in_dir]
        
        return cif_files
    
    @staticmethod
    def read_cif_content(file_path):
        """
        Read the contents of a CIF file.
        
        Args:
            file_path (str): CIF file path
            
        Returns:
            list: List containing each line of the CIF file
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return [line.strip() for line in f if line.strip()]
        except Exception as e:
            print(f"Error reading CIF file: {e}")
            return None
    
    @staticmethod
    def display_file_menu(cif_files):
        """
        Display menu of found CIF files and get user selection.
        
        Args:
            cif_files (dict): Dictionary in the format {folder_name: [CIF file paths]}
            
        Returns:
            tuple: (selected folder name, selected file path) or list of file paths for multiple selection
        """
        if not cif_files:
            return None, None
        
        # Create a flat list of all files with their numbers
        all_files = []
        file_counter = 1
        
        for folder, paths in cif_files.items():
            for path in paths:
                all_files.append({
                    'number': file_counter,
                    'folder': folder,
                    'path': path,
                    'filename': os.path.basename(path)
                })
                file_counter += 1
        
        while True:
            choice = input("\nEnter your choice: ").strip().lower()
            
            if choice == 'x':
                return 'exit', None
            elif choice == '0':
                return 'all', None
            
            # Handle multiple file selection (e.g., "3,4,7")
            elif ',' in choice:
                try:
                    selected_numbers = [int(x.strip()) for x in choice.split(',')]
                    selected_files = []
                    
                    for num in selected_numbers:
                        if 1 <= num <= len(all_files):
                            selected_files.append(all_files[num - 1])
                        else:
                            print(f"Invalid file number: {num}")
                            continue
                    
                    if selected_files:
                        return 'multiple', selected_files
                    else:
                        print("No valid files selected. Please try again.")
                        continue
                        
                except ValueError:
                    print("Invalid format. Use comma-separated numbers (e.g., 3,4,7)")
                    continue
            
            # Handle single file selection
            else:
                try:
                    file_num = int(choice)
                    if 1 <= file_num <= len(all_files):
                        selected_file = all_files[file_num - 1]
                        return selected_file['folder'], selected_file['path']
                    else:
                        print(f"Invalid file number. Please enter a number between 1 and {len(all_files)}")
                        continue
                except ValueError:
                    print("Invalid input. Please enter a number, comma-separated numbers, '0' for all, or 'x' to exit")
                    continue 