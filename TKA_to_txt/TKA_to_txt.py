import os

def replace_extension(folder_path):
    # Walk through the folder and its subfolders
    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            # Check if the file has a .TKA extension
            if file_name.endswith('.TKA'):
                # Create the new file name with .txt extension
                new_file_name = file_name[:-4] + '.txt'
                # Get the full paths for the old and new file names
                old_file_path = os.path.join(root, file_name)
                new_file_path = os.path.join(root, new_file_name)
                # Rename the file
                os.rename(old_file_path, new_file_path)
                print(f'Renamed: {old_file_path} -> {new_file_path}')




# Specify the folder path
folder_path = 'C:/Users/arist/Desktop/WS24_25/Fortgeschrittenenpraktikum/Versuche/T01 (Detektoren & Strahlenschutz)/Daten'

# Call the function to replace extensions
replace_extension(folder_path)


