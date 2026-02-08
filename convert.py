import os
import subprocess
import argparse

def convert_pdbqt_to_pdb(input_folder, output_folder):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Iterate through all pdbqt files in the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".pdbqt"):
            input_file = os.path.join(input_folder, file_name)
            output_file = os.path.join(output_folder, file_name.replace(".pdbqt", ".pdb"))
            
            # Use OpenBabel to convert pdbqt to pdb
            try:
                subprocess.run(["obabel", input_file, "-O", output_file], check=True)
                print(f"Converted: {input_file} -> {output_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error converting {input_file}: {e}")
    
    print("Conversion complete.")

def main():
    # Initialize argparse
    parser = argparse.ArgumentParser(description="Convert .pdbqt files to .pdb files")
    
    # Add command-line arguments
    parser.add_argument("input_folder", type=str, help="Path to the folder containing .pdbqt files")
    parser.add_argument("output_folder", type=str, help="Path to the folder where .pdb files will be saved")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the conversion function
    convert_pdbqt_to_pdb(args.input_folder, args.output_folder)

if __name__ == "__main__":
    main()

