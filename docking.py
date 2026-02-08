import argparse
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from vina import Vina
import os
import numpy as np
from multiprocessing import Pool, cpu_count

# Set up argparse to handle command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Dock ligands to a protein using AutoDock Vina and RDKit")
    parser.add_argument('ligand1', type=str, help="Path to the ligand file in SDF format (can contain multiple ligands)")
    parser.add_argument('ligand2', type=str, help="Path to the second ligand in PDBQT format")
    parser.add_argument('protein', type=str, help="Path to the protein in PDBQT format")
    parser.add_argument('output', type=str, default="Ligand_vina_out", help="Output folder for docked poses")
    parser.add_argument('n_poses', type=int, default=10, help="Number of poses to dock (default: 10)")
    parser.add_argument('box_size', type=float, nargs=3, default=[20, 20, 20],
                        help="Box size for docking as height, width, and length (default: 20, 20, 20)")
    parser.add_argument('exhaustiveness', type=int, default=32, help="Exhaustiveness of docking (default: 32)")
    parser.add_argument('score_output', type=str, default="docking_scores.txt", help="File to save docking scores")
    parser.add_argument('max_scores', type=int, default=10, help="Maximum number of docking scores to select ")
    return parser.parse_args()

# Docking function to be used in multiprocessing
def dock_ligand(ligand_info):
    ligand_mol, ligand_id, ligand2_path, protein_path, box_size, exhaustiveness, n_poses, output_folder = ligand_info
    
    try:
        # Initialize Vina object
        v = Vina(sf_name='vina')
        
        # Load protein content from file
        with open(protein_path, 'r') as protein_file:
            protein_content = protein_file.read()

        # Set receptor
        v.set_receptor(protein_path)
        
        # Calculate centroid of ligand
        centroid = ComputeCentroid(ligand_mol.GetConformer())
        
        # Set ligand2 for docking (assumed to be the same for all ligands)
        v.set_ligand_from_file(ligand2_path)
        
        # Compute Vina maps
        v.compute_vina_maps(center=[centroid.x, centroid.y, centroid.z], box_size=box_size)
        
        # Perform energy minimization and docking
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        
        # Save the docking poses
        combined_file = os.path.join(output_folder, f"{ligand_id}.pdbqt")
        with open(combined_file, 'w') as output_combined:
            docked_poses = v.poses(n_poses=n_poses).split("ENDMDL")
            for pose_num, pose in enumerate(docked_poses):
                if pose.strip():
                    output_combined.write(f"MODEL {pose_num + 1}\n")
                    output_combined.write(protein_content + "\n")
                    output_combined.write(pose.strip() + "\nENDMDL\n")
        
        # Extract docking scores
        docking_scores = v.energies()
        docking_scores = np.array(docking_scores)
        result_scores = []
        
        if docking_scores.size > 0:
            for pose_num, score in enumerate(docking_scores):
                if score.size > 0:
                    pose_id = f"{ligand_id}_pose_{pose_num + 1}"
                    first_score = score[0]
                    result_scores.append((pose_id, first_score))
        
        return result_scores
    
    except Exception as e:
        print(f"Error docking ligand {ligand_id}: {e}")
        return []

# Main docking process with multiprocessing
def main():
    args = parse_arguments()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Load ligands from the SDF file (multiple ligands possible)
    ligands = Chem.SDMolSupplier(args.ligand1, removeHs=False)
    if ligands is None or len(ligands) == 0:
        raise ValueError("No ligands found in the SDF file.")

    # Prepare ligand info for multiprocessing (pass only necessary information)
    ligand_info = []
    for i, lig in enumerate(ligands):
        if lig is None:
            print(f"Skipping invalid ligand at index {i}.")
            continue
        
        ligand_id = lig.GetProp('_Name') if lig.HasProp('_Name') else f"ligand_{i}"
        ligand_info.append((lig, ligand_id, args.ligand2, args.protein, args.box_size, args.exhaustiveness, args.n_poses, args.output))
    
    # Start multiprocessing
    with Pool(processes=cpu_count()) as pool:
        all_results = pool.map(dock_ligand, ligand_info)

    # Flatten all results into a single list
    all_scores = [score for sublist in all_results for score in sublist]

    # Sort the scores from lowest to highest
    all_scores.sort(key=lambda x: x[1])

    # Select the top `max_scores` poses
    top_scores = all_scores[:args.max_scores]

    # Write selected docking scores to the output file
    with open(args.score_output, 'w') as score_file:
        for pose_id, score in top_scores:
            score_file.write(f"{pose_id}\t{score:.3f}\n")
            print(f"{pose_id} docking score: {score:.3f}")

if __name__ == '__main__':
    main()

