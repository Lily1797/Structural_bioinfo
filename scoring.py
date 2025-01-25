import os
import math
from collections import defaultdict
import argparse

# Define valid residues and valid pairs
VALID_RESIDUES = {"A", "U", "G", "C"}
VALID_PAIRS = {"AA", "AU", "AC", "AG", "CC", "CG", "CU", "GG", "GU", "UU"}

# Read the scores from the previous training (interaction profiles)
def read_scores(input_dir):
    """Read score files for each base pair."""
    scores = {}
    for pair in VALID_PAIRS:
        pair_file = os.path.join(input_dir, f"{pair}.txt")
        if os.path.exists(pair_file):
            with open(pair_file, 'r') as f:
                scores[pair] = [float(line.strip()) for line in f.readlines()]
        else:
            print(f"Warning: {pair}.txt not found.")
    return scores

# Parse the PDB file to extract C3' atom coordinates for all residues
def parse_pdb(file_path):
    """Parse a PDB file to extract residues with their C3' atom coordinates."""
    residues = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and "C3'" in line:
                chain = line[21].strip()
                res_id = int(line[22:26].strip())
                res_name = line[17:20].strip()
                if res_name in VALID_RESIDUES:
                    x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
                    residues.append((chain, res_id, res_name, (x, y, z)))
    return residues

# Compute interatomic distances and apply scoring
def compute_distances_and_scores(residues, scores):
    """Compute the interatomic distances and apply the scoring function."""
    total_score = 0
    distances = defaultdict(lambda: [0] * 20)

    for i, (chain1, res_id1, res_name1, coord1) in enumerate(residues):
        for j, (chain2, res_id2, res_name2, coord2) in enumerate(residues):
            if chain1 == chain2 and abs(res_id1 - res_id2) >= 3:  # Same chain, separated by ≥3 residues
                dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))
                bin_index = int(dist)
                if bin_index < 20:
                    pair = ''.join(sorted([res_name1, res_name2]))
                    if pair in VALID_PAIRS:
                        distances[pair][bin_index] += 1

    # Calculate total score using the observed and reference frequencies
    for pair, bins in distances.items():
        if pair in scores:
            total_counts = sum(bins)
            pair_score = 0
            for i, count in enumerate(bins):
                observed_freq = count / total_counts if total_counts > 0 else 0
                ref_freq = 1 / 20
                score = -math.log(observed_freq / ref_freq) if observed_freq > 0 else 10
                pair_score += min(score, 10)  # Cap score at 10

            total_score += pair_score

    return total_score

# Evaluate a set of RNA predicted structures
def evaluate_structures(pdb_files, score_dir):
    """Evaluate the RNA structures by calculating the total Gibbs free energy."""
    scores = read_scores(score_dir)
    total_gibbs_energy = 0
    for pdb_file in pdb_files:
        residues = parse_pdb(pdb_file)
        score = compute_distances_and_scores(residues, scores)
        total_gibbs_energy += score

    # Return the average Gibbs free energy
    num_files = len(pdb_files)
    if num_files > 0:
        return total_gibbs_energy / num_files
    else:
        return 0

# Evaluate all puzzles in the RNA-Puzzles dataset
def evaluate_all_puzzles(base_dir, score_dir):
    """Evaluate all RNA puzzles in the dataset."""
    results = {}
    for puzzle_dir in os.listdir(base_dir):
        puzzle_path = os.path.join(base_dir, puzzle_dir)
        if os.path.isdir(puzzle_path):
            pdb_files = [os.path.join(puzzle_path, f) for f in os.listdir(puzzle_path) if f.endswith(".pdb")]
            # Evaluate structures for this puzzle
            gibbs_energy = evaluate_structures(pdb_files, score_dir)
            results[puzzle_dir] = gibbs_energy
    return results

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Evaluate Gibbs free energy for RNA puzzles in a dataset."
    )
    parser.add_argument(
        "base_dir",
        type=str,
        help="Path to the base directory containing the RNA puzzles.",
    )
    parser.add_argument(
        "score_dir",
        type=str,
        help="Path to the directory containing the score files.",
    )

    # Parse arguments
    args = parser.parse_args()

    # Evaluate all puzzles in the dataset
    puzzle_results = evaluate_all_puzzles(args.base_dir, args.score_dir)

    # Print the results (Gibbs free energy for each puzzle)
    for puzzle, gibbs_energy in puzzle_results.items():
        print(f"Puzzle {puzzle}: Estimated Gibbs Free Energy = {gibbs_energy:.4f}")
