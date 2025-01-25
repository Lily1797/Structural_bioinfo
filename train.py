import os
import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import numpy as np

# Define valid residues and valid pairs
VALID_RESIDUES = {"A", "U", "G", "C"}
VALID_PAIRS = {"AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"}

def initialize_counts():
    """Initialize a list of zeros for distance bins (20 bins)."""
    return np.zeros(20)

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

def compute_distances(residues):
    """Compute interatomic distances for valid residue pairs separated by ≥3 positions."""
    distances = defaultdict(initialize_counts)
    for i, (chain1, res_id1, res_name1, coord1) in enumerate(residues):
        for j, (chain2, res_id2, res_name2, coord2) in enumerate(residues):
            if chain1 == chain2 and abs(res_id1 - res_id2) >= 3:  # Same chain, separated by ≥3 residues
                dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))
                bin_index = int(dist)
                if bin_index < 20:  # Only consider distances up to 20 Å
                    pair = ''.join(sorted([res_name1, res_name2]))
                    if pair in VALID_PAIRS:  # Only count valid pairs
                        distances[pair][bin_index] += 1
    return distances

def aggregate_counts(all_distances):
    """Aggregate counts across all PDB files."""
    combined_counts = defaultdict(lambda: [0] * 20)
    for distances in all_distances:
        for pair, bins in distances.items():
            if pair in VALID_PAIRS:  # Only aggregate valid pairs
                for bin_idx, count in enumerate(bins):
                    combined_counts[pair][bin_idx] += count
    return combined_counts

def calculate_scores(combined_counts):
    """Calculate scores for each base pair and distance bin."""
    scores = {}
    for pair, counts in combined_counts.items():
        total_counts = sum(counts)
        pair_scores = []
        for i, count in enumerate(counts):
            observed_freq = count / total_counts if total_counts > 0 else 0
            ref_freq = 1 / 20  # Uniform reference frequency for now
            score = -math.log(observed_freq / ref_freq) if observed_freq > 0 else 10
            pair_scores.append(min(score, 10))  # Cap score at 10
        scores[pair] = pair_scores
    return scores


def write_scores(scores, output_dir):
    """Write scores to files, one per base pair."""
    os.makedirs(output_dir, exist_ok=True)
    for pair in VALID_PAIRS:
        if pair in scores and any(score > 0 for score in scores[pair]):  # Check if the pair has any valid data
            with open(os.path.join(output_dir, f"{pair}.txt"), 'w') as f:
                for score in scores[pair]:
                    f.write(f"{score:.4f}\n")


def process_single_file(pdb_file):
    """Process a single PDB file to compute distances."""
    residues = parse_pdb(pdb_file)
    return compute_distances(residues)

def process_pdb_files_parallel(pdb_files, output_dir):
    """Process multiple PDB files in parallel and compute scores."""
    with ProcessPoolExecutor() as executor:
        all_distances = list(executor.map(process_single_file, pdb_files))
    combined_counts = aggregate_counts(all_distances)
    scores = calculate_scores(combined_counts)
    write_scores(scores, output_dir)

if __name__ == "__main__":
    # Define input and output directories
    pdb_directory = "/Structural_bioinfo/data"
    output_directory = "/Structural_bioinfo/output"

    # List all PDB files in the directory
    pdb_files = [os.path.join(pdb_directory, f) for f in os.listdir(pdb_directory) if f.endswith(".pdb")]

    # Process PDB files and generate output
    process_pdb_files_parallel(pdb_files, output_directory)
