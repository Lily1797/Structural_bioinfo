import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

# Define valid base pairs
VALID_PAIRS = ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]

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

def interpolate_scores(score_list):
    """Interpolate the scores within each distance bin."""
    # Create the distance bins (1 to 20 Å)
    distance_bins = list(range(1, 21))
    interpolated_scores = []

    for i in range(len(score_list) - 1):
        x1, y1 = distance_bins[i], score_list[i]
        x2, y2 = distance_bins[i + 1], score_list[i + 1]

        # Interpolation for distances between x1 and x2
        num_points = 10  # 10 points between each bin
        for j in range(num_points):
            x = x1 + (x2 - x1) * j / (num_points - 1)  # Linearly interpolate between x1 and x2
            y = y1 + (y2 - y1) * j / (num_points - 1)  # Linearly interpolate between y1 and y2
            interpolated_scores.append((x, y))

    return interpolated_scores

def plot_interaction_profiles(scores, output_dir):
    """Plot interaction profiles for each base pair."""
    os.makedirs(output_dir, exist_ok=True)

    # Set up the plot for each base pair
    for pair, score_list in scores.items():
        plt.figure(figsize=(8, 6))

        # Interpolate the scores
        interpolated_data = interpolate_scores(score_list)

        # Plot the interaction profile (Distance vs Score)
        x_vals = [x for x, _ in interpolated_data]
        y_vals = [y for _, y in interpolated_data]

        plt.plot(x_vals, y_vals, marker='o', linestyle='-', label=pair)
        plt.xlabel("Distance (Å)", fontsize=12)
        plt.ylabel("Pseudo-Energy Score", fontsize=12)
        plt.title(f"Interaction Profile for Base Pair {pair}", fontsize=14)
        plt.grid(True)
        plt.legend()

        # Save the plot to a file
        plt.savefig(os.path.join(output_dir, f"{pair}_interaction_profile.png"))
        plt.close()

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Plot interaction profiles for RNA base pairs from score files."
    )
    parser.add_argument(
        "input_directory",
        type=str,
        help="Path to the directory containing score files.",
    )
    parser.add_argument(
        "output_directory",
        type=str,
        help="Path to the directory where plots will be saved.",
    )

    # Parse arguments
    args = parser.parse_args()

    # Read scores from the input directory
    scores = read_scores(args.input_directory)

    # Plot and save the figures
    plot_interaction_profiles(scores, args.output_directory)
