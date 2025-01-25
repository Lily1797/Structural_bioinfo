# Structural_bioinfo
Creation of an objective function for the RNA folding problem

### Context
For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold among the astronomically large number of possible conformations. The native fold being the
one with the lowest Gibbs free energy, the objective function should be an estimator of this energy.

For this practical course, we implemented 3 Python scripts to:
- train the objective function, using interatomic distance distributions that are computed from a dataset of known (i.e., experimentally determined) 3D structures;
- plot the scoring profiles, i.e. the score (or estimated Gibbs free energy) as a function of the interatomic distance;
- use the objective function to evaluate predicted structures from the RNA-Puzzles dataset

## Training
#### Data:
We retrieved a dataset of known 3D structures from the [RCSB Protein Data Bank](https://www.rcsb.org/) with the following
QUERY: Polymer Entity Type = "RNA" AND Polymer Entity Sequence Length > 2000 AND Experimental Method = "X-RAY DIFFRACTION", and saved to [data directory](data/).
#### Train the objective function
To analyze RNA structures in PDB format, we first used the *train_all.py* script to compute interatomic distances between residue pairs.
```
usage: script_name.py [-h] pdb_directory output_directory

Process PDB files to compute pseudo-energy scores for RNA base pairs.

positional arguments:
  pdb_directory      Path to the directory containing PDB files.
  output_directory   Path to the directory where output scores will be saved.

optional arguments:
  -h, --help         Show this help message and exit.
```
The scores for each residue pair were saved into separate text files in the [output directory](output/).
#### Plot the scoring profiles
Subsequently, we used the *plot.py* script to generate interaction profile plots for RNA base pairs based on computed scores.
```
usage: script_name.py [-h] input_directory output_directory

Plot interaction profiles for RNA base pairs from score files.

positional arguments:
  input_directory    Path to the directory containing score files.
  output_directory   Path to the directory where plots will be saved.

optional arguments:
  -h, --help         Show this help message and exit.
```
This script reads pseudo-energy score files for each valid base pair from the output directory. In the training step, a single score value was assigned to each distance bin. However, the actual score values for distances within the same bin can vary (e.g., in the bin [6-7], the score at 6.1 may differ from the score at 6.9). To address this, the script performs linear interpolation between the bins, generating a smooth interaction profile curve. It then creates plots of "Distance (Å)" vs. "Pseudo-Energy Score" for each base pair, saving each plot as a PNG file in the [plot directory](plot/).

## Scoring
Predicted structures: [RNA puzzles](https://github.com/RNA-Puzzles/standardized_dataset.git).
We downloaded all 23 structures and stored them in the [puzzles directory](standardized_dataset/). 
To evaluate the predicted structures from the RNA-Puzzles dataset, we used the *scoring.py* script, which shares some similarities with the training script. This script calculates all interatomic distances for a given structure, applying the same thresholds (20 Å and 𝑖, 𝑖+4). For each distance, a scoring value is determined through linear interpolation. By summing all these scores, the script computes the estimated Gibbs free energy of the evaluated RNA conformation.
```
usage: script_name.py [-h] base_dir score_dir

Evaluate Gibbs free energy for RNA puzzles in a dataset.

positional arguments:
  base_dir    Path to the base directory containing the RNA puzzles.
  score_dir   Path to the directory containing the score files.

optional arguments:
  -h, --help  Show this help message and exit.
```
Results
```
Puzzle rp02: Estimated Gibbs Free Energy = 1372.3060
Puzzle rp11: Estimated Gibbs Free Energy = 1042.2684
Puzzle rp10: Estimated Gibbs Free Energy = 787.1789
Puzzle rp12: Estimated Gibbs Free Energy = 777.0650
Puzzle rp05: Estimated Gibbs Free Energy = 692.4727
Puzzle rp07: Estimated Gibbs Free Energy = 715.1541
Puzzle rp09: Estimated Gibbs Free Energy = 890.7312
Puzzle rp18: Estimated Gibbs Free Energy = 982.4398
Puzzle rp03: Estimated Gibbs Free Energy = 941.8519
Puzzle rp04: Estimated Gibbs Free Energy = 767.4005
Puzzle rp20: Estimated Gibbs Free Energy = 1095.2514
Puzzle rp08: Estimated Gibbs Free Energy = 827.3233
Puzzle rp15: Estimated Gibbs Free Energy = 1037.9824
Puzzle rp17: Estimated Gibbs Free Energy = 901.8072
Puzzle rp14_free: Estimated Gibbs Free Energy = 944.9496
Puzzle rp13: Estimated Gibbs Free Energy = 999.4704
Puzzle rp06: Estimated Gibbs Free Energy = 722.3957
Puzzle rp01: Estimated Gibbs Free Energy = 1276.2236
Puzzle rp24: Estimated Gibbs Free Energy = 867.8210
Puzzle rp19: Estimated Gibbs Free Energy = 1232.4318
Puzzle rp16_TBA: Estimated Gibbs Free Energy = 0.0000
Puzzle rp14_bound: Estimated Gibbs Free Energy = 987.3209
Puzzle .git: Estimated Gibbs Free Energy = 0.0000
Puzzle rp21: Estimated Gibbs Free Energy = 950.9333
```

### Directory Structure
```
Structural_bioinfo/
├── data/
│   ├── structure1.pdb
│   ├── structure2.pdb
│   ├── ...
├── output/
│   ├── AA.txt
│   ├── AU.txt
│   ├── ...
├── plot/
│   ├── AA_interaction_profile.png
│   ├── AU_interaction_profile.png
│   ├── ...
├── standardized_dataset/
│   ├── puzzle_1/
│   │   ├── structure1.pdb
│   │   ├── structure2.pdb
│   ├── puzzle_2/
│   │   ├── structure1.pdb
│   │   ├── structure2.pdb
│   ├── ...
```
