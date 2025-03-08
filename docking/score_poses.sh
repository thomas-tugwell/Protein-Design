#!/bin/bash

# Gain access to rosetta on the cluster
module load rosetta 

# Define output file for scores
output_file="rosetta_scores.txt"
echo "PDB_File Total_Score" > $output_file

# Loop through all PDBs
for pdb in ringplane_translation_*.pdb; do
    echo "Scoring $pdb..."

    # Run Rosetta scoring
    score_jd2.mpi.linuxgccrelease -s "$pdb" -score:weights score12 > "${pdb%.pdb}_score.log"

    # Extract the total score
    score=$(grep "pose " "${pdb%.pdb}_score.log" | awk '{print $NF}')

    # Save result
    echo "$pdb $score" >> $output_file
done

echo "<U+2705> Scoring complete! "
