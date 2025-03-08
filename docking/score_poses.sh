#!/bin/bash

# NOTE: DO NOT TRUST THE RANKINGS IMPLICITY, YOU SHOULD STILL CHECK STRUCTURES

module load rosetta # load rosetta

# Define output files
score_file="score.sc"   # Rosetta generates this file
sorted_output="rosetta_scores_sorted.txt"

# Ensure the output file starts fresh
echo "<U+1F52C> Running Rosetta scoring for all structures..."
rm -f $score_file  # Remove old score.sc if it exists

# Run Rosetta scoring on all PDBs
score_jd2.mpi.linuxgccrelease -in:file:s ringplane_translation_*.pdb -score:weights score12 -out:file:scorefile $score_file

echo "<U+2705> Scoring complete! Now extracting and sorting scores..."

# Parse score.sc and extract total_score, fa_rep (clashes), and hbond_bb_sc (backbone H-bonds)
awk '
    BEGIN { OFS="\t"; print "PDB_File", "Total_Score", "fa_rep", "hbond_bb_sc" }
    NR > 2 { print $NF, $2, $11, $14 }  # Extract PDB name (last col), total_score, fa_rep, hbond_bb_sc
' $score_file | sort -k2,2n > $sorted_output

echo "<U+2705> Sorting complete! Results saved in $sorted_output"
echo "<U+1F4CA> Top 5 lowest-scoring structures:"
head -6 $sorted_output  # Show top 5 results
