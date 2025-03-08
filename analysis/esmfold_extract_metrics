import re
import argparse

# Set up command line argument parsing
parser = argparse.ArgumentParser(description="Process log file to find top 10 sequences with highest pLDDT, pTM, and combined score.")
parser.add_argument("file", help="Path to the log file to process")
args = parser.parse_args()

# Initialize lists to store sequences and scores
plddt_sequences = []
ptm_sequences = []
combined_sequences = []

# Regex pattern to extract sequence, pLDDT, and pTM
pattern = re.compile(r"sequence_(\d+).+pLDDT (\d+\.\d+).+pTM (\d+\.\d+)")

# Open and read the log file
with open(args.file, 'r') as file:
    for line in file:
        match = pattern.search(line)
        if match:
            seq_id = match.group(1)
            plddt = float(match.group(2))
            ptm = float(match.group(3))

            # Normalize pLDDT to a scale of 0-1
            normalized_plddt = plddt / 100

            # Calculate combined score (sum or average)
            combined_score = (normalized_plddt + ptm) / 2  # Average of normalized pLDDT and pTM

            # Add sequences and scores to their respective lists
            plddt_sequences.append((seq_id, plddt))
            ptm_sequences.append((seq_id, ptm))
            combined_sequences.append((seq_id, combined_score))

# Sort the lists in descending order by score
plddt_sequences.sort(key=lambda x: x[1], reverse=True)
ptm_sequences.sort(key=lambda x: x[1], reverse=True)
combined_sequences.sort(key=lambda x: x[1], reverse=True)

# Print the top 10 sequences for each category
print("Top 10 Sequences by pLDDT:")
for seq, score in plddt_sequences[:10]:
    print(f"Sequence {seq} with pLDDT {score}")

print("\nTop 10 Sequences by pTM:")
for seq, score in ptm_sequences[:10]:
    print(f"Sequence {seq} with pTM {score}")

print("\nTop 10 Sequences by Combined Score:")
for seq, score in combined_sequences[:10]:
    print(f"Sequence {seq} with combined score {score}")
