#!/usr/bin/env python3
import re
import sys

def parse_log(logfile):
    """
    Parse the log file and return a list of dictionaries for each query.
    Each dictionary will contain:
      - 'sequence': the query sequence name
      - 'pLDDT': pLDDT value from the rank_001 line (as float)
      - 'pTM': pTM value from the rank_001 line (as float)
    """
    # Regex patterns:
    # Pattern to capture the query line: sequence is the text between ": " and " (length"
    query_pattern = re.compile(r"Query\s+\d+/\d+:\s+(.+?)\s+\(length")
    # Pattern to capture the rank_001 line: we expect something like:
    # rank_001_alphafold2_ptm_model_3_seed_000 pLDDT=87.2 pTM=0.32
    rank_pattern = re.compile(r"rank_001_\S+\s+pLDDT=([\d.]+)\s+pTM=([\d.]+)")

    entries = []
    current_sequence = None
    current_metrics = None

    with open(logfile, "r") as f:
        for line in f:
            line = line.strip()
            # Check for query start
            q_match = query_pattern.search(line)
            if q_match:
                # If we were processing a previous query and already captured metrics, save it.
                if current_sequence and current_metrics:
                    entry = {
                        "sequence": current_sequence,
                        "pLDDT": current_metrics["pLDDT"],
                        "pTM": current_metrics["pTM"]
                    }
                    entries.append(entry)
                # Start a new query block.
                current_sequence = q_match.group(1)
                current_metrics = {}  # reset metrics for this query
                continue

            # Check for rank_001 line within the current query block.
            r_match = rank_pattern.search(line)
            if r_match and current_sequence:
                current_metrics["pLDDT"] = float(r_match.group(1))
                current_metrics["pTM"] = float(r_match.group(2))
                # We assume one rank_001 per query; once found, we can continue.
                continue

    # Capture the last query block if present
    if current_sequence and current_metrics:
        entry = {
            "sequence": current_sequence,
            "pLDDT": current_metrics["pLDDT"],
            "pTM": current_metrics["pTM"]
        }
        entries.append(entry)

    return entries

def write_rankings(entries, outfile):
    # Sort entries by pTM descending (highest to lowest)
    sorted_entries = sorted(entries, key=lambda x: x["pTM"], reverse=True)
    with open(outfile, "w") as f:
        for entry in sorted_entries:
            f.write(f"{entry['sequence']}: pTM={entry['pTM']} pLDDT={entry['pLDDT']}\n")
    print(f"Rankings written to {outfile}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python parse_rankings.py <logfile>")
        sys.exit(1)
    logfile = sys.argv[1]
    entries = parse_log(logfile)
    if not entries:
        print("No query entries found in the log.")
        sys.exit(1)
    write_rankings(entries, "rankings.txt")

if __name__ == "__main__":
    main()

