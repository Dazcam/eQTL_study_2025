import re
import csv
import json
import subprocess
import os
import argparse
from datetime import datetime

def parse_snakemake_log(log_file):
    """Extracts job details and run information from a Snakemake log file."""
    jobs = []
    with open(log_file, 'r') as f:
        lines = f.readlines()

    current_job = {}
    for i, line in enumerate(lines):
        if line.startswith("Reason:"):
            # Extract the run information from the Reason line
            match = re.search(r'scanpy_(.+?)\.html', line)
            if match:
                current_job['Run'] = match.group(1)
        if "Submitted batch job" in line:
            # Extract Job ID
            match = re.search(r"Submitted batch job (\d+)", line)
            if match:
                current_job['JOB_ID'] = match.group(1)
            # Extract and convert Date from previous lines (if available)
            for prev_line in reversed(lines[:i]):
                if prev_line.startswith("["):
                    date_match = re.search(r"\[(.+?)\]", prev_line)
                    if date_match:
                        try:
                            # Convert to YYYY-MM-DD format
                            date = datetime.strptime(date_match.group(1), "%a %b %d %H:%M:%S %Y")
                            current_job['Date'] = date.strftime("%Y-%m-%d")
                        except ValueError:
                            current_job['Date'] = date_match.group(1)
                        break
            jobs.append(current_job)
            current_job = {}  # Reset for the next job
    return jobs

def fetch_seff_output(job_id):
    """Runs the seff command for a given job ID and parses its output."""
    try:
        seff_output = subprocess.check_output(["seff", job_id], text=True)
        details = {}
        for line in seff_output.splitlines():
            if line.startswith("Cores per node:"):
                details['Cores'] = int(line.split(":")[1].strip())
            elif line.startswith("Nodes:"):
                details['Nodes'] = int(line.split(":")[1].strip())
            elif line.startswith("State:"):
                details['State'] = line.split(":")[1].strip()
            elif line.startswith("CPU Efficiency:"):
                efficiency = re.search(r"([\d\.]+)%", line)
                if efficiency:
                    details['CPU_efficiency_pct'] = float(efficiency.group(1))
            elif line.startswith("Memory Efficiency:"):
                efficiency = re.search(r"([\d\.]+)%", line)
                if efficiency:
                    details['Mem_efficiency_pct'] = float(efficiency.group(1))
            elif line.startswith("Memory Utilized:"):
                mem_utilized = re.search(r"([\d\.]+) GB", line)
                if mem_utilized:
                    details['Mem_utilized_GB'] = float(mem_utilized.group(1))
            elif line.startswith("Job Wall-clock time:"):
                elapsed_time = line.split(":", 1)[1].strip()
                details['Elapsed_time'] = elapsed_time
        return details
    except subprocess.CalledProcessError:
        print(f"Failed to fetch seff output for job ID {job_id}")
        return {}

def combine_data(snakemake_jobs, run_details=None):
    """Combines Snakemake job details with seff output."""
    combined_data = []
    for job in snakemake_jobs:
        job_id = job.get('JOB_ID')
        if job_id:
            seff_details = fetch_seff_output(job_id)
            combined_data.append({**job, **seff_details, 'run_details': run_details})
        else:
            combined_data.append({**job, 'run_details': run_details})  # Add Snakemake details if no job ID
    return combined_data

def write_csv(data, output_file):
    """Writes combined data to a CSV file, appending if the file exists."""
    file_exists = os.path.exists(output_file)
    with open(output_file, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=[
            'Date', 'JOB_ID', 'Run', 'Cores', 'Nodes', 'State', 
            'Elapsed_time', 'CPU_efficiency_pct', 'Mem_utilized_GB', 
            'Mem_efficiency_pct', 'run_details'
        ])
        if not file_exists:
            writer.writeheader()  # Write headers only if the file is new
        writer.writerows(data)

def write_json(data, output_file):
    """Writes combined data to a JSON file."""
    with open(output_file, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Parse Snakemake logs and fetch seff output.")
    parser.add_argument("input", help="Path to the Snakemake log file.")
    parser.add_argument("output", help="Path to the output file (.csv or .json).")
    parser.add_argument("--details", help="Optional additional details for the 'run_details' column.", default=None)
    args = parser.parse_args()

    # Parse the log file
    snakemake_jobs = parse_snakemake_log(args.input)

    # Combine data with seff output and optional details
    combined_data = combine_data(snakemake_jobs, run_details=args.details)

    # Write to the appropriate output format
    if args.output.endswith(".csv"):
        write_csv(combined_data, args.output)
    elif args.output.endswith(".json"):
        write_json(combined_data, args.output)
    else:
        print("Unsupported file format. Please use .csv or .json.")

if __name__ == "__main__":
    main()

