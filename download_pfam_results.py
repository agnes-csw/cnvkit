#!/usr/bin/env python3
"""
Script to download all PFAM/HMMER results from EBI HMMER web service.
Automatically fetches all domain predictions for each transcript.

USAGE:
    python download_pfam_results.py

    Or with a custom job ID:
    python download_pfam_results.py YOUR_JOB_ID

The script will create a 'pfam_results' directory with:
  - all_results.json: Complete raw API response
  - domain_summary.tsv: Summary table of all hits
  - all_domains.tsv: Detailed domain predictions
  - domain_architectures.json: Architecture data (if available)
  - Individual JSON files per hit (optional)
"""

import requests
import json
import time
import sys
import argparse
from pathlib import Path

# Default job ID (replace with your own or pass as argument)
DEFAULT_JOB_ID = "50b4593e-293e-4875-b211-1454cba3cbb8"

# API endpoints
BASE_URL = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"

def get_result_url(job_id):
    return f"{BASE_URL}/result/{job_id}"

def get_architecture_url(job_id):
    return f"{BASE_URL}/architecture/{job_id}"

# Output directory
OUTPUT_DIR = Path("pfam_results")


def fetch_json(url, params=None, retries=3):
    """Fetch JSON data from URL with error handling and retries."""
    headers = {
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }

    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, params=params, timeout=60)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                wait_time = 2 ** attempt
                print(f"  Attempt {attempt + 1} failed, retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"Error fetching {url}: {e}")
                return None


def download_all_results(job_id, output_dir, save_individual=False):
    """Download all HMMER results for the job."""
    output_dir.mkdir(exist_ok=True)

    result_url = get_result_url(job_id)
    print(f"Fetching main results from job: {job_id}")
    print(f"URL: {result_url}")

    # First, get the main results
    results = fetch_json(result_url)

    if results is None:
        print("Failed to fetch main results. Trying alternative approach...")
        # Try the score endpoint format
        alt_url = f"https://www.ebi.ac.uk/Tools/hmmer/results/{job_id}/score"
        print(f"Note: The web interface URL is: {alt_url}")
        return None

    # Save the main results JSON
    main_results_file = output_dir / "all_results.json"
    with open(main_results_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved main results to: {main_results_file}")

    # Parse and extract individual transcript results
    if 'results' in results:
        hits = results.get('results', {}).get('hits', [])
    elif 'hits' in results:
        hits = results.get('hits', [])
    else:
        hits = results if isinstance(results, list) else []

    print(f"Found {len(hits)} hits/transcripts")

    # Create a summary TSV file
    summary_file = output_dir / "domain_summary.tsv"
    with open(summary_file, 'w') as f:
        f.write("transcript_id\taccession\tname\tdescription\tevalue\tscore\tbias\tdomains\n")

        for i, hit in enumerate(hits):
            # Extract hit information (structure may vary)
            if isinstance(hit, dict):
                acc = hit.get('acc', hit.get('accession', 'N/A'))
                name = hit.get('name', 'N/A')
                desc = hit.get('desc', hit.get('description', 'N/A'))
                evalue = hit.get('evalue', 'N/A')
                score = hit.get('score', 'N/A')
                bias = hit.get('bias', 'N/A')
                domains = hit.get('domains', [])
                num_domains = len(domains) if isinstance(domains, list) else 'N/A'

                f.write(f"{name}\t{acc}\t{name}\t{desc}\t{evalue}\t{score}\t{bias}\t{num_domains}\n")

                # Save individual hit details (optional)
                if save_individual:
                    safe_name = name.replace('/', '_').replace('\\', '_')[:50]
                    hit_file = output_dir / f"hit_{i+1:04d}_{safe_name}.json"
                    with open(hit_file, 'w') as hf:
                        json.dump(hit, hf, indent=2)

    print(f"Saved summary to: {summary_file}")

    # Also try to get domain architecture data
    print("\nFetching domain architecture data...")
    arch_url = get_architecture_url(job_id)
    arch_data = fetch_json(arch_url)
    if arch_data:
        arch_file = output_dir / "domain_architectures.json"
        with open(arch_file, 'w') as f:
            json.dump(arch_data, f, indent=2)
        print(f"Saved architecture data to: {arch_file}")

    return results, hits


def create_domains_table(hits, output_dir):
    """Create a detailed table of all domain predictions."""
    if not hits:
        return

    domains_file = output_dir / "all_domains.tsv"

    with open(domains_file, 'w') as f:
        f.write("sequence_name\tsequence_acc\tdomain_name\tdomain_acc\tdomain_desc\t")
        f.write("seq_start\tseq_end\thmm_start\thmm_end\tievalue\tcevalue\tscore\n")

        for hit in hits:
            if not isinstance(hit, dict):
                continue

            seq_name = hit.get('name', 'N/A')
            seq_acc = hit.get('acc', hit.get('accession', 'N/A'))

            domains = hit.get('domains', [])
            for dom in domains:
                if isinstance(dom, dict):
                    dom_name = dom.get('alihmmacc', dom.get('name', 'N/A'))
                    dom_acc = dom.get('alihmmname', dom.get('acc', 'N/A'))
                    dom_desc = dom.get('alihmmdesc', dom.get('desc', 'N/A'))
                    seq_from = dom.get('alisqfrom', dom.get('seq_from', 'N/A'))
                    seq_to = dom.get('alisqto', dom.get('seq_to', 'N/A'))
                    hmm_from = dom.get('alihmmfrom', dom.get('hmm_from', 'N/A'))
                    hmm_to = dom.get('alihmmto', dom.get('hmm_to', 'N/A'))
                    ievalue = dom.get('ievalue', 'N/A')
                    cevalue = dom.get('cevalue', 'N/A')
                    score = dom.get('bitscore', dom.get('score', 'N/A'))

                    f.write(f"{seq_name}\t{seq_acc}\t{dom_name}\t{dom_acc}\t{dom_desc}\t")
                    f.write(f"{seq_from}\t{seq_to}\t{hmm_from}\t{hmm_to}\t{ievalue}\t{cevalue}\t{score}\n")

    print(f"Saved detailed domains table to: {domains_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Download PFAM/HMMER results from EBI HMMER web service',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python download_pfam_results.py
    python download_pfam_results.py YOUR-JOB-UUID-HERE
    python download_pfam_results.py 50b4593e-293e-4875-b211-1454cba3cbb8 -o my_results
        """
    )
    parser.add_argument('job_id', nargs='?', default=DEFAULT_JOB_ID,
                        help='HMMER job UUID (default: %(default)s)')
    parser.add_argument('-o', '--output', type=Path, default=OUTPUT_DIR,
                        help='Output directory (default: %(default)s)')
    parser.add_argument('--individual', action='store_true',
                        help='Save individual JSON files per hit')

    args = parser.parse_args()

    job_id = args.job_id
    output_dir = args.output

    print("=" * 60)
    print("HMMER/PFAM Results Downloader")
    print("=" * 60)
    print(f"\nJob ID: {job_id}")
    print(f"Output directory: {output_dir.absolute()}\n")

    result = download_all_results(job_id, output_dir, args.individual)

    if result is not None:
        results, hits = result
        create_domains_table(hits, output_dir)
        print("\n" + "=" * 60)
        print("Download complete!")
        print(f"Check the '{output_dir}' directory for all files.")
        print("=" * 60)
        print("\nFiles created:")
        print("  - all_results.json: Complete raw API response")
        print("  - domain_summary.tsv: Summary of all transcript hits")
        print("  - all_domains.tsv: Detailed domain predictions")
        print("  - domain_architectures.json: Architecture data (if available)")
    else:
        print("\n" + "=" * 60)
        print("Could not fetch results via API.")
        print("\nAlternative approaches:")
        print("1. Run this script on your local machine (not behind a proxy)")
        print("2. Use the browser console script (download_pfam_browser.js)")
        print("3. Use curl commands manually (shown below)")
        print("=" * 60)

        result_url = get_result_url(job_id)
        arch_url = get_architecture_url(job_id)

        print("\nTry these curl commands manually:")
        print(f"\ncurl -H 'Accept: application/json' '{result_url}' -o results.json")
        print(f"\ncurl -H 'Accept: application/json' '{arch_url}' -o architecture.json")


if __name__ == "__main__":
    main()
