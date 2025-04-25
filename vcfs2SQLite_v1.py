#!/usr/bin/env python3
import os
import subprocess
import glob
import pandas as pd
import re
import sys
import csv
import sqlite3
from datetime import datetime
import time

def run_command(command, description=None):
    """Execute a shell command and handle errors."""
    if description:
        print(f"Running: {description}")
    try:
        result = subprocess.run(command, shell=True, check=True, 
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}")
        print(f"Error message: {e.stderr}")
        sys.exit(1)

def process_vcf_files(output_dir="."):
    """Process VCF files using bcftools (equivalent to bash script)."""
    print("Step 1: Processing VCF files with bcftools")
    
    # Create file list
    vcf_files = glob.glob("*.vcf.gz")
    if not vcf_files:
        print("No .vcf.gz files found in the current directory.")
        sys.exit(1)
    
    with open("filelist.txt", "w") as f:
        for file in vcf_files:
            f.write(f"{file}\n")
    
    # Index all VCF files
    for vcf_file in vcf_files:
        run_command(f"bcftools index {vcf_file}", f"Indexing {vcf_file}")
    
    # Merge VCF files
    run_command("bcftools merge --file-list filelist.txt -Oz -o merged.vcf.gz", 
              "Merging VCF files")
    
    # Convert to TSV
    run_command("zcat merged.vcf.gz | grep -v '##' > merged.vcf_table.tsv",
              "Converting to TSV")
    
    print("VCF processing completed.")
    return os.path.join(output_dir, "merged.vcf_table.tsv")

def transform_tsv_to_csv(input_file, output_dir="."):
    """Transform TSV file to CSV with R-like manipulations."""
    print(f"Step 2: Transforming {input_file} to CSV")
    
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Step 1: Rename the first column
    df.rename(columns={df.columns[0]: "CHROM"}, inplace=True)
    
    # Step 2: Remove columns 6 to 9 (0-indexed in Python)
    df = df.drop(df.columns[5:9], axis=1)
    
    # Step 3: Replace './.:.:.:.:.' with '0:0:0:0:0' in all cells
    df = df.replace('./.:.:.:.:.',  '0:0:0:0:0')
    
    # Step 4: Keep only the first part of the column names from column 6 to the last
    new_cols = df.columns[:5].tolist()
    for col in df.columns[5:]:
        new_cols.append(col.split('_')[0])
    df.columns = new_cols
    
    # Step 5: Remove 'X' from column names from column 6 to the last
    new_cols = df.columns[:5].tolist()
    for col in df.columns[5:]:
        new_cols.append(col.replace('X', ''))
    df.columns = new_cols
    
    # Write to CSV
    output_file = os.path.join(output_dir, "merged.vcf_table_stage1.1.csv")
    df.to_csv(output_file, index=False)
    print(f"Transformation completed. Output saved to {output_file}")
    return output_file

def count_csv_rows(file_path):
    """Count the number of rows in a CSV file (excluding header)"""
    with open(file_path, 'r') as f:
        return sum(1 for _ in f) - 1  # Subtract 1 for header

def create_database(output_db):
    """Create the SQLite database and table schema"""
    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()
    
    # Create the variants table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS variants (
        file_id TEXT,
        CHROM TEXT,
        POS INTEGER,
        ID TEXT,
        REF TEXT,
        ALT TEXT,
        GT TEXT,
        AD TEXT,
        DP INTEGER,
        GQ INTEGER,
        PL TEXT,
        PRIMARY KEY (file_id, CHROM, POS)
    )
    ''')
    
    # Create indexes
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variants(CHROM, POS)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_id ON variants(ID)')
    
    conn.commit()
    return conn

def transform_and_load_data(input_csv, output_db, batch_size=10000, log_frequency=50000):
    """Transform CSV data and load into SQLite database"""
    start_time = datetime.now()
    
    # Get file size for reporting
    file_size = os.path.getsize(input_csv)
    print(f"Processing file: {input_csv} ({file_size / (1024 * 1024 * 1024):.2f} GB)")
    
    # Count total rows for progress estimation
    print("Estimating total rows (this may take a while for large files)...")
    try:
        # Try to get the row count from wc command (faster)
        result = subprocess.run(['wc', '-l', input_csv], capture_output=True, text=True)
        total_csv_rows = int(result.stdout.split()[0]) - 1  # Subtract header
        print(f"Total rows to process: {total_csv_rows:,}")
    except:
        print("Could not use wc command, falling back to manual count")
        # Fallback to manual counting
        total_csv_rows = None
    
    # Create database and connection
    conn = create_database(output_db)
    cursor = conn.cursor()
    
    # Open CSV file
    with open(input_csv, 'r') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)  # Get header row
        
        # First 5 columns are standard genomic fields
        sample_ids = header[5:]
        
        print(f"Found {len(sample_ids)} sample IDs")
        
        # Process rows in batches
        total_rows = 0
        total_variants = 0
        batch = []
        last_log_time = time.time()
        
        for row_num, row in enumerate(reader, 1):
            basic_data = row[:5]  # CHROM, POS, ID, REF, ALT
            
            # Process each sample
            for i, sample_id in enumerate(sample_ids):
                if i + 5 >= len(row):  # Ensure index is within bounds
                    continue
                    
                cell_value = row[i + 5]
                
                # Skip empty values or values indicating no data
                if not cell_value or cell_value == "0:0:0:0:0" or cell_value == ".:.:.:.:.":
                    continue
                
                # Split the cell value by colon
                parts = cell_value.split(':')
                
                # If we don't have exactly 5 parts, skip or handle appropriately
                if len(parts) != 5:
                    # Possible cases: missing data, different format
                    continue
                
                # Extract individual components
                gt, ad, dp, gq, pl = parts
                
                # Convert types where needed, handle potential parsing errors
                try:
                    dp = int(dp) if dp and dp.isdigit() else None
                    gq = int(gq) if gq and gq.isdigit() else None
                except ValueError:
                    dp = None
                    gq = None
                
                # Create a row for this variant
                try:
                    pos_value = int(basic_data[1]) if basic_data[1].isdigit() else 0
                    batch.append((
                        sample_id,       # file_id
                        basic_data[0],   # CHROM
                        pos_value,       # POS
                        basic_data[2],   # ID
                        basic_data[3],   # REF
                        basic_data[4],   # ALT
                        gt,              # GT
                        ad,              # AD
                        dp,              # DP
                        gq,              # GQ
                        pl               # PL
                    ))
                except Exception as e:
                    print(f"Error processing row {row_num}, sample {sample_id}: {e}")
                    print(f"Row data: {basic_data}")
                    continue
                
                total_variants += 1
                
                # Insert in batches
                if len(batch) >= batch_size:
                    try:
                        cursor.executemany(
                            'INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                            batch
                        )
                        conn.commit()
                        batch = []
                    except sqlite3.Error as e:
                        print(f"SQLite error: {e}")
                        conn.rollback()
            
            total_rows += 1
            
            # Log progress based on time or row count
            current_time = time.time()
            if total_rows % log_frequency == 0 or (current_time - last_log_time) >= 60:
                elapsed = (datetime.now() - start_time).total_seconds()
                last_log_time = current_time
                
                # Calculate progress percentage if we know total rows
                progress_msg = ""
                if total_csv_rows:
                    pct_complete = (total_rows / total_csv_rows) * 100
                    progress_msg = f"{pct_complete:.1f}% complete, "
                
                rows_per_sec = total_rows / elapsed if elapsed > 0 else 0
                
                print(f"Processed {total_rows:,} rows, {total_variants:,} variants "
                      f"({progress_msg}{rows_per_sec:.1f} rows/sec)")
    
    # Insert any remaining rows
    if batch:
        try:
            cursor.executemany(
                'INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                batch
            )
            conn.commit()
        except sqlite3.Error as e:
            print(f"SQLite error when inserting final batch: {e}")
            conn.rollback()
    
    # Log final statistics
    elapsed = (datetime.now() - start_time).total_seconds()
    print(f"\nProcessing complete:")
    print(f"- Processed {total_rows:,} CSV rows")
    print(f"- Generated {total_variants:,} variant records")
    print(f"- Total time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"- Average speed: {total_rows/elapsed:.1f} rows/sec")
    
    # Get database size
    db_size_mb = os.path.getsize(output_db) / (1024 * 1024)
    print(f"- Database size: {db_size_mb:.2f} MB")
    
    # Set pragmas for query optimization
    cursor.execute('PRAGMA journal_mode = WAL')
    cursor.execute('PRAGMA synchronous = NORMAL')
    cursor.execute('PRAGMA cache_size = -10000')
    cursor.execute('ANALYZE')
    
    conn.close()
    return output_db

def main():
    # Default configuration
    output_dir = "."
    output_db = os.path.join(output_dir, "variants.db")
    batch_size = 10000
    log_frequency = 50000
    
    print("Starting VCF processing pipeline...")
    
    # Step 1: Process VCF files (bash script equivalent)
    tsv_file = process_vcf_files(output_dir)
    
    # Step 2: Transform TSV to CSV (R script equivalent)
    csv_file = transform_tsv_to_csv(tsv_file, output_dir)
    
    # Step 3: Transform and load data into SQLite database
    transform_and_load_data(csv_file, output_db, batch_size, log_frequency)
    
    print(f"Pipeline completed. Final database: {output_db}")

if __name__ == "__main__":
    main()