# MCAAT Benchmarking Tool
A comprehensive benchmarking tool designed to evaluate the performance of the [Metagenomic CRISPR Array Analysis Tool (MCAAT)](https://github.com/RNABioInfo/mcaat) using artificial sequencing data generated from known CRISPR-containing genomes.

## Overview
This tool automates the process of:
1. Downloading genomes with known CRISPR arrays from NCBI
2. Generating artificial sequencing reads using InSilicoSeq
3. Running MCAAT on the synthetic data
4. Comparing results against known CRISPR locations for accuracy assessment

## System Requirements
- **Operating System**: Linux (tested on Ubuntu 20.04+)
- **PostgreSQL**: Version 12 or higher
- **Memory**: Minimum 4GB RAM (8GB+ recommended)
- **Storage**: At least 10GB free space for genome downloads and read generation

## Setup Instructions
### 1. Install Required Dependencies
#### PostgreSQL
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install postgresql postgresql-contrib

# Start PostgreSQL service
sudo systemctl start postgresql
sudo systemctl enable postgresql
```

#### Additional Tools
```bash
# Install unzip
sudo apt-get install unzip

# Install NCBI Datasets CLI
# Download from: https://github.com/ncbi/datasets
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets
sudo mv datasets /usr/local/bin/

# Install InSilicoSeq
pip install InSilicoSeq
# or
conda install -c bioconda insilicoseq
```

### 2. Database Setup
#### Create Database and User
```bash
# Switch to postgres user
sudo -u postgres psql

# In PostgreSQL shell:
CREATE DATABASE crispr_cas_db;
CREATE USER crispr_user WITH PASSWORD 'your_password';
GRANT ALL PRIVILEGES ON DATABASE crispr_cas_db TO crispr_user;
\q
```

#### Download and Import CRISPR-Cas Database
```bash
# Download the SQL dump from CRISPR-Cas database
wget https://crisprcas.i2bc.paris-saclay.fr/dumps/latest_dump.sql.gz
gunzip latest_dump.sql.gz

# Import the database
psql -U crispr_user -d crispr_cas_db -f latest_dump.sql
```

### 3. Generate Genome-CRISPR Mapping File
Connect to your PostgreSQL database and run the following SQL command:

```bash
psql -U crispr_user -d crispr_cas_db
```

```sql
COPY (
    SELECT s.id, s.length, s.description, s.genbank, COUNT(*) as crispr_count, 
           string_agg('(' || c.start || '_' || c.length || ')', ';') as crispr_region
    FROM crisprlocus c
    JOIN (
        SELECT sequence.id, sequence.length, sequence.description, strain.genbank
        FROM sequence
        JOIN strain ON strain.id = sequence.strain
    ) as s ON c.sequence = s.id
    GROUP BY s.id, s.length, s.description, s.genbank
    HAVING COUNT(*) > 0
    ORDER BY crispr_count DESC
) TO '/tmp/genome_crispr_combination.csv'
WITH CSV HEADER;
```

Move the generated file to your working directory:
```bash
mv /tmp/genome_crispr_combination.csv ./genome_crispr_combination.csv
```

### 4. Compile MCAAT
Follow the instructions from the [MCAAT repository](https://github.com/RNABioInfo/mcaat) to compile the tool.
You'll need the binary executable for the benchmarking process.

## Usage
### Basic Syntax
```bash
./run_benchmark.sh <mcaat_binary> <genome_crispr_combination.csv> --ids <id1> [id2 ...] [OPTIONS]
```

### Required Arguments
- `mcaat_binary`: Path to the compiled MCAAT executable
- `genome_crispr_combination.csv`: The CSV file generated during setup
- `--ids`: Space-separated list of genome IDs from the CSV file, or a number for random sampling

### Optional Arguments
- `--ram <size>`: Memory limit (e.g., 4G, 500M). Default: 80% of system RAM
- `--threads <count>`: Number of threads to use. Default: 2
- `--output_folder <path>`: Output directory for results. Default: "benchmark"
- `--quick`: Enable quick mode (uses reduced genome regions around CRISPR arrays)

### Example Usage
#### Benchmark Specific Genomes
```bash
# Benchmark specific genome IDs
./run_benchmark.sh ./mcaat genome_crispr_combination.csv --ids 12345 67890 --threads 4 --ram 8G

# Quick benchmark with reduced genome size
./run_benchmark.sh ./mcaat genome_crispr_combination.csv --ids 12345 --quick --output_folder quick_test
```

#### Random Sampling
```bash
# Benchmark 5 random genomes
./run_benchmark.sh ./mcaat genome_crispr_combination.csv --ids 5 --threads 8

# Quick benchmark of 10 random genomes
./run_benchmark.sh ./mcaat genome_crispr_combination.csv --ids 10 --quick --ram 16G
```

#### Production Benchmarking
```bash
# Comprehensive benchmark with high resources
./run_benchmark.sh ./mcaat genome_crispr_combination.csv --ids 50 --threads 16 --ram 32G --output_folder production_benchmark
```

### Output Structure
The tool creates the following directory structure:
```
benchmark/
├── genomes/           # Downloaded genome FASTA files
├── reads/             # Artificially generated sequencing reads
├── crispr_sequences/  # Expected CRISPR sequences extracted from genomes
└── results/           # MCAAT benchmark results and performance metrics
```

### Performance Considerations
- **Quick Mode**: Reduces genome size to regions around CRISPR arrays, significantly faster but may miss some edge cases
- **Memory Usage**: Adjust `--ram` based on your system capabilities
- **Thread Count**: Optimal thread count is typically equal to your CPU cores
- **Storage**: Each genome with reads can require 1-5GB of storage

## Troubleshooting
### Common Issues
1. **PostgreSQL Connection Error**
   - Ensure PostgreSQL is running: `sudo systemctl status postgresql`
   - Check database credentials and permissions
2. **NCBI Datasets Download Failures**
   - Verify internet connection
   - Check if the GenBank accession exists
   - Some genomes may be restricted or moved
3. **InSilicoSeq Errors**
   - Ensure sufficient disk space
   - Check that the genome FASTA file is valid
4. **Memory Issues**
   - Reduce `--ram` parameter
   - Use `--quick` mode for large genomes
   - Close other applications to free memory

### Getting Help
- Check the [MCAAT repository](https://github.com/RNABioInfo/mcaat) for tool-specific issues
- Ensure all dependencies are properly installed and accessible in your PATH
- Verify that your PostgreSQL database contains the expected tables and data

## Citation
If you use this benchmarking tool in your research, please cite the original MCAAT paper and mention this benchmarking framework.

## License
This tool is provided under the same license as MCAAT. Please refer to the main repository
