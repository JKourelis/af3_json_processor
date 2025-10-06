# AlphaFold 3 JSON Processor

A lightweight Python tool for processing CSV files containing AlphaFold 3 job specifications and converting them to JSON format suitable for AlphaFold 3 Server upload.

## Features

- **Multi-sequence support**: Handle up to 10 sequences per job (protein, DNA, RNA, ligand, ion)
- **PTM validation**: Post-translational modification support with residue-specific validation
- **Token calculation**: Automatic calculation of job token costs with 5000 token limit
- **Batch processing**: Group jobs into configurable JSON file sizes (default: 30 jobs per file)
- **Error reporting**: Comprehensive validation and error reporting
- **Seed management**: Automatic random seed generation for reproducibility
- **Lightweight**: Only requires pandas - no Excel dependencies!

## Installation

1. Clone the repository:
```
git clone https://github.com/JKourelis/af3_json_processor
cd af3_json_processor
```

2. Install required dependencies:
```
pip install -r requirements.txt
```

**Note:** Only pandas is required - no Excel dependencies needed!

## Usage

### Command Line

```
python af3_json_processor.py input.csv [options]
```

**Arguments:**
- `input.csv`: CSV file with job specifications

**Options:**
- `--reference-file FILE`: Reference CSV file (default: alphafold_reference.csv)
- `--jobs-per-json N`: Number of jobs per JSON file (default: 30)
- `--output-dir DIR`: Output directory (default: same as input CSV)

### Python Module

```python
from af3_json_processor import AlphaFoldJobProcessor

# Initialize processor (uses alphafold_reference.csv by default)
processor = AlphaFoldJobProcessor()

# Process CSV file
output_csv = processor.process_csv('input_jobs.csv')
```

## Input CSV Format

The input CSV must contain the following columns:

**Required columns:**
- `jobname`: Unique identifier for the job (will be sanitized according to AlphaFold rules)
- `seed`: Random seed (optional, will generate if empty - range: 0 to 2,147,483,647)

**Sequence columns (up to 10 sequences per job):**
- `seq1_name`, `seq1_type`, `seq1_copies`, `seq1`, `seq1_ptms`
- `seq2_name`, `seq2_type`, `seq2_copies`, `seq2`, `seq2_ptms`
- ... up to `seq10_name`, `seq10_type`, `seq10_copies`, `seq10`, `seq10_ptms`

**Note:** The `seq*_ptms` columns are used for:
- Protein post-translational modifications (PTMs) and glycans
- DNA modifications
- RNA modifications

**Sequence types:**
- `protein`: Amino acid sequences (supports PTMs via separate column)
- `dna`: DNA sequences
- `rna`: RNA sequences
- `ligand`: Ligand CCD codes (ATP, ADP, etc.)
- `ion`: Ion CCD codes (MG, ZN, CA, etc.)

### Modifications Specification Format

#### Protein PTMs and Glycans
For protein sequences, PTMs and glycans are specified in the `seq*_ptms` columns:

**PTMs:**
```
CCD_SEP:8;CCD_TPO:12;CCD_PTR:25
```

**Example:**
- `CCD_SEP:8` = Phosphoserine at position 8 (position 8 must be Serine)
- `CCD_MLY:15` = Methyllysine at position 15 (position 15 must be Lysine)

**Glycans:**
```
NAG:2;NAG(NAG)(BMA):8;BMA:10
```

**Example:**
- `NAG:2` = N-Acetylglucosamine at position 2 (position 2 must be N, T, or S)
- `NAG(NAG)(BMA):8` = Branched glycan at position 8 (position 8 must be N, T, or S)

**Combined PTMs and Glycans:**
```
NAG:2;CCD_SEP:8;BMA:10;CCD_TPO:12
```

#### DNA Modifications
For DNA sequences, modifications use the same `seq*_ptms` columns:

```
CCD_5CM:5;CCD_6OG:12
```

**Example:**
- `CCD_5CM:5` = 5-Methylcytosine at position 5 (position 5 must be Cytosine)
- `CCD_6OG:12` = 8-Oxoguanine at position 12 (position 12 must be Guanine)

#### RNA Modifications
For RNA sequences, modifications also use the `seq*_ptms` columns:

```
CCD_PSU:10;CCD_5MC:15
```

**Example:**
- `CCD_PSU:10` = Pseudouridine at position 10 (position 10 must be Uridine)
- `CCD_5MC:15` = 5-Methylcytosine at position 15 (position 15 must be Cytosine)

**Format Rules for All Modifications:**
- Modification codes must match those in the reference CSV file
- Positions are 1-based
- Multiple modifications separated by semicolons
- The script validates that the target base/residue at each position matches the modification requirement

### Jobname Sanitization Rules

Jobnames are automatically sanitized according to AlphaFold requirements:

**Allowed characters**: letters, numbers, spaces, dashes, underscores, colons
**Automatic conversions (when downloaded)**:
- Convert to lowercase: `MyJob` → `myjob`
- Replace dashes with underscores: `my-job` → `my_job`
- Remove invalid characters: `job@test#` → `jobtest`
- Truncate to 128 characters maximum

**Jobs will be excluded if**:
- Missing jobname
- Jobname becomes empty after sanitization (e.g., `@#$%`)

## Reference CSV File Format

The reference CSV file (`alphafold_reference.csv`) must contain columns:
- `Type`: Entry type (PTM, Ligand, Ion, Glycan, DNA Mod, RNA Mod)
- `CCD Code`: Chemical Component Dictionary code
- `Name`: Descriptive name (use "a" and "b" instead of α and β symbols)
- `Target Residue`: For modifications, the target amino acid/base (e.g., SER, THR, DC, U)
- `Molecular Formula`: Chemical formula (informational)
- `All Atom Count`: Total number of atoms (informational)
- `Heavy Atom Count`: Number of heavy atoms used for token calculation

**File Location:**
- Place `alphafold_reference.csv` in the same directory as `af3_json_processor.py`
- Or specify custom location with `--reference-file`

## Output Files

### JSON Files
- Named as `{input_filename}_0001.json`, `{input_filename}_0002.json`, etc.
- Each contains up to 30 jobs (configurable)
- Ready for AlphaFold server upload

### Processed CSV
- Named as `{input_filename}_processed.csv`
- Contains original data plus:
  - `token_count`: Calculated token cost (using heavy atoms for ligands/ions/modifications)
  - `status`: Processing status (OK/ERROR/EXCLUDED)
  - `json_file`: Assigned JSON file name
  - `errors`: Error messages if any
- Empty columns are automatically removed for cleaner output

## Token Calculation

- **Protein**: 1 token per amino acid + PTM/glycan heavy atom counts
- **DNA**: 1 token per nucleotide + DNA modification heavy atom counts
- **RNA**: 1 token per nucleotide + RNA modification heavy atom counts
- **Ligand/Ion**: Heavy atom count from reference file
- **Maximum**: 5000 tokens per job
- **Note**: Heavy atoms (non-hydrogen) are used for more accurate token estimation

## Error Handling

The processor validates:
- **Sequence characters**: Only accepts standard allowed characters for each sequence type
- **PTM/Glycan validation**: PTM and glycan codes exist in reference data, target residues match sequence (N/T/S for glycans), positions are within bounds
- **DNA/RNA modification validation**: Modification codes exist, target bases match sequence, positions are within bounds
- **CCD code validation**: Ligand/ion codes exist and are used with correct sequence types
- **Cross-category detection**: Identifies when ligand codes are used as ions, PTM codes as ligands, etc.
- **Token limits**: Jobs do not exceed 5000 tokens
- **Sequence types**: Valid sequence types are specified

### Validation Details

**Protein Sequences:**
- **Allowed**: Only the 20 standard amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
- **Rejected**: Any other characters including X, Z, B, numbers, symbols

**DNA Sequences:**
- **Allowed**: Only A, T, C, G
- **Rejected**: Any other characters including ambiguous IUPAC codes (N, R, Y, etc.)

**RNA Sequences:**
- **Allowed**: Only A, U, C, G
- **Rejected**: Any other characters including ambiguous IUPAC codes

**Ligands/Ions/PTMs/Glycans:**
- **Allowed**: Only CCD codes defined in the reference CSV file
- **Cross-validation**: Prevents using ligand codes as ions and vice versa

## Example

### Simple Usage
```bash
# Basic usage (looks for alphafold_reference.csv automatically)
python af3_json_processor.py my_jobs.csv

# With custom options
python af3_json_processor.py my_jobs.csv --jobs-per-json 50 --output-dir results
```

### Input CSV
```csv
jobname,seq1_name,seq1_type,seq1_copies,seq1,seq1_ptms,seq2_name,seq2_type,seq2_copies,seq2,seq2_ptms,seed
ValidComplex,ProteinA,protein,1,MKTAYIAKQRQISFVKS,"NAG:2;CCD_SEP:8",ATP_ligand,ligand,2,ATP,,12345678
ErrorJob,BadProtein,protein,1,MKTXAYIAK,,WrongIon,ion,1,UNKNOWN_ION,,87654321
```

**Note:** The second job will generate errors:
- Invalid/unknown amino acid X at position 4  
- Unknown ion code "UNKNOWN_ION"

### Generated JSON
```json
[{
  "name": "ValidComplex",
  "modelSeeds": [12345678],
  "sequences": [
    {
      "proteinChain": {
        "sequence": "MKTAYIAKQRQISFVKS",
        "count": 1,
        "useStructureTemplate": true,
        "glycans": [
          {"residues": "NAG", "position": 2}
        ],
        "modifications": [
          {"ptmType": "CCD_SEP", "ptmPosition": 8}
        ]
      }
    },
    {
      "ligand": {
        "ligand": "CCD_ATP",
        "count": 2
      }
    }
  ],
  "dialect": "alphafoldserver",
  "version": 1
}]
```

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Support

For issues and questions:
1. Check existing issues on GitHub
2. Create a new issue with detailed description
3. Include sample input files if possible