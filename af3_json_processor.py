#!/usr/bin/env python3
"""
AlphaFold 3 JSON Processor

This script processes CSV input files containing AlphaFold 3 job specifications and outputs:
1. JSON files suitable for AlphaFold 3 Server upload (grouped by configurable batch size)
2. CSV file mapping jobs to their output JSON files with token calculations

Features:
- Supports protein, DNA, RNA, ligand, and ion sequences
- PTM (Post-Translational Modification) validation and application via separate columns
- Token cost calculation (max 5000 per job)
- Automatic seed number generation
- Batch processing with configurable JSON file sizes
- Uses CSV reference file (no Excel dependencies needed)

GitHub: https://github.com/JKourelis/af3_json_processor
"""

import pandas as pd
import json
import os
import sys
import random
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import argparse

class AlphaFoldJobProcessor:
    """Main class for processing AlphaFold jobs"""
    
    def __init__(self, reference_file: str = None, jobs_per_json: int = 30):
        """
        Initialize the processor
        
        Args:
            reference_file: Path to Excel file with CCD codes and PTM information (optional)
            jobs_per_json: Number of jobs to include per JSON file (default: 30)
        """
        self.jobs_per_json = jobs_per_json
        
        # Use default reference file if not provided
        if reference_file is None:
            reference_file = "alphafold_reference.csv"
        
        self.reference_data = self._load_reference_data(reference_file)
        
        # Create lookup dictionaries for quick access
        self.ptm_lookup = self._create_ptm_lookup()
        self.ligand_lookup = self._create_ligand_lookup()
        self.ion_lookup = self._create_ion_lookup()
        self.dna_mod_lookup = self._create_dna_mod_lookup()
        self.rna_mod_lookup = self._create_rna_mod_lookup()
        self.glycan_lookup = self._create_glycan_lookup()
        
        # Only accept these standard characters - anything else is invalid
        self.amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
        self.dna_bases = set('ATCG')
        self.rna_bases = set('AUCG')
        
        # Create amino acid mapping for PTM validation
        self.aa_3to1 = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
    
    def _load_reference_data(self, reference_file: str) -> pd.DataFrame:
        """Load reference data from CSV file"""
        try:
            # Get the directory where this script is located
            script_dir = os.path.dirname(os.path.abspath(__file__))
            
            # If reference_file is just a filename, look in script directory
            if not os.path.dirname(reference_file):
                reference_file = os.path.join(script_dir, reference_file)
            
            if os.path.exists(reference_file):
                df = pd.read_csv(reference_file)
                print(f"Loaded reference data from {reference_file}: {len(df)} entries")
                return df
            else:
                raise FileNotFoundError(f"Reference file '{reference_file}' not found in script directory")
        except Exception as e:
            raise Exception(f"Error loading reference file {reference_file}: {e}")
    
    def _create_ptm_lookup(self) -> Dict[str, Dict[str, Any]]:
        """Create PTM lookup dictionary"""
        ptm_data = self.reference_data[self.reference_data['Type'] == 'PTM']
        lookup = {}
        for _, row in ptm_data.iterrows():
            if pd.notna(row['CCD Code']):
                lookup[row['CCD Code']] = {
                    'name': row['Name'],
                    'target_residue': row['Target Residue'],
                    'atom_count': int(row['Heavy Atom Count']) if pd.notna(row['Heavy Atom Count']) else 0
                }
        return lookup
    
    def _create_ligand_lookup(self) -> Dict[str, Dict[str, Any]]:
        """Create ligand lookup dictionary"""
        ligand_data = self.reference_data[self.reference_data['Type'] == 'Ligand']
        lookup = {}
        for _, row in ligand_data.iterrows():
            if pd.notna(row['CCD Code']):
                lookup[row['CCD Code']] = {
                    'name': row['Name'],
                    'atom_count': int(row['Heavy Atom Count']) if pd.notna(row['Heavy Atom Count']) else 1
                }
        return lookup
    
    def _create_ion_lookup(self) -> Dict[str, Dict[str, Any]]:
        """Create ion lookup dictionary"""
        ion_data = self.reference_data[self.reference_data['Type'] == 'Ion']
        lookup = {}
        for _, row in ion_data.iterrows():
            if pd.notna(row['CCD Code']):
                lookup[row['CCD Code']] = {
                    'name': row['Name'],
                    'atom_count': int(row['Heavy Atom Count']) if pd.notna(row['Heavy Atom Count']) else 1
                }
        return lookup

    def _create_dna_mod_lookup(self) -> Dict[str, Dict[str, Any]]:
        """Create DNA modification lookup dictionary"""
        dna_mod_data = self.reference_data[self.reference_data['Type'] == 'DNA Mod']
        lookup = {}
        for _, row in dna_mod_data.iterrows():
            if pd.notna(row['CCD Code']):
                lookup[row['CCD Code']] = {
                    'name': row['Name'],
                    'target_base': row['Target Residue'] if pd.notna(row['Target Residue']) else 'NA',
                    'atom_count': int(row['Heavy Atom Count']) if pd.notna(row['Heavy Atom Count']) else 0
                }
        return lookup

    def _create_rna_mod_lookup(self) -> Dict[str, Dict[str, Any]]:
        """Create RNA modification lookup dictionary"""
        rna_mod_data = self.reference_data[self.reference_data['Type'] == 'RNA Mod']
        lookup = {}
        for _, row in rna_mod_data.iterrows():
            if pd.notna(row['CCD Code']):
                lookup[row['CCD Code']] = {
                    'name': row['Name'],
                    'target_base': row['Target Residue'] if pd.notna(row['Target Residue']) else 'NA',
                    'atom_count': int(row['Heavy Atom Count']) if pd.notna(row['Heavy Atom Count']) else 0
                }
        return lookup

    def _create_glycan_lookup(self) -> Dict[str, Dict[str, Any]]:
        """Create glycan lookup dictionary"""
        glycan_data = self.reference_data[self.reference_data['Type'] == 'Glycan']
        lookup = {}
        for _, row in glycan_data.iterrows():
            if pd.notna(row['CCD Code']):
                lookup[row['CCD Code']] = {
                    'name': row['Name'],
                    'target_residue': row['Target Residue'] if pd.notna(row['Target Residue']) else 'NA',
                    'atom_count': int(row['Heavy Atom Count']) if pd.notna(row['Heavy Atom Count']) else 0
                }
        return lookup

    def _validate_sequence_characters(self, sequence: str, seq_type: str) -> List[str]:
        """
        Validate sequence characters - only accept standard allowed characters
        
        Args:
            sequence: The sequence string to validate
            seq_type: Type of sequence (protein, dna, rna)
            
        Returns:
            List of validation error messages
        """
        errors = []
        sequence = sequence.upper()
        
        if seq_type.lower() == 'protein':
            invalid_chars = []
            
            for i, char in enumerate(sequence, 1):
                if char.isalpha() and char not in self.amino_acids:
                    invalid_chars.append(f"{char}@{i}")
                elif not char.isalpha() and not char.isspace():
                    invalid_chars.append(f"{char}@{i}")
            
            if invalid_chars:
                errors.append(f"Invalid/unknown amino acids: {', '.join(invalid_chars[:10])}" + 
                            ("..." if len(invalid_chars) > 10 else ""))
        
        elif seq_type.lower() == 'dna':
            invalid_chars = []
            
            for i, char in enumerate(sequence, 1):
                if char.isalpha() and char not in self.dna_bases:
                    invalid_chars.append(f"{char}@{i}")
                elif not char.isalpha() and not char.isspace():
                    invalid_chars.append(f"{char}@{i}")
            
            if invalid_chars:
                errors.append(f"Invalid/unknown DNA bases: {', '.join(invalid_chars[:10])}" + 
                            ("..." if len(invalid_chars) > 10 else ""))
        
        elif seq_type.lower() == 'rna':
            invalid_chars = []
            
            for i, char in enumerate(sequence, 1):
                if char.isalpha() and char not in self.rna_bases:
                    invalid_chars.append(f"{char}@{i}")
                elif not char.isalpha() and not char.isspace():
                    invalid_chars.append(f"{char}@{i}")
            
            if invalid_chars:
                errors.append(f"Invalid/unknown RNA bases: {', '.join(invalid_chars[:10])}" + 
                            ("..." if len(invalid_chars) > 10 else ""))
        
        return errors
    
    def _validate_ccd_code(self, code: str, code_type: str) -> List[str]:
        """
        Validate CCD codes for ligands and ions
        
        Args:
            code: The CCD code to validate
            code_type: Type of code ('ligand' or 'ion')
            
        Returns:
            List of validation error messages
        """
        errors = []
        code = code.strip()
        
        if code_type.lower() == 'ligand':
            # Remove CCD_ prefix if present for lookup
            lookup_code = code.replace('CCD_', '') if code.startswith('CCD_') else code
            
            if lookup_code not in self.ligand_lookup:
                # Check if it might be in other categories
                if lookup_code in self.ion_lookup:
                    errors.append(f"'{code}' is an ion, not a ligand. Use seq_type='ion'")
                elif any(lookup_code in ptm_code for ptm_code in self.ptm_lookup.keys()):
                    errors.append(f"'{code}' appears to be a PTM code, not a ligand")
                else:
                    errors.append(f"Unknown ligand CCD code: '{code}'. Check reference file")
        
        elif code_type.lower() == 'ion':
            if code not in self.ion_lookup:
                # Check if it might be in other categories
                if code in self.ligand_lookup:
                    errors.append(f"'{code}' is a ligand, not an ion. Use seq_type='ligand'")
                elif any(code in ptm_code for ptm_code in self.ptm_lookup.keys()):
                    errors.append(f"'{code}' appears to be a PTM code, not an ion")
                else:
                    errors.append(f"Unknown ion CCD code: '{code}'. Check reference file")
        
        return errors
    
    def _parse_ptms(self, ptm_string: str, sequence: str) -> Tuple[List[Dict], List[str]]:
        """
        Parse PTM specifications from separate column
        Format: "CCD_SEP:5;CCD_TPO:12;CCD_PTR:25"

        Args:
            ptm_string: Semicolon-separated PTM specifications
            sequence: The protein sequence to validate against

        Returns:
            tuple: (ptm_list, errors)
        """
        errors = []
        ptms = []

        if pd.isna(ptm_string) or str(ptm_string).strip() == '' or str(ptm_string).lower() == 'nan':
            return ptms, errors

        ptm_string = str(ptm_string).strip()
        
        # Parse PTM entries separated by semicolon
        ptm_entries = [entry.strip() for entry in ptm_string.split(';') if entry.strip()]
        
        for entry in ptm_entries:
            if ':' not in entry:
                errors.append(f"Invalid PTM format: '{entry}'. Use 'CCD_CODE:position'")
                continue
                
            ptm_code, pos_str = entry.split(':', 1)
            ptm_code = ptm_code.strip()
            
            try:
                position = int(pos_str.strip())
                
                # Validate PTM code exists
                if ptm_code not in self.ptm_lookup:
                    # Check if it might be a ligand or ion code
                    lookup_code = ptm_code.replace('CCD_', '') if ptm_code.startswith('CCD_') else ptm_code
                    if lookup_code in self.ligand_lookup:
                        errors.append(f"'{ptm_code}' is a ligand, not a PTM")
                    elif lookup_code in self.ion_lookup:
                        errors.append(f"'{ptm_code}' is an ion, not a PTM")
                    else:
                        errors.append(f"Unknown PTM code: '{ptm_code}'. Check reference file")
                    continue
                
                # Validate position is within sequence
                if position < 1 or position > len(sequence):
                    errors.append(f"PTM position {position} out of range for sequence length {len(sequence)}")
                    continue
                
                # Validate target residue
                target_residue_3letter = self.ptm_lookup[ptm_code]['target_residue']
                actual_residue = sequence[position - 1].upper()  # Convert to 0-based
                
                if target_residue_3letter != 'NA':
                    # Convert 3-letter code to 1-letter
                    target_residue_1letter = self.aa_3to1.get(target_residue_3letter.upper())
                    
                    if target_residue_1letter and actual_residue != target_residue_1letter:
                        errors.append(f"PTM {ptm_code} targets {target_residue_3letter}({target_residue_1letter}) but found {actual_residue} at position {position}")
                        continue
                
                ptms.append({
                    'ptmType': ptm_code,
                    'ptmPosition': position
                })
                
            except ValueError:
                errors.append(f"Invalid PTM position: '{pos_str}' in '{entry}'")

        return ptms, errors

    def _parse_glycans(self, glycan_string: str, sequence: str) -> Tuple[List[Dict], List[str]]:
        """
        Parse glycan specifications from the ptms column
        Format: "NAG:35;NAG(NAG)(BMA):8;BMA:10"

        Args:
            glycan_string: Semicolon-separated glycan specifications
            sequence: The protein sequence to validate against

        Returns:
            tuple: (glycan_list, errors)
        """
        errors = []
        glycans = []

        if pd.isna(glycan_string) or str(glycan_string).strip() == '' or str(glycan_string).lower() == 'nan':
            return glycans, errors

        glycan_string = str(glycan_string).strip()

        # Parse glycan entries separated by semicolon
        glycan_entries = [entry.strip() for entry in glycan_string.split(';') if entry.strip()]

        for entry in glycan_entries:
            if ':' not in entry:
                continue  # Skip if not a glycan entry (might be a PTM)

            # Split at the last colon to handle branched notation
            last_colon = entry.rfind(':')
            residues = entry[:last_colon].strip()
            pos_str = entry[last_colon+1:].strip()

            # Check if this looks like a glycan (contains glycan codes)
            # Extract base glycan code (before any parentheses)
            base_glycan = residues.split('(')[0] if '(' in residues else residues

            # Skip if not a glycan code
            if base_glycan not in self.glycan_lookup:
                continue

            try:
                position = int(pos_str)

                # Validate position is within sequence
                if position < 1 or position > len(sequence):
                    errors.append(f"Glycan position {position} out of range for sequence length {len(sequence)}")
                    continue

                # Validate target residue (N, T, or S)
                actual_residue = sequence[position - 1].upper()
                if actual_residue not in ['N', 'T', 'S']:
                    errors.append(f"Glycan at position {position} requires N/T/S but found {actual_residue}")
                    continue

                # Validate all glycan codes in the branched notation
                # Extract all glycan codes from the notation
                import re
                glycan_codes = re.findall(r'[A-Z]{3}', residues)
                for code in glycan_codes:
                    if code not in self.glycan_lookup:
                        errors.append(f"Unknown glycan code: '{code}' in '{residues}'")
                        break
                else:
                    # All codes valid, add the glycan
                    glycans.append({
                        'residues': residues,
                        'position': position
                    })

            except ValueError:
                continue  # Skip if position is not a number

        return glycans, errors

    def _separate_ptms_and_glycans(self, ptm_string: str, sequence: str) -> Tuple[List[Dict], List[Dict], List[str]]:
        """
        Separate PTMs and glycans from the combined ptm_string

        Returns:
            tuple: (ptms, glycans, errors)
        """
        ptms = []
        glycans = []
        errors = []

        if pd.isna(ptm_string) or str(ptm_string).strip() == '' or str(ptm_string).lower() == 'nan':
            return ptms, glycans, errors

        ptm_string = str(ptm_string).strip()
        entries = [entry.strip() for entry in ptm_string.split(';') if entry.strip()]

        for entry in entries:
            if ':' not in entry:
                errors.append(f"Invalid format: '{entry}'. Use 'CODE:position'")
                continue

            # Check if it's a glycan or PTM
            last_colon = entry.rfind(':')
            code_part = entry[:last_colon].strip()
            pos_str = entry[last_colon+1:].strip()

            # Extract base code
            base_code = code_part.split('(')[0] if '(' in code_part else code_part

            # Check if it's a glycan
            if base_code in self.glycan_lookup:
                # Process as glycan
                glyc, glyc_err = self._parse_glycans(entry, sequence)
                glycans.extend(glyc)
                errors.extend(glyc_err)
            elif base_code in self.ptm_lookup or base_code.startswith('CCD_'):
                # Process as PTM
                ptm, ptm_err = self._parse_ptms(entry, sequence)
                ptms.extend(ptm)
                errors.extend(ptm_err)
            else:
                errors.append(f"Unknown modification code: '{base_code}'")

        return ptms, glycans, errors

    def _parse_nucleic_mods(self, mod_string: str, sequence: str, seq_type: str) -> Tuple[List[Dict], List[str]]:
        """
        Parse DNA/RNA modification specifications from separate column
        Similar format to PTMs: "CCD_5CM:5;CCD_6OG:12"

        Args:
            mod_string: Semicolon-separated modification specifications
            sequence: The DNA/RNA sequence to validate against
            seq_type: Type of sequence ('dna' or 'rna')

        Returns:
            tuple: (mod_list, errors)
        """
        errors = []
        mods = []

        if pd.isna(mod_string) or str(mod_string).strip() == '' or str(mod_string).lower() == 'nan':
            return mods, errors

        mod_string = str(mod_string).strip()

        # Choose the appropriate lookup based on sequence type
        if seq_type.lower() == 'dna':
            mod_lookup = self.dna_mod_lookup
            mod_type = 'DNA modification'
            valid_bases = self.dna_bases
        elif seq_type.lower() == 'rna':
            mod_lookup = self.rna_mod_lookup
            mod_type = 'RNA modification'
            valid_bases = self.rna_bases
        else:
            return mods, errors

        # Parse modification entries separated by semicolon
        mod_entries = [entry.strip() for entry in mod_string.split(';') if entry.strip()]

        for entry in mod_entries:
            if ':' not in entry:
                errors.append(f"Invalid {mod_type} format: '{entry}'. Use 'CCD_CODE:position'")
                continue

            mod_code, pos_str = entry.split(':', 1)
            mod_code = mod_code.strip()

            try:
                position = int(pos_str.strip())

                # Validate modification code exists
                if mod_code not in mod_lookup:
                    errors.append(f"Unknown {mod_type} code: '{mod_code}'. Check reference file")
                    continue

                # Validate position is within sequence
                if position < 1 or position > len(sequence):
                    errors.append(f"{mod_type} position {position} out of range for sequence length {len(sequence)}")
                    continue

                # Validate target base if specified
                target_base = mod_lookup[mod_code]['target_base']
                actual_base = sequence[position - 1].upper()  # Convert to 0-based

                if target_base != 'NA':
                    # For DNA mods, target might be like 'DC' for cytosine
                    # Extract just the base letter
                    if len(target_base) > 1 and target_base.startswith('D'):
                        target_base = target_base[1]  # Remove 'D' prefix

                    if actual_base != target_base:
                        errors.append(f"{mod_type} {mod_code} targets {target_base} but found {actual_base} at position {position}")
                        continue

                mods.append({
                    'modificationType': mod_code,  # Correct key name for DNA/RNA mods
                    'basePosition': position
                })

            except ValueError:
                errors.append(f"Invalid {mod_type} position: '{pos_str}' in '{entry}'")

        return mods, errors
    
    def _sanitize_jobname(self, jobname: str) -> Tuple[str, List[str]]:
        """
        Sanitize jobname according to AlphaFold rules
        Allowed characters: letters, numbers, spaces, dashes, underscores, colons
        When downloading, everything converts to lowercase with dashes→underscores

        Args:
            jobname: Raw jobname string

        Returns:
            tuple: (sanitized_jobname, errors)
        """
        errors = []

        if pd.isna(jobname) or str(jobname).strip() == '':
            errors.append("Missing jobname")
            return '', errors

        original_jobname = str(jobname)

        # Convert to lowercase
        sanitized = original_jobname.lower()

        # Replace dashes with underscores
        sanitized = sanitized.replace('-', '_')

        # Keep only allowed characters: a-z, 0-9, space, underscore, colon
        allowed_pattern = re.compile(r'[^a-z0-9 _:]')
        sanitized = allowed_pattern.sub('', sanitized)

        # Truncate to 128 characters
        if len(sanitized) > 128:
            sanitized = sanitized[:128]
            errors.append(f"Jobname truncated from {len(original_jobname)} to 128 characters")

        # Check if sanitized name is empty or whitespace only
        if not sanitized.strip():
            errors.append(f"Jobname '{original_jobname}' became empty after sanitization")
            return '', errors

        # Report if jobname was changed
        if sanitized != original_jobname:
            # Don't report as error, just a note - jobname sanitization is normal
            pass

        return sanitized, errors
    
    def _generate_seed(self, seed_value: str) -> int:
        """Generate a random seed if not provided - uses 32-bit signed integer range (AlphaFold Server requirement)"""
        if pd.isna(seed_value) or str(seed_value).strip() == '':
            return random.randint(0, 2147483647)  # 2^31 - 1 (32-bit signed max)
        try:
            seed = int(seed_value)
            # Validate seed is within AlphaFold Server's allowed range (32-bit signed integer)
            if seed < 0 or seed > 2147483647:
                return random.randint(0, 2147483647)
            return seed
        except:
            return random.randint(0, 2147483647)
    
    def _calculate_tokens(self, sequence: str, seq_type: str, ptms: List[Dict] = None) -> Tuple[int, List[str]]:
        """
        Calculate token count for a sequence
        
        Returns:
            tuple: (token_count, errors)
        """
        errors = []
        tokens = 0
        
        if seq_type.lower() == 'protein':
            # Count amino acids
            for char in sequence.upper():
                if char in self.amino_acids:
                    tokens += 1
                elif char.isalpha():  # Non-standard amino acid
                    tokens += 1  # Still count as 1 token
            
            # Add PTM tokens
            if ptms:
                for ptm in ptms:
                    ptm_code = ptm['ptmType']
                    if ptm_code in self.ptm_lookup:
                        tokens += self.ptm_lookup[ptm_code]['atom_count']
                    else:
                        errors.append(f"Unknown PTM for token calculation: {ptm_code}")
        
        elif seq_type.lower() == 'dna':
            for char in sequence.upper():
                if char in self.dna_bases:
                    tokens += 1
                elif char.isalpha():
                    tokens += 1  # Count non-standard bases as 1 token

            # Add DNA modification tokens
            if ptms:  # Using ptms parameter which contains modifications
                for mod in ptms:
                    mod_code = mod.get('modificationType', mod.get('ptmType', ''))
                    if mod_code in self.dna_mod_lookup:
                        tokens += self.dna_mod_lookup[mod_code]['atom_count']
                    else:
                        errors.append(f"Unknown DNA mod for token calculation: {mod_code}")

        elif seq_type.lower() == 'rna':
            for char in sequence.upper():
                if char in self.rna_bases:
                    tokens += 1
                elif char.isalpha():
                    tokens += 1  # Count non-standard bases as 1 token

            # Add RNA modification tokens
            if ptms:  # Using ptms parameter which contains modifications
                for mod in ptms:
                    mod_code = mod.get('modificationType', mod.get('ptmType', ''))
                    if mod_code in self.rna_mod_lookup:
                        tokens += self.rna_mod_lookup[mod_code]['atom_count']
                    else:
                        errors.append(f"Unknown RNA mod for token calculation: {mod_code}")
        
        elif seq_type.lower() == 'ligand':
            # For ligands, use the CCD code to look up atom count
            ligand_code = sequence.strip()
            lookup_code = ligand_code.replace('CCD_', '') if ligand_code.startswith('CCD_') else ligand_code
            if lookup_code in self.ligand_lookup:
                tokens += self.ligand_lookup[lookup_code]['atom_count']
            else:
                errors.append(f"Unknown ligand CCD code: {ligand_code}")
                tokens += 1  # Default to 1 token
        
        elif seq_type.lower() == 'ion':
            # For ions, use the CCD code to look up atom count
            ion_code = sequence.strip()
            if ion_code in self.ion_lookup:
                tokens += self.ion_lookup[ion_code]['atom_count']
            else:
                errors.append(f"Unknown ion CCD code: {ion_code}")
                tokens += 1  # Default to 1 token
        
        return tokens, errors
    
    def _process_sequence(self, name: str, seq_type: str, copies: int, sequence: str, ptm_string: str = '') -> Tuple[Dict, int, List[str]]:
        """
        Process a single sequence and return JSON structure
        
        Returns:
            tuple: (sequence_dict, total_tokens, errors)
        """
        errors = []
        total_tokens = 0
        
        if pd.isna(sequence) or str(sequence).strip() == '':
            return None, 0, []
        
        sequence = str(sequence).strip()
        copies = int(copies) if pd.notna(copies) and str(copies).strip() != '' else 1
        
        if seq_type.lower() == 'protein':
            # Validate sequence characters
            char_errors = self._validate_sequence_characters(sequence, seq_type)
            errors.extend(char_errors)

            # Parse PTMs and glycans from the ptm_string column
            ptms, glycans, mod_errors = self._separate_ptms_and_glycans(ptm_string, sequence)
            errors.extend(mod_errors)

            # Calculate tokens (include both PTMs and glycans)
            seq_tokens, token_errors = self._calculate_tokens(sequence, seq_type, ptms)
            errors.extend(token_errors)

            # Add glycan tokens
            for glycan in glycans:
                import re
                glycan_codes = re.findall(r'[A-Z]{3}', glycan['residues'])
                for code in glycan_codes:
                    if code in self.glycan_lookup:
                        seq_tokens += self.glycan_lookup[code]['atom_count']

            total_tokens = seq_tokens * copies

            seq_dict = {
                "proteinChain": {
                    "sequence": sequence,
                    "count": copies,
                    "useStructureTemplate": True
                }
            }

            if glycans:
                seq_dict["proteinChain"]["glycans"] = glycans
            if ptms:
                seq_dict["proteinChain"]["modifications"] = ptms
        
        elif seq_type.lower() == 'dna':
            # Validate sequence characters
            char_errors = self._validate_sequence_characters(sequence, seq_type)
            errors.extend(char_errors)

            # Parse DNA modifications from separate column
            mods, mod_errors = self._parse_nucleic_mods(ptm_string, sequence, 'dna')
            errors.extend(mod_errors)

            seq_tokens, token_errors = self._calculate_tokens(sequence, seq_type, mods)
            errors.extend(token_errors)
            total_tokens = seq_tokens * copies

            seq_dict = {
                "dnaSequence": {
                    "sequence": sequence,
                    "count": copies
                }
            }

            if mods:
                seq_dict["dnaSequence"]["modifications"] = mods
        
        elif seq_type.lower() == 'rna':
            # Validate sequence characters
            char_errors = self._validate_sequence_characters(sequence, seq_type)
            errors.extend(char_errors)

            # Parse RNA modifications from separate column
            mods, mod_errors = self._parse_nucleic_mods(ptm_string, sequence, 'rna')
            errors.extend(mod_errors)

            seq_tokens, token_errors = self._calculate_tokens(sequence, seq_type, mods)
            errors.extend(token_errors)
            total_tokens = seq_tokens * copies

            seq_dict = {
                "rnaSequence": {
                    "sequence": sequence,
                    "count": copies
                }
            }

            if mods:
                seq_dict["rnaSequence"]["modifications"] = mods
        
        elif seq_type.lower() == 'ligand':
            ligand_code = sequence.strip()
            
            # Validate ligand CCD code
            ccd_errors = self._validate_ccd_code(ligand_code, 'ligand')
            errors.extend(ccd_errors)
            
            seq_tokens, token_errors = self._calculate_tokens(sequence, seq_type)
            errors.extend(token_errors)
            total_tokens = seq_tokens * copies
            
            seq_dict = {
                "ligand": {
                    "ligand": f"CCD_{ligand_code}" if not ligand_code.startswith('CCD_') else ligand_code,
                    "count": copies
                }
            }
        
        elif seq_type.lower() == 'ion':
            ion_code = sequence.strip()
            
            # Validate ion CCD code
            ccd_errors = self._validate_ccd_code(ion_code, 'ion')
            errors.extend(ccd_errors)
            
            seq_tokens, token_errors = self._calculate_tokens(sequence, seq_type)
            errors.extend(token_errors)
            total_tokens = seq_tokens * copies
            
            seq_dict = {
                "ion": {
                    "ion": ion_code,
                    "count": copies
                }
            }
        
        else:
            errors.append(f"Unknown sequence type: {seq_type}")
            return None, 0, errors
        
        return seq_dict, total_tokens, errors
    
    def _process_job(self, row: pd.Series) -> Tuple[Dict, int, List[str]]:
        """Process a single job row from CSV"""
        errors = []
        total_tokens = 0
        sequences = []
        
        # Sanitize jobname first
        raw_jobname = row.get('jobname', '')
        sanitized_jobname, jobname_errors = self._sanitize_jobname(raw_jobname)
        errors.extend(jobname_errors)
        
        # If jobname is empty after sanitization, this job will be excluded
        if not sanitized_jobname:
            return None, 0, errors
        
        # Generate seed
        seed = self._generate_seed(row.get('seed', ''))
        
        # Process up to 10 sequences
        for i in range(1, 11):
            seq_name = row.get(f'seq{i}_name')
            seq_type = row.get(f'seq{i}_type')
            seq_copies = row.get(f'seq{i}_copies')
            seq_sequence = row.get(f'seq{i}')
            seq_ptms = row.get(f'seq{i}_ptms', '')  # PTM column
            
            # Skip if no sequence name or sequence
            if pd.isna(seq_name) or pd.isna(seq_sequence):
                continue
            
            seq_dict, seq_tokens, seq_errors = self._process_sequence(
                str(seq_name), str(seq_type), seq_copies, str(seq_sequence), str(seq_ptms)
            )
            
            if seq_dict:
                sequences.append(seq_dict)
                total_tokens += seq_tokens
            
            errors.extend(seq_errors)
        
        # Create job structure with sanitized jobname
        job = {
            "name": sanitized_jobname,  # Use sanitized jobname
            "modelSeeds": [seed],
            "sequences": sequences,
            "dialect": "alphafoldserver",
            "version": 1
        }
        
        return job, total_tokens, errors
    
    def process_csv(self, input_csv: str, output_dir: str = None) -> str:
        """
        Process the input CSV file and generate JSON files and output CSV

        Args:
            input_csv: Path to input CSV file
            output_dir: Output directory (defaults to same directory as input CSV)

        Returns:
            Path to output CSV file
        """
        if output_dir is None:
            output_dir = os.path.dirname(input_csv)
            # If dirname returns empty string (file in current dir), use current directory
            if not output_dir:
                output_dir = '.'

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Read input CSV
        try:
            df = pd.read_csv(input_csv)
            print(f"Processing {len(df)} jobs from {input_csv}")
        except Exception as e:
            raise Exception(f"Error reading CSV file {input_csv}: {e}")
        
        # Process each job
        processed_jobs = []
        output_rows = []
        excluded_jobs = 0
        
        for idx, row in df.iterrows():
            job, tokens, errors = self._process_job(row)
            
            # Prepare output row (copy original data)
            output_row = row.to_dict()
            
            # Sanitize jobname for output CSV
            raw_jobname = row.get('jobname', '')
            sanitized_jobname, jobname_errors = self._sanitize_jobname(raw_jobname)
            
            # Update jobname in output CSV to sanitized version
            output_row['jobname'] = sanitized_jobname
            output_row['token_count'] = tokens
            output_row['status'] = 'OK'
            output_row['json_file'] = ''
            output_row['errors'] = ''
            
            # Check if job was excluded due to invalid jobname
            if job is None:
                excluded_jobs += 1
                output_row['status'] = 'EXCLUDED'
                output_row['errors'] = '; '.join(errors)
                output_rows.append(output_row)
                continue
            
            # Check for other errors
            if errors:
                output_row['status'] = 'ERROR'
                output_row['errors'] = '; '.join(errors)
            elif tokens > 5000:
                output_row['status'] = 'ERROR'
                output_row['errors'] = f'Token count ({tokens}) exceeds maximum (5000)'
            else:
                processed_jobs.append((job, idx))
            
            output_rows.append(output_row)
        
        # Group jobs into JSON files
        base_name = Path(input_csv).stem
        json_files_created = []
        
        for i in range(0, len(processed_jobs), self.jobs_per_json):
            batch = processed_jobs[i:i + self.jobs_per_json]
            json_filename = f"{base_name}_{i//self.jobs_per_json + 1:04d}.json"
            json_path = os.path.join(output_dir, json_filename)
            
            # Create JSON structure
            json_data = [job for job, _ in batch]
            
            # Save JSON file
            with open(json_path, 'w') as f:
                json.dump(json_data, f, indent=2)
            
            json_files_created.append(json_filename)
            
            # Update output rows with JSON filename
            for job, original_idx in batch:
                output_rows[original_idx]['json_file'] = json_filename
            
            print(f"Created {json_filename} with {len(batch)} jobs")
        
        # Save output CSV
        output_csv = os.path.join(output_dir, f"{base_name}_processed.csv")
        output_df = pd.DataFrame(output_rows)

        # Remove columns that are entirely empty (all NaN or empty strings)
        output_df = output_df.dropna(axis=1, how='all')
        # Also remove columns that only contain empty strings
        for col in output_df.columns:
            if output_df[col].dtype == 'object':
                if output_df[col].apply(lambda x: str(x).strip() == '' or pd.isna(x)).all():
                    output_df = output_df.drop(columns=[col])

        output_df.to_csv(output_csv, index=False)
        
        print(f"\nProcessing complete:")
        print(f"- Created {len(json_files_created)} JSON files")
        print(f"- Successfully processed {len(processed_jobs)} jobs")
        print(f"- Excluded {excluded_jobs} jobs (invalid jobnames)")
        print(f"- Errors in {len(output_rows) - len(processed_jobs) - excluded_jobs} jobs")
        print(f"- Output CSV: {output_csv}")
        
        return output_csv


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(description='Process AlphaFold jobs from CSV to JSON')
    parser.add_argument('input_csv', help='Input CSV file with job specifications')
    parser.add_argument('--reference-file', default='alphafold_reference.csv',
                       help='CSV file with CCD codes and PTM information (default: alphafold_reference.csv)')
    parser.add_argument('--jobs-per-json', type=int, default=30, 
                       help='Number of jobs per JSON file (default: 30)')
    parser.add_argument('--output-dir', help='Output directory (default: same as input CSV)')
    
    args = parser.parse_args()
    
    try:
        processor = AlphaFoldJobProcessor(args.reference_file, args.jobs_per_json)
        output_csv = processor.process_csv(args.input_csv, args.output_dir)
        print(f"\nSuccess! Output saved to: {output_csv}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
