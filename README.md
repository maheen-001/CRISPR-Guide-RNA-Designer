# CRISPR gRNA Designer

A Python/Flask web application for identifying and ranking candidate CRISPR guide RNAs (gRNAs) from a DNA sequence.

The project combines bioinformatics sequence processing with a web-based interface to identify potential guide sequences, evaluate PAM compatibility, calculate basic sequence characteristics, and rank candidate guides using a rule-based scoring heuristic.

> **Project status:** Completed educational/portfolio project

## Overview

CRISPR guide RNA design involves identifying DNA sequences that are compatible with a particular CRISPR-associated enzyme while considering sequence characteristics that may affect guide quality.

This project provides a simple computational workflow for:

1. Accepting a DNA sequence through text input or FASTA upload
2. Validating and normalizing the sequence
3. Searching for candidate gRNAs on both DNA strands
4. Identifying compatible PAM sequences
5. Classifying canonical and non-canonical PAMs
6. Calculating GC content
7. Checking for simple repetitive sequence patterns
8. Estimating potential off-target matches
9. Ranking candidate gRNAs using a rule-based scoring system
10. Exporting results as CSV or PDF

The application is intended as a learning and portfolio project demonstrating the use of programming and data-processing techniques in a biological context.

## Technologies

* **Python**
* **Flask**
* **HTML/CSS**
* **Git / GitHub**
* **ReportLab**
* **CSV processing**

## Supported CRISPR Enzymes

The current implementation supports several Cas enzymes with different guide lengths and PAM requirements:

| Enzyme | Guide Length | PAM Pattern | PAM Location |
| ------ | -----------: | ----------- | ------------ |
| SpCas9 |        20 bp | NGG, NAG    | Downstream   |
| SaCas9 |        21 bp | NNGRRT      | Downstream   |
| Cas12a |        23 bp | TTTV        | Upstream     |

PAM matching supports IUPAC nucleotide codes such as `N`, `R`, and `V`.

## How It Works

### 1. Sequence Input

The application accepts DNA sequences through:

* Direct text input
* FASTA file upload

Input sequences are normalized by removing FASTA headers and whitespace and converting the sequence to uppercase.

The sequence is then validated to ensure that it contains only the standard DNA bases:

`A`, `T`, `C`, and `G`.

### 2. PAM Detection

The application scans the input sequence for PAM patterns associated with the selected Cas enzyme.

For example, SpCas9 currently supports:

* `NGG` as a canonical PAM
* `NAG` as a non-canonical PAM

PAM patterns are matched using IUPAC nucleotide rules.

### 3. Guide Identification

For every compatible PAM site, the application extracts the corresponding guide sequence based on the selected enzyme's guide length and PAM orientation.

Candidate guides are searched on both:

* Forward strand
* Reverse-complement strand

Each candidate stores information including its sequence, PAM, position, strand, enzyme, and PAM classification.

### 4. gRNA Scoring

Candidate guides are assigned a rule-based score starting at 100.

The current heuristic considers:

* GC content
* Guide position
* Repetitive nucleotide patterns
* PAM classification
* Potential off-target matches

The score is then used to rank candidate guides from highest to lowest.

### 5. Off-Target Checking

The current implementation performs a basic sequence-similarity search using Hamming distance.

Potential matches are identified when a DNA window differs from the guide by no more than a specified number of nucleotides.

The current default is:

```text
Maximum mismatches: 2
```

Both the input sequence and its reverse complement are searched.

This is a simplified educational implementation and is **not intended to replace established off-target prediction tools**.

## Results

The web interface displays ranked candidate guides along with information such as:

* Guide sequence
* PAM
* PAM classification
* GC content
* Number of potential off-target matches
* Score
* Position
* Strand
* Enzyme

The application can also export results as:

* CSV
* PDF

## Project Structure

```text
crispr-grna-designer/
│
├── app.py
├── crispr_logic.py
├── templates/
│   └── index.html
├── static/
│   └── ...
├── README.md
└── ...
```

## Example Workflow

```text
DNA sequence
      │
      ▼
Normalize & Validate
      │
      ▼
Select Cas Enzyme
      │
      ▼
Search Both DNA Strands
      │
      ▼
Identify PAM Sites
      │
      ▼
Extract Candidate gRNAs
      │
      ▼
Calculate Sequence Features
      │
      ▼
Estimate Potential Off-Targets
      │
      ▼
Calculate Rule-Based Score
      │
      ▼
Rank Candidates
      │
      ├──► Display in Web Interface
      ├──► Export CSV
      └──► Export PDF
```

## Current Limitations

This project is still under development. Some current limitations include:

* The scoring system is a simplified heuristic rather than a validated biological prediction model.
* Off-target analysis currently uses basic Hamming-distance matching.
* The application does not currently account for the full range of biological factors used by established gRNA design tools.
* The current implementation is designed primarily for learning and experimentation rather than experimental or clinical use.
* Additional testing and validation are still needed.

## What I Learned

This project has given me experience with:

* Applying programming concepts to biological data
* Processing and validating DNA sequences
* Working with reverse complements and sequence matching
* Implementing pattern matching using IUPAC nucleotide codes
* Designing rule-based scoring algorithms
* Building a Python/Flask web application
* Handling file uploads and user input
* Generating CSV and PDF reports
* Organizing application logic separately from the web interface
* Using Git and GitHub for version control

## Disclaimer

This is an educational and portfolio project. The results generated by this application should not be treated as experimentally validated recommendations or used for clinical or therapeutic decision-making.
