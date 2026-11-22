# Maestro

**Heterogeneity-aware network analysis of molecular expression data**

Maestro identifies coordinated expression patterns across molecular features and evaluates their association with a phenotype. It uses the **Duo** correlation measure to construct an expression-state network, identifies connected components as candidate modules, and tests those modules in separate discovery and validation cohorts.

Unlike conventional co-expression networks, Maestro represents each analyte with two nodes: one for **high expression** and one for **low expression**. This enables it to distinguish high–high, high–low, low–high, and low–low associations, including patterns observed in subsets of samples.

> **Scope:** Network edges indicate correlations between expression states. They do not, by themselves, establish direct gene regulation or causal interactions.

## Contents

- [Overview](#overview)
- [Software requirements](#software-requirements)
- [Installation](#installation)
- [Package structure](#package-structure)
- [Input data](#input-data)
- [Analysis workflow](#analysis-workflow)
  - [1. Prepare and normalize data](#1-prepare-and-normalize-data)
  - [2. Create discovery and validation cohorts](#2-create-discovery-and-validation-cohorts)
  - [3. Construct the Duo network](#3-construct-the-duo-network)
  - [4. Retain the highest-weight edges](#4-retain-the-highest-weight-edges)
  - [5. Identify modules](#5-identify-modules)
  - [6. Test discovery associations](#6-test-discovery-associations)
  - [7. Validate significant modules](#7-validate-significant-modules)
- [Command reference](#command-reference)
- [Output files](#output-files)
- [Configuration and reproducibility](#configuration-and-reproducibility)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [License](#license)

## Overview

Maestro is designed for expression datasets in which a coordinated molecular signature may be present in only some individuals. Duo evaluates four distinct expression-state relationships for each analyte pair:

| Analyte A | Analyte B | Relationship |
| --- | --- | --- |
| High | High | High–high |
| High | Low | High–low |
| Low | High | Low–high |
| Low | Low | Low–low |

The resulting network has **two nodes per analyte**. Statistically supported relationships form weighted edges. After the strongest edges are selected, a breadth-first search (BFS) finds connected components, each of which defines a candidate expression module.

The network is constructed **without using case/control labels**. Modules are subsequently tested for association with the phenotype in the discovery cohort; significant discovery patterns are evaluated in a held-out validation cohort.

**Workflow:**

```text
Expression matrix
       │
       ▼
Data cleaning and normalization
       │
       ▼
Discovery / validation split
       │
       ▼
Duo network construction  ──►  Weighted GML network
       │
       ▼
keepHi: strongest edges
       │
       ▼
BFS: connected components  ──►  Candidate modules
       │
       ▼
d2d: discovery associations
       │
       ▼
validate: held-out associations
```

## Software requirements

- A Linux or Unix-like system with a C++ compiler (`g++`)
- GNU Make
- Standard shell utilities (`tar`, `wc`, `head`, etc.)
- Sufficient memory and CPU time for the chosen dataset; Duo network construction and module association testing can be resource-intensive for large datasets

The analysis tools are command-line programs written in C++. No Python environment is required to build or run the core pipeline.

## Installation

Clone the repository and extract the versioned source package:

```bash
git clone https://github.com/climers/Maestro.git
cd Maestro
tar -xzf Maestro_v1.0.tar.gz
cd Maestro_v1.0
```

Compile the core analysis programs:

```bash
for dir in normalize duo keepHi bfs d2d validate; do
    make -C "$dir"
done
```

Compile the optional data-preparation utilities:

```bash
for dir in tools/transpose tools/extractRows tools/dataRetainer; do
    make -C "$dir"
done
```

Each executable is written to its own directory (for example, `duo/duo` or `validate/validate`). Run the commands below from the extracted `Maestro_v1.0` directory, adjusting input paths to your files.

## Package structure

```text
Maestro_v1.0/
├── normalize/             Normalize expression values
├── duo/                   Construct a Duo correlation network
├── keepHi/                Select highest-weight network edges
├── bfs/                   Identify connected components
├── d2d/                   Test discovery associations
├── validate/              Evaluate patterns in held-out samples
└── tools/
    ├── transpose/         Transpose expression matrices
    ├── extractRows/       Select samples using identifier lists
    ├── dataRetainer/      Inspect and reduce missing data
    ├── replaceEmpty/      Replace empty fields
    └── rmSingletons/      Remove isolated nodes from GML files
```

Each compiled program is run independently; input and output files link the stages together.

## Input data

Maestro expects **whitespace-delimited numerical expression matrices** with:

- **Rows:** analytes (genes, transcripts, miRNAs, or other measured features)
- **Columns:** individuals or samples
- **Optional header rows:** sample metadata
- **Optional header columns:** analyte identifiers

Example with **one header row** and **one header column**:

```text
AnalyteID   Sample_01   Sample_02   Sample_03
Gene_A      5.21        4.13        6.07
Gene_B      2.84        3.09        2.45
Gene_C      8.17        7.62        8.45
```

The `numAnalytes` and `numSamples` arguments count **data rows and data columns only**; header rows and columns are specified separately.

For an end-to-end analysis, prepare the following sample groups:

| Input matrix | Used for |
| --- | --- |
| Combined **discovery cases + controls** | Duo network construction |
| **Discovery cases** | Discovery association testing |
| **Discovery controls** | Discovery association testing |
| **Validation cases** | Independent validation |
| **Validation controls** | Independent validation |

All five matrices must retain the **same analytes in the same order**. Discovery and validation cohorts must contain **different samples**. The package does not require phenotype labels in the matrix used for Duo network construction.

## Analysis workflow

The examples below assume that the package is compiled and commands are run from `Maestro_v1.0/`. Substitute your own file paths and numeric dimensions for the uppercase arguments.

### 1. Prepare and normalize data

Clean the expression matrix, identify the relevant case and control samples, and address missing values. The included `transpose`, `extractRows`, and `dataRetainer` programs can assist with these tasks.

Normalize the **combined case/control matrix** before splitting it into the discovery and validation groups:

```bash
./normalize/normalize \
  cases_controls.txt NUM_ANALYTES NUM_SAMPLES HEADER_ROWS HEADER_COLS \
  0 1 cases_controls_normalized.txt
```

This program processes each analyte row, winsorizes extreme values using its interquartile range, and scales the values into the specified target range (here, `0` to `1`).

### 2. Create discovery and validation cohorts

Divide the normalized samples into nonoverlapping discovery and validation cohorts. Where available, account for important covariates (such as age and sex) when making the split.

Extract four matrices: discovery cases, discovery controls, validation cases, and validation controls. Combine the **discovery cases and controls** into a single matrix while preserving analyte order. The discovery combined file is used only to construct the network; the separate case/control matrices are used for association testing.

Suggested filenames:

```text
discovery_combined.txt
discovery_cases.txt
discovery_controls.txt
validation_cases.txt
validation_controls.txt
```

**No separate discretization step is needed.** In the distributed Duo implementation, high- and low-expression thresholds are determined within the network-construction program.

### 3. Construct the Duo network

```bash
./duo/duo \
  discovery_combined.txt discovery.gml \
  NUM_DISCOVERY_SAMPLES NUM_ANALYTES HEADER_COLS HEADER_ROWS 0.7
```

`0.7` is an illustrative minimum Duo score for retaining candidate relationships. The program evaluates expression-state associations, applies permutation-based edge screening, and writes a weighted network in **GML** format.

The network contains `2 × NUM_ANALYTES` nodes: one HIGH and one LOW state for every analyte. Use the generated `*.duo.log` file to find the number of retained edges for the next step.

### 4. Retain the highest-weight edges

```bash
./keepHi/keepHi \
  discovery.gml NUM_NODES NUM_DUO_EDGES 1000 discovery_top.gml
```

- `NUM_NODES` must equal **twice the number of analytes**.
- `NUM_DUO_EDGES` is the edge count reported by Duo.
- `1000` is the **requested** number of edges to retain.

**Tied edge weights are retained**, so the reduced network may contain more than 1,000 edges. Check the `*.keepHi.log` file for the actual count and cutoff.

### 5. Identify modules

```bash
./bfs/bfs discovery_top.gml discovery_top1k.bfs
```

BFS identifies the connected components of the reduced graph. Each non-singleton component becomes a candidate multi-analyte expression pattern. The `*.bfs.log` file reports the number of isolated nodes, two-node components, components with three or more nodes, and the largest component size.

### 6. Test discovery associations

```bash
./d2d/d2d \
  discovery_cases.txt discovery_controls.txt \
  NUM_DISC_CASES NUM_DISC_CONTROLS NUM_ANALYTES HEADER_ROWS HEADER_COLS \
  discovery_results.txt discovery_top1k.bfs
```

`d2d` tests whether each selected module's expression pattern differs between discovery cases and controls. Its outputs include pattern definitions, case/control counts, association statistics, and significant patterns.

The `discovery_top1k.bfs` argument specifies which network modules to test. The `*.sig` file provides the significant discovery patterns and their selected thresholds for validation.

### 7. Validate significant modules

```bash
./validate/validate \
  validation_cases.txt validation_controls.txt \
  NUM_VAL_CASES NUM_VAL_CONTROLS NUM_ANALYTES HEADER_ROWS HEADER_COLS \
  discovery_results.txt.sig validation_results.txt discovery_top1k.bfs
```

`validate` assesses the **same discovery-selected patterns and thresholds** in the independent validation cohort. Do not use validation results to redefine modules or select new thresholds.

## Command reference

All commands below use paths relative to `Maestro_v1.0/`.

| Program | Command syntax |
| --- | --- |
| `normalize` | `./normalize/normalize INPUT N_ANALYTES N_SAMPLES H_ROWS H_COLS MIN MAX OUTPUT` |
| `duo` | `./duo/duo INPUT OUTPUT.gml N_SAMPLES N_ANALYTES H_COLS H_ROWS THRESHOLD` |
| `keepHi` | `./keepHi/keepHi INPUT.gml N_NODES N_EDGES N_KEEP OUTPUT.gml` |
| `bfs` | `./bfs/bfs INPUT.gml OUTPUT.bfs` |
| `d2d` | `./d2d/d2d CASES CONTROLS N_CASES N_CONTROLS N_ANALYTES H_ROWS H_COLS OUTPUT [MODULES.bfs]` |
| `validate` | `./validate/validate CASES CONTROLS N_CASES N_CONTROLS N_ANALYTES H_ROWS H_COLS DISCOVERY.sig OUTPUT [MODULES.bfs]` |

**Header argument order matters:** `duo` uses **header columns before header rows**. `normalize`, `d2d`, and `validate` use **header rows before header columns**.

The `.bfs` argument is optional for `d2d` and `validate`; include it when evaluating the modules produced by the network workflow. When the argument is not used, associations for each analyte, one at a time, will be computed. 

## Output files

| Output | Description |
| --- | --- |
| `*.gml` | Full or reduced weighted network in GML format |
| `*.duo.log` | Duo run configuration, permutation details, and edge count |
| `*.keepHi.log` | Edge-selection cutoff and number of retained edges |
| `*.bfs` | Network components used for downstream testing |
| `*.bfs.log` | Counts and sizes of connected components |
| `discovery_results.txt` | Discovery pattern association results |
| `discovery_results.txt.sig` | Significant discovery patterns and thresholds |
| `validation_results.txt` | Held-out association results |
| `validation_results.txt.sig` | Significant validation results |
| `*.log` | Execution parameters and summary statistics produced by each program |

Association outputs may report cases and controls meeting a pattern, odds ratios with 95% confidence intervals, Youden's J, and p-values. These are **statistical expression associations**; they require further biological validation before being interpreted as mechanisms or clinical biomarkers.

## Configuration and reproducibility

Record enough information to reproduce each run:

- Dataset version, feature/sample counts, header counts, and missing-value handling
- Sample IDs and discovery/validation assignment
- Normalization range and preprocessing decisions
- Duo threshold, permutation settings, and retained edge count
- Exact commands, package version, and output log files

The source header `duo/duo.h` contains configuration constants, including `NUMRANDTRIALS`, which determines the number of permutation trials used for edge validation. Rebuild Duo with `make -C duo` after changing its header configuration.

Large analyses can be submitted as batch jobs on a computing cluster. Record the software version and hardware resources when comparing results across systems.

## Troubleshooting

**The program reports unexpected matrix dimensions.** Confirm that sample and analyte counts exclude all headers. Check that every data row has the same number of fields and that the header-row and header-column arguments are in the correct order for the chosen program.

**No or very few edges are produced.** Check the Duo threshold and inspect `*.duo.log`. A lower threshold may admit additional edges, but thresholds should be recorded and selected consistently with the analysis plan.

**The reduced network has more edges than requested.** `keepHi` retains edges tied at its cutoff. The log reports the actual number retained.

**An input value exceeds the supported range.** Some programs enforce configurable numeric bounds in their source headers. Check the relevant error message and input data before changing such bounds, then recompile if necessary.

**A computation is terminated for insufficient memory.** Larger numbers of analytes and larger network modules can substantially increase resource requirements, especially during `d2d` and `validate`. Inspect system/job logs and use adequate compute resources.

**No patterns replicate in validation.** A completed validation run may produce an empty `*.sig` file. Check the companion `*.log` and full output before interpreting the result; an empty significant-results file is not the same as a missing run.

## References

- Climer, S. (2021). [Connecting the dots: The boons and banes of network modeling](https://doi.org/10.1016/j.patter.2021.100374). *Patterns*.
- Climer, S., et al. (2020). [Synchronized genetic activities in Alzheimer's brains revealed by heterogeneity-capturing network analysis](https://doi.org/10.1101/2020.01.28.923730).

## License

Maestro is released under the [GNU General Public License v3.0](LICENSE).
