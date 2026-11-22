# Maestro

**Network-based discovery of coordinated molecular expression patterns**

Maestro identifies groups of analytes with coordinated expression patterns, including relationships present in only a subset of samples. It uses **Duo correlation** to build an expression-state network, identifies connected components as candidate modules, and tests their association with a phenotype in independent discovery and validation cohorts.

Each analyte is represented by two network nodes, corresponding to **high** and **low** expression. Duo evaluates high–high, high–low, low–high, and low–low relationships separately. Network edges represent expression correlations, not necessarily direct regulatory interactions.

## Installation

**Requirements:** Linux or another Unix-like environment, a C++ compiler (`g++`), GNU Make, and standard command-line utilities. Large datasets may require substantial memory and compute resources.

```bash
git clone https://github.com/climers/Maestro.git
cd Maestro
tar -xzf Maestro_v1.0.tar.gz
cd Maestro_v1.0

for dir in tools/transpose tools/extractRows tools/dataRetainer normalize duo keepHi bfs d2d validate; do
    make -C "$dir"
done
```

The compiled executables are located in their corresponding directories. The `tools/` directory also contains utilities for preparing expression matrices.

## Input data

Maestro uses whitespace-delimited expression matrices with **analytes in rows** and **samples in columns**. Rows and columns containing identifiers or metadata are specified with header counts.

Example with one header row and one header column:

```text
AnalyteID    Sample_1    Sample_2    Sample_3
Gene_A       5.21        4.13        6.07
Gene_B       2.84        3.09        2.45
Gene_C       8.17        7.62        8.45
```

Record the number of analytes, samples, header rows, and header columns. **Data counts exclude headers.** All matrices within an analysis must contain the same analytes in the same order.

Prepare the following datasets:

| Dataset | Purpose |
| --- | --- |
| Combined discovery cases and controls | Duo network construction |
| Discovery cases | Module association testing |
| Discovery controls | Module association testing |
| Validation cases | Independent validation |
| Validation controls | Independent validation |

Discovery and validation samples must not overlap. Case/control labels are **not used to construct the network**; they are used to assess module associations afterward.

## Analysis workflow

### 1. Prepare and normalize expression data

Clean the expression matrix, identify the relevant cases and controls, and check for missing values. The included `tools/transpose`, `tools/extractRows`, and `tools/dataRetainer` utilities support data preparation and sample selection.

Normalize the cleaned matrix to the range **0–1**:

```bash
./normalize/normalize input.txt NUM_ANALYTES NUM_SAMPLES HEADER_ROWS HEADER_COLS 0 1 normalized.txt
```

Then create nonoverlapping discovery and validation cohorts, preserving analyte order. Combine the normalized discovery cases and controls into one matrix for network construction. Duo determines high- and low-expression states internally.

### 2. Construct the Duo network

```bash
./duo/duo discovery_combined.txt discovery.gml NUM_DISCOVERY_SAMPLES NUM_ANALYTES HEADER_COLS HEADER_ROWS 0.7
```

`0.7` is the minimum Duo correlation threshold in this example. The output is a weighted network in GML format with **two nodes per analyte**. The accompanying `.duo.log` reports the number of edges produced.

### 3. Retain the highest-weight edges

```bash
./keepHi/keepHi discovery.gml NUM_NODES NUM_DUO_EDGES 1000 discovery_top.gml
```

Set `NUM_NODES` to `2 × NUM_ANALYTES` and `NUM_DUO_EDGES` to the edge count reported by Duo. `1000` is the requested number of edges. **Edges tied at the cutoff are retained**, so the output may contain more than 1,000 edges.

### 4. Identify network modules

```bash
./bfs/bfs discovery_top1k.gml discovery_top1k.bfs
```

Breadth-first search identifies connected components, which are treated as candidate expression modules. The `.bfs.log` reports singleton nodes, two-node components, larger components, and the largest module size.

### 5. Test associations in discovery data

```bash
./d2d/d2d discovery_cases.txt discovery_controls.txt NUM_DISC_CASES NUM_DISC_CONTROLS NUM_ANALYTES HEADER_ROWS HEADER_COLS discovery_results.txt discovery_top1k.bfs
```

`d2d` evaluates whether the network modules distinguish cases from controls. Use the **normalized expression matrices**. The `.sig` output contains significant discovery patterns and their thresholds.

### 6. Validate significant patterns

```bash
./validate/validate validation_cases.txt validation_controls.txt NUM_VAL_CASES NUM_VAL_CONTROLS NUM_ANALYTES HEADER_ROWS HEADER_COLS discovery_results.txt.sig validation_results.txt discovery_top1k.bfs
```

Validation evaluates discovery-significant modules using the held-out cohort and the previously selected thresholds. Do not use validation results to redefine the modules or tune discovery thresholds.

**Argument order:** `duo` takes `HEADER_COLS` before `HEADER_ROWS`. `normalize`, `d2d`, and `validate` take `HEADER_ROWS` before `HEADER_COLS`.

## Outputs

| File | Description |
| --- | --- |
| `*.gml` | Weighted Duo network or reduced network |
| `*.bfs` | Connected components (modules) |
| `*.log` | Parameters, execution details, and summary statistics |
| `discovery_results.txt` | Discovery module association results |
| `discovery_results.txt.sig` | Significant discovery patterns and thresholds |
| `validation_results.txt` | Association results in held-out samples |
| `validation_results.txt.sig` | Significant validated patterns |

Association results include the number of cases and controls exhibiting a pattern, odds ratios with confidence intervals, Youden's J, and p-values. A statistically associated module is a candidate expression signature; further biological or experimental work is required to establish its mechanism.

## Reproducibility

Retain the input data version, sample assignments, normalization settings, command-line arguments, Duo threshold, permutation settings, and output log files. The number of permutation trials for network construction can be configured in `duo/duo.h` (`NUMRANDTRIALS`); recompile after changing it.

## References

- Climer, S. (2021). [Connecting the dots: The boons and banes of network modeling](https://doi.org/10.1016/j.patter.2021.100374). *Patterns*.
- Climer, S., et al. (2020). [Synchronized genetic activities in Alzheimer's brains revealed by heterogeneity-capturing network analysis](https://doi.org/10.1101/2020.01.28.923730).

## License

Maestro is distributed under the [GNU General Public License v3.0](LICENSE).
