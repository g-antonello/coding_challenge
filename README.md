# MAC Genotype & Resistance Caller 🧬

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A521.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-Enabled-blue.svg)](https://www.docker.com/)
[![Clinical Bioinformatics](https://img.shields.io/badge/Bioinformatics-Clinical-red.svg)]()

This Nextflow (DSL2) pipeline is designed for high-throughput genotyping and antibiotic resistance mutation identification in **Mycobacterium avium complex (MAC)** and other bacterial species. Specifically, it targets highly characterized coordinates in the **23S rRNA (rrl) gene** (equivalent to *Escherichia coli* positions **A2058** and **A2059**), which are major hotspots for mutations conferring resistance to macrolides (e.g., clarithromycin).

By leveraging a custom species-coordinate reference table and state-of-the-art variant calling tools (GATK Mutect2), the pipeline automatically aligns raw reads, calls variants, maps them to standard positions, and generates an actionable genotyping report for clinical or research samples.

---

## Quickstart

Ensure you have [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/) (or Apptainer) installed on your system. Docker is enabled by default to run all pipeline steps in a pre-configured, self-contained environment.

### 1. Prepare Input FASTQs
Copy your fastqs in a directory called `input/`, the pipeline will look for them there.

pseudocode:

```bash
mkdir input
cp path/to/fastq/folder/* input/.
```

### 2. Run the Pipeline
The pipeline supports execution profiles for running on local workstations or Slurm-managed high-performance computing (HPC) clusters.

#### Run Locally
For running on a personal computer or single workstation:
```bash
nextflow run main.nf -profile local
```

#### Run on a Slurm Cluster
For executing in an HPC environment:
```bash
nextflow run main.nf -profile slurm
```

---

## Pipeline Architecture

The workflow consists of five key automated steps:

```
  [Raw Reads] ----> QC & Trimming (fastp)
                           |
                           v
                     Alignment (minimap2) <--- [Indexed Reference Database]
                           |
                           v
                    Variant Calling (GATK Mutect2)
                           |
                           v
                Genotype Reporting (Rscript) <--- [Species Accession Reference Table]
                           |
                           v
                    [Variant calling Reports]
```

1. **Database Indexing (`INDEX_REF_DB`)**: Prepares reference databases by indexing FASTA sequences with `samtools faidx` and creating a sequence dictionary with `gatk CreateSequenceDictionary`.
2. **Quality Control (`QC_FASTP`)**: Trims adapters and filters out low-quality reads using `fastp` to ensure high-accuracy downstream mapping.
3. **Sequence Alignment (`ALIGN_MINIMAP2`)**: Maps reads against indexed reference genomes using `minimap2` and generates sorted, indexed BAM files via `samtools`.
4. **Variant Calling (`VARIANT_CALLING_GATK`)**: Calls variants directly from alignments using `GATK Mutect2`, which is optimized for bacterial species and somatic/subclonal variant detection.
5. **Report Generation (`REPORT`)**: An R script parsing the output VCF files, extracting coverage and allele frequency values, translating the variant positions to standard *E. coli* 23S rRNA numbers, and outputting an actionable CSV summary.

---

## Modifiable Parameters

All configuration options are defined in `nextflow.config` file. 
You can override these parameters at runtime by passing `--<param_name> <value>` 
to the `nextflow run` command.
Example, if the thresholding needs to change

```bash
nextflow run main.nf -profile local --min_read_depth 20
```

### Key Configuration Parameters

| Parameter | Default Value | Description |
| :--- | :--- | :--- |
| `params.input_dir` | `"input"` | Directory where input read files are stored. |
| `params.reads` | `"${params.input_dir}/*.{fastq,fq,fasta,fa}.gz"` | Glob pattern identifying raw sequence reads. |
| `params.outdir` | `"results"` | Output directory where processed data and reports are saved. |
| `params.reference` | `"data/ref_23S_database.fna"` | Path to the database of 23S rRNA FASTA sequences. |
| `params.table_references` | `"data/ref_genes_accessions.tsv"` | TSV mapping reference FASTA headers to E. coli equivalent positions. |
| `params.min_read_depth` | `10` | Minimum sequencing depth required at the variant position for a confident call. |
| `params.min_allele_freq` | `0.05` | Minimum allele frequency required to declare a resistance mutation. |
| `params.min_map_quality` | `20` | Minimum mapping quality threshold used during variant calling. |
| `params.min_base_quality` | `23` | Minimum base quality score (Q-score) required to consider a base. |

> **Base quality threshold selection:** The default `min_base_quality = 23` 
represents a base call error rate of $\approx 0.5\%$. This is chosen because we 
require at least 10x higher base call accuracy than the minimum allele frequency 
(`min_allele_freq = 0.05` or $5\%$) to confidently differentiate mutations from 
sequencing noise (derived via $-10 \log_{10}(5 \times 10^{-3}) \approx 23$).
---

## Output Directory Structure

Upon completion, results are organized inside the designated output directory (default: `results/`):

```
results/
├── reads_qc/
│   ├── <sample_id>.trimmed.fastq.gz   # High-quality trimmed sequencing reads
│   ├── <sample_id>.fastp.html          # QC HTML report for visual inspection
│   └── <sample_id>.fastp.json          # QC metrics in JSON format
├── alignments/
│   ├── <sample_id>.sorted.bam         # Coordinate-sorted BAM alignment file
│   └── <sample_id>.sorted.bam.bai     # BAM index file
├── variants/
│   └── <sample_id>.vcf                 # GATK Mutect2 variant calling file
└── reports/
    └── <sample_id>_report.csv          # genotyping and resistance report csv
```

### Understanding the Genotyping Report

The generated `<sample_id>_report.csv` contains detailed information on key hotspots:

* **`sample_id`**: Name of the processed sample.
* **`gene`**: Accession ID of the matched reference gene.
* **`ecoli_pos`**: Equivalent position in *E. coli* 23S rRNA (e.g., 2058 or 2059).
* **`ref_pos`**: Actual position in the species-specific reference FASTA database.
* **`ref_base` / `alt_base`**: Reference base vs. the called mutant base.
* **`alt_depth`**: Count of high-quality reads supporting the mutant allele.
* **`alt_freq`**: Allele frequency of the mutant allele (range `0.0 - 1.0`).
* **`mutation_detected`**: Formatted mutation description (e.g. `A2068G`).
* **`decision`**:
  * `PASS`: High-confidence variant. Meets both depth (`> min_read_depth`) and frequency (`> min_allele_freq`) thresholds.
  * `WARN`: Marginal call. Meets only one of the threshold metrics (e.g., high-frequency but low overall read coverage).
  * `FAIL`: Low-confidence call. Both depth and frequency fail to meet thresholds.
* **`notes`**: Displays `VARIANT_NOT_FOUND` if wildtype/no mutation is called.

---

## Reference Coordinate Mapping

The pipeline maps mutations from any custom database header back to *E. coli* numbering. 
This information is kept in `data/ref_genes_accessions.tsv`. 

The mapping database now includes:
* **E. coli 23S rRNA** (`E.coli_23S_rRNA`): Positions **2058** and **2059** map to `2058` and `2059`.
* **M. intracellulare** (`M.intracellulare_23S_rRNA`): Positions **2058** and **2059** map to **2268** and **2269**.

To extend the pipeline to support additional species or genes:
1. Append the new reference sequence to `data/ref_23S_database.fna`.
2. Add a corresponding entry in `data/ref_genes_accessions.tsv` specifying the fasta header name in `ID_in_fasta_reference_db` along with the positions and reference bases equivalent to *E. coli* A2058 and A2059.

---

## Containerization & Reproducibility

* **Docker container**: Runs fully containerized using the pre-built `giacomoa/coding_challenge:1.0` image containing all alignment, variant calling, and reporting tools.
* **Error Handling**: Configured with `process.maxErrors = 3` to automatically retry transiently failing processes up to three times with exponential backoff configurations if needed.
* **Hardware Allocations**:
  * Quality Control (`QC_FASTP`): 4 CPUs, 8 GB Memory
  * Alignment (`ALIGN_MINIMAP2`): 2 CPUs, 4 GB Memory
  * Variant Calling (`VARIANT_CALLING_GATK`): 2 CPUs, 4 GB Memory
  * Reporting (`REPORT`): 1 CPU, 2 GB Memory

# Next steps

* Generalize the pipeline to accept fasta files as input too. this means no QC applied in those cases
* take into account 23S copy number: how does that alter variant calling?
* pipeline unit tests missing. An important one
* improve reporting script structure