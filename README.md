# BHBfinder

BHBfinder is a pipeline for identifying Bulge–Helix–Bulge (BHB)-mediated splicing products and RNA ligation products from RNA-seq reads. It can be executed either in batch mode or step-by-step mode using Python commands or a shell script.

BHBfinder consists of five main modules. `normal_alignment` performs conventional read alignment with HISAT2. `spliced_alignment` identifies candidate hybrid reads through BLASTN alignment. `junction` detects BHB motifs and processes junctions. `set1` integrates junctions from all datasets of a species to generate Set 1. `set2` selects junctions for Set 2 and classifies ligation products into splicing, non-spliced ligation, and low-confidence products.


## Requirements

The following software should be installed and accessible through the system `PATH`.

- blast >= 2.12.0
- hisat2 = 2.2.1
- samtools = 1.5
- bedtools =  v2.26.0
- sra-tools  = 3.0.10
- Python 3.7+
- Python packages: `biopython`, `pyyaml`

```bash
pip install biopython pyyaml
```

You can also install the core tools with conda:

```bash
conda install hisat2=2.2.1 -c bioconda
conda install samtools=1.5 -c bioconda
conda install bedtools=2.26.0 -c bioconda
conda install blast=2.12.0 -c bioconda
conda install sra-tools=3.0.10 -c bioconda
conda install biopython=1.84 -c bioconda
conda install pyyaml=6.0.3 -c conda-forge
```

## Input

The correct file names and directory structure are essential for successful execution of BHBfinder. All files and directories for a species should be organized within a single working directory, with different species analyzed in separate working directories. For each species, the genome sequence and annotation files (.fna and .gff) should be placed in the `genome/` folder. Each RNA-seq dataset in FASTA format should be placed in a separate folder named after the corresponding data file without the file extension. The `config.yaml` file defines parameters for batch execution and should be located in the working directory. The `commands.sh` script enables step-by-step execution of individual modules.

```txt
work_dir/
├── config.yaml
├── commands.sh
├── genome/
│   ├── species.fna
│   └── species.gff
├── sample1/
│   └── sample1.fasta
└── sample2/
    └── sample2.fasta
```

Example `config.yaml`:

```yaml
species: S_solfataricus
genome_dir: genome
genome_file: GCF_000007005.1_ASM700v1_genomic.fna
gff_file: genomic.gff
specific: True
samples:
  - Ssol_data1
  - Ssol_data2

threads: 8
blast_threads: 20
```

`genome_dir` can be an absolute path or a path relative to the config file.

## Run with a Python command in batch mode

From the working directory:

```bash
python /path/to/BHBfinder.py --processing -c config.yaml
```

This is the same as:

```bash
python /path/to/BHBfinder.py --step all -c config.yaml
```

Config fields can also be overwritten from the command line:

```bash
python /path/to/BHBfinder.py --step all -c config.yaml \
  --species S_solfataricus \
  --genome_dir genome \
  --genome_file GCF_000007005.1_ASM700v1_genomic.fna \
  --gff_file genomic.gff \
  --samples Ssol_data1 Ssol_data2 \
  --threads 8 \
  --blast_threads 20
```

## Run with a Python command in step mode

Each processing step can also be run separately.

```bash
python /path/to/BHBfinder.py --step normal_alignment -c config.yaml
python /path/to/BHBfinder.py --step spliced_alignment -c config.yaml
python /path/to/BHBfinder.py --step junction -c config.yaml
python /path/to/BHBfinder.py --step set1 -c config.yaml
python /path/to/BHBfinder.py --step set2 -c config.yaml
```

## Run with a shell script

The shell script `commands.sh` allows users to execute individual BHBfinder modules with maximal control and flexibility. Before running the pipeline, configure the parameters at the beginning of the script according to the input data, reference genome, and desired analysis settings. To skip specific processing steps, comment out the corresponding commands.


```bash
bash commands.sh
```


## Files in working directory：

| File                                                         | Origin              | Content                                          | Usage                                       |
| ------------------------------------------------------------ | ------------------- | ------------------------------------------------ | ------------------------------------------- |
| `config.yaml`                                                | user                | configuration file                               | input for BHBfinder.py                      |
| `commands.sh`                                                | user                | configuration and commands                       | shell script for BHBfinder                  |
| `genome/species.fna`                                         | user                | genome sequence                                  |                                             |
| `genome/species.gff`                                         | user                | genome annotation                                |                                             |
| `genome/hisat2_*.ht2`                                        | `normal_alignment`  | HISAT2 alignment index                           | input for HISAT2                            |
| `genome/short_blast.*`                                       | `spliced_alignment` | BLAST alignment index                            | input for BLAST                             |
| `sample/sample.fasta`                                        | user                | read 1 data in fasta format                      | input for alignment                         |
| `sample/sample_normal_sorted.bam`                            | `normal_alignment`  | alignment file for normal reads                  | intermediate file for extraction            |
| `sample/sample_normal_m.bam`, `sample/sample_normal_m.bam.bai` | `normal_alignment`  | normal alignment file for properly paired  reads | IGV view and genome coverage                |
| `sample/sample_normal.bedgraph`                              | `normal_alignment`  | coverage of normal reads                         | input for calculating ligation rate         |
| `sample/sample_hisat2_normal_alignment.log`                  | `normal_alignment`  | HISAT2 log file                                  |                                             |
| `sample/sample_R1.tab`                                       | `spliced_alignment` | output from BLASTN alignment                     | input for junction identification           |
| `sample/sample_circ.bed`                                     | `junction`          | junction of individual hybrid in bed format      |                                             |
| `sample/sample_circ.out`                                     | `junction`          | junction of individual hybrid                    | intermediate results in junction processing |
| `sample/sample_junction.bed`                                 | `junction`          | processed junctions                              | merged into Set 1                           |
| `sample/sample_circ_R1.fasta`, `sample/sample_circ_R2.fasta` | `junction`          | two parts of hybrid reads as a fake read pair    | input for HISAT2 alignment                  |
| `sample/sample_BHB_sorted.bam`, `sample/sample_BHB_sorted.bam.bai` | `junction`          | sorted alignment file for junctions              | IGV view                                    |
| `sample/sample_BHB_C_sorted.bam`, `sample/sample_BHB_C_sorted.bam.bai` | `junction`          | alignment file for circular junctions            | IGV view                                    |
| `sample/sample_BHB_E_sorted.bam`, `sample/sample_BHB_E_sorted.bam.bai` | `junction`          | alignment file for linear junctions              | IGV view                                    |
| `SPECIES_Set1.bed`                                           | `set1`              | Set 1, merged junctions from all datasets        | processed into Set 2                        |
| `SPECIES_Set2.bed`                                           | `set2`              | Set 2, selected junctions with annotation        | final results of BHBfinder                  |
| `pipeline.log`                                               | BHBfinder           | log file                                         |                                             |




## Main Output

```txt
work_dir/
├── SPECIES_Set1.bed
├── SPECIES_Set2.bed
├── sample1/
│   ├── sample1_normal_m.bam
│   ├── sample1_R1.tab
│   ├── sample1_circ_R1.fasta
│   ├── sample1_circ_R2.fasta
│   ├── sample1_BHB_C_sorted.bam
│   ├── sample1_BHB_E_sorted.bam
│   └── sample1_junction.bed
└── pipeline.log
```

## Demo

```bash
cd demo
python ../BHBfinder.py --step all -c config.yaml
```

The demo reads are derived from SRA accessions `SRR12455257` and `SRR12455258`.

## Contact

- Author: Songlin Wu, 903277191@qq.com
- Author: Keqiong Ye, yekeqiong@ibp.ac.cn

## Citation

If you use BHBfinder in your work, please cite:

Songlin Wu, Yunxiao Yan, Lingfei Liang, Min Yue and Keqiong Ye. Large-scale analysis of archaeal transcriptomes reveals expanded roles of splicing in RNA decay and processing. (2026) To be published.
