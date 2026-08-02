#!/usr/bin/env python3
import argparse
import yaml
import subprocess
import logging
import shutil
from pathlib import Path
from datetime import datetime
import sys

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('pipeline.log'),
            logging.StreamHandler()
        ]
    )


def load_config(config_path):
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        config = yaml.safe_load(f)
    # Convert string paths to Path objects
    paths = ['genome_dir']
    for path_key in paths:
        if path_key in config:
            path_value = Path(config[path_key])
            if not path_value.is_absolute():
                path_value = config_path.parent / path_value
            config[path_key] = path_value.resolve()
            # .resolve() for full path.

    return config


def str2bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).lower()
    if value in ['true', 't', 'yes', 'y', '1']:
        return True
    if value in ['false', 'f', 'no', 'n', '0']:
        return False
    raise argparse.ArgumentTypeError("Boolean value expected")


def validate_config(config):
    required_keys = ['species', 'genome_dir', 'genome_file', 'gff_file', 'samples']
    missing = [k for k in required_keys if k not in config]
    if missing:
        raise ValueError(f"Missing required config keys: {missing}")
    # Check gn
    if not config['genome_dir'].exists():
        raise FileNotFoundError(f"Genome directory not found: {config['genome_dir']}")
    if not (config['genome_dir'] / config['genome_file']).exists():
        raise FileNotFoundError(f"Genome file not found: {config['genome_dir'] / config['genome_file']}")
    if not (config['genome_dir'] / config['gff_file']).exists():
        raise FileNotFoundError(f"GFF file not found: {config['genome_dir'] / config['gff_file']}")
    if not config['samples']:
        raise ValueError("No samples found in config")


def run_command(cmd, cwd=None):
    """Execute shell command with error handling"""
    try:
        logging.info(f"Running command: {cmd}")
        subprocess.run(
            cmd, 
            shell=True, 
            check=True, 
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e}\nOutput:\n{e.output}")
        raise


def check_tools(tools):
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise FileNotFoundError(f"Required software not found in PATH: {', '.join(missing)}")


def check_file(file):
    if not Path(file).exists():
        raise FileNotFoundError(f"File not found: {file}")


def check_sample_fastas(samples):
    for sample in samples:
        check_file(Path(sample) / f"{sample}.fasta")


def get_paths(config):
    find_circ_script = get_script_path("find_circrna_v3.py")
    merge_script = get_script_path('merge_junction.py')
    get_set2_script = get_script_path('get_set2.py')
    paths = {
        'genome_dir': config['genome_dir'],
        'genome': config['genome_dir'] / config['genome_file'],
        'gff': config['genome_dir'] / config['gff_file'],
        'hisat_index': config['genome_dir'] / "hisat2",
        'blast_db': config['genome_dir'] / "short_blast",
        'parent_dir': str(Path.cwd()),
        #
        'script': find_circ_script,
        'merge_script': merge_script,
        'get_set2_script': get_set2_script
    }
    return paths


def check_database_status(paths):
    """Check the status of genome index database"""
    status = {'hisat2': False, 'blast': False}

    # hisat2 index
    hisat_files = list(paths['hisat_index'].parent.glob(f"{paths['hisat_index'].name}*.ht2"))
    status['hisat2'] = len(hisat_files) > 0

    # blast index
    blast_files = list(paths['blast_db'].parent.glob(f"{paths['blast_db'].name}.nhr"))
    status['blast'] = len(blast_files) > 0

    return status


def build_hisat_database(paths):
    """Build HISAT2 index database"""
    genome_file = paths['genome']
    genome_dir = paths['genome_dir']

    if not check_database_status(paths)['hisat2']:
        logging.info("Building HISAT2 index...")
        hisat_cmd = f"hisat2-build {genome_file} {paths['hisat_index']}"
        run_command(hisat_cmd, genome_dir )
    else:
        logging.info("HISAT2 index already exists, skipping build")


def build_blast_database(paths):
    """Build BLAST index database"""
    genome_file = paths['genome']
    genome_dir = paths['genome_dir']

    if not check_database_status(paths)['blast']:
        logging.info("Building BLAST database...")
        blast_cmd = f"makeblastdb -in {genome_file} -dbtype nucl -out {paths['blast_db']}"
        run_command(blast_cmd, genome_dir)
    else:
        logging.info("BLAST database already exists, skipping build")


def normal_alignment(samples, config, paths):
    """
    None-spliced alignment with HISAT2.
    Input: sample/sample.fasta and genome fasta.
    Output: sample/sample_normal_m.bam and sample/sample_normal.bedgraph.
    """
    check_tools(['hisat2-build', 'hisat2', 'samtools', 'bedtools'])
    check_sample_fastas(samples)
    build_hisat_database(paths)

    for sample in samples:
        sample_dir = Path(sample)
        logging.info(f"Processing sample: {sample}")

        ### HISAT2 alignment
        hisat_cmd = f"""hisat2 -p {config.get('threads', 4)} \
            -f --no-spliced-alignment \
            -x {paths['hisat_index']} \
            -U {sample}.fasta \
            -S {sample}_normal.sam
        """
        run_command(hisat_cmd, sample_dir)

        ### SAMtools processing
        samtools_cmds = [
            f"samtools sort {sample}_normal.sam > {sample}_normal_sorted.bam",
            f"samtools flagstat {sample}_normal_sorted.bam > {sample}_hisat2_normal_alignment.log",
            f"samtools view -h -b -F 4 {sample}_normal_sorted.bam > {sample}_normal_m.bam",
            f"samtools index {sample}_normal_m.bam"
        ]
        for cmd in samtools_cmds:
            run_command(cmd, sample_dir)

        ### Bedtools processing
        # 0-coor-based, all none-zero position
        bedtools_cmd = f"bedtools genomecov -dz -ibam {sample}_normal_m.bam > {sample}_normal.bedgraph"
        run_command(bedtools_cmd, sample_dir)

        ### Remove intermediate files
        rm_cmd = f"rm {sample}_normal.sam"
        run_command(rm_cmd, sample_dir)


def spliced_alignment(samples, config, paths):
    """
    Spliced alignment with BLAST.
    Input: sample/sample.fasta and genome fasta.
    Output: sample/sample_R1.tab.
    """
    check_tools(['makeblastdb', 'blastn'])
    check_sample_fastas(samples)
    build_blast_database(paths)

    for sample in samples:
        sample_dir = Path(sample)
        logging.info(f"Processing sample: {sample}")

        ### BLAST alignment
        # Time consuming
        blast_cmd = f"""blastn -db {paths['blast_db']} \
            -num_threads {config.get('blast_threads', 20)} \
            -task "blastn-short" \
            -evalue 1e-3 -strand both \
            -query {sample}.fasta \
            -out {sample}_R1.tab \
            -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
        """
        run_command(blast_cmd, sample_dir)


def find_junction(samples, config, paths):
    """
    Find BHB junctions.
    Input: sample/sample_R1.tab and sample/sample_normal.bedgraph.
    Output: sample/sample_junction.bed and sample/sample_BHB*.bam.
    """
    check_tools(['hisat2', 'samtools', 'awk'])
    check_sample_fastas(samples)
    if not check_database_status(paths)['hisat2']:
        raise FileNotFoundError("HISAT2 index not found, please run --step normal_alignment first")

    for sample in samples:
        sample_dir = Path(sample)
        check_file(sample_dir / f"{sample}_R1.tab")
        check_file(sample_dir / f"{sample}_normal.bedgraph")
        logging.info(f"Processing sample: {sample}")

        ### Python analysis steps
        analysis_cmds = [
            f"python {paths['script']} -a circ -basename {sample} -r1 {sample}.fasta -tab1 {sample}_R1.tab",
            f"python {paths['script']} -a junction --anti -basename {sample} -g {paths['genome']} -gff {paths['gff']} -c {sample}_circ.out -j1 {sample}_circ_R1.fasta -j2 {sample}_circ_R2.fasta -bedgraph {sample}_normal.bedgraph"
        ]
        for cmd in analysis_cmds:
            run_command(cmd, sample_dir)

        ### Hisat alignment for igv view
        igv_cmds = [
            f"hisat2  -f --no-spliced-alignment -x {paths['hisat_index']} -1 {sample}_BHB_R1.fasta -2 {sample}_BHB_R2.fasta -S {sample}_BHB.sam",
            f"samtools sort {sample}_BHB.sam > {sample}_BHB_sorted.bam; samtools index {sample}_BHB_sorted.bam",
            "awk '$1~/^@/ || $1~/_E$/ {print}' " +  f"{sample}_BHB.sam | samtools sort > {sample}_BHB_E_sorted.bam; samtools index {sample}_BHB_E_sorted.bam",
            "awk '$1~/^@/ || $1~/_C$/ {print}' " +  f"{sample}_BHB.sam | samtools sort > {sample}_BHB_C_sorted.bam; samtools index {sample}_BHB_C_sorted.bam",
            f"samtools index {sample}_BHB_sorted.bam"
        ]
        for cmd in igv_cmds:
            run_command(cmd, sample_dir)

        ### Remove intermediate files
        rm_cmd = f"rm {sample}_BHB.sam"
        run_command(rm_cmd, sample_dir)


def get_script_path(script_name):
    """Full path of scripts"""
    current_dir = Path(__file__).resolve().parent
    script_path = current_dir / "scripts" / script_name

    if not script_path.exists():
        print(f"File not exist {script_path}")
        sys.exit(1)

    return script_path


def generate_example_config():
    """Config example"""
    example = """# Config file (analysis_config.yaml)
species: S_solfataricus
genome_dir: genome
genome_file: GCF_000007005.1_ASM700v1_genomic.fna
gff_file: genomic.gff
specific: True  # If the input data is strand-specific or not, True or False
samples:
  - Ssol_data1
  - Ssol_data2

# Optional params
threads: 8                  # Hisat2 CPU number
blast_threads: 20           # Blast CPU number
"""
    with open("analysis_config.example.yaml", "w") as f:
        f.write(example)
    print("Make config file: analysis_config.example.yaml")


def merge_junction(samples, species, paths):
    """Merge multiple sample files"""
    junction_file_list = []
    parent_dir = paths['parent_dir']
    for sample in samples:
        junction_file = Path(sample) / f"{sample}_junction.bed"
        check_file(junction_file)
        junction_file_list.append(str(junction_file))
    merge_cmd = f"python {paths['merge_script']} --beds {' '.join(junction_file_list)} -o {species}_Set1.bed"
    run_command(merge_cmd, parent_dir)


def get_set2(config, paths):
    set1_file = f"{config['species']}_Set1.bed"
    check_file(Path(paths['parent_dir']) / set1_file)
    specific = config.get('specific', True)
    set2_cmd = f"python {paths['get_set2_script']} --species {config['species']} --set1 {set1_file} --specific {specific}"
    run_command(set2_cmd, paths['parent_dir'])
    return


def main():
    parser = argparse.ArgumentParser(description="CircRNA Analysis Pipeline")
    parser.add_argument("-c", "--config", type=Path, help="Path to YAML config file")
    parser.add_argument("--step",
                        choices=['all', 'normal_alignment', 'spliced_alignment', 'junction', 'set1', 'set2'],
                        default='all',
                        help="Run one step only. Default: all")
    parser.add_argument("--processing", action="store_true",
                        help='Run the whole pipeline. Same as --step all')
    parser.add_argument("--generate-config", action="store_true",
                        help="Generate config and exit")
    parser.add_argument("--species", help="Species name, overwrite config")
    parser.add_argument("--genome_dir", type=Path, help="Genome directory, overwrite config")
    parser.add_argument("--genome_file", help="Genome fasta file name, overwrite config")
    parser.add_argument("--gff_file", help="Genome gff file name, overwrite config")
    parser.add_argument("--specific", type=str2bool, help="If the input data is strand-specific, overwrite config")
    parser.add_argument("--samples", nargs='+', help="Sample names, overwrite config")
    parser.add_argument("--threads", type=int, help="Hisat2 CPU number, overwrite config")
    parser.add_argument("--blast_threads", type=int, help="Blast CPU number, overwrite config")
    args = parser.parse_args()

    if args.generate_config:
        generate_example_config()
        return

    setup_logging()

    try:
        ### configuring
        logging.info("Config...")
        if args.config is None:
            config = {}
        else:
            config = load_config(args.config)

        if args.species is not None:
            config['species'] = args.species
        if args.genome_dir is not None:
            genome_dir = args.genome_dir
            if not genome_dir.is_absolute():
                genome_dir = Path.cwd() / genome_dir
            config['genome_dir'] = genome_dir.resolve()
        if args.genome_file is not None:
            config['genome_file'] = args.genome_file
        if args.gff_file is not None:
            config['gff_file'] = args.gff_file
        if args.specific is not None:
            config['specific'] = args.specific
        if args.samples is not None:
            config['samples'] = args.samples
        if args.threads is not None:
            config['threads'] = args.threads
        if args.blast_threads is not None:
            config['blast_threads'] = args.blast_threads

        validate_config(config)
        paths = get_paths(config)

        ### run steps
        steps = ['normal_alignment', 'spliced_alignment', 'junction', 'set1', 'set2'] if args.processing or args.step == 'all' else [args.step]
        start_time = datetime.now()
        for step in steps:
            if step == 'normal_alignment':
                logging.info("none-spliced alignemnt with Hisat2")
                normal_alignment(config['samples'], config, paths)
            elif step == 'spliced_alignment':
                logging.info("Spliced alignment with blast")
                spliced_alignment(config['samples'], config, paths)
            elif step == 'junction':
                logging.info("Find BHB junctions")
                find_junction(config['samples'], config, paths)
            elif step == 'set1':
                logging.info("Processing Set 1...")
                merge_junction(config['samples'], config['species'], paths)
            elif step == 'set2':
                logging.info("Processing Set 2...")
                get_set2(config, paths)
        duration = datetime.now() - start_time
        logging.info(f"Pipeline completed successfully in {duration}")
        return

    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        exit(1)


if __name__ == "__main__":
    main()
