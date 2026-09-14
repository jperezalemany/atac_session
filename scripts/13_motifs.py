from pybedtools import BedTool
import pandas as pd
from pathlib import Path
import subprocess
import logging

logging.basicConfig(level=logging.INFO)

raw_dir = "/data/classes/seg_epig2026/data/atac_session/00_raw_data/"
meme_path = "/data/classes/seg_epig2026/software/meme-5.5.9/"
fasta = f"{raw_dir}/TAIR10_chr_all.fasta"
genome = f"{raw_dir}/TAIR10_chr_all.fasta.fai"
motif_db = f"{raw_dir}/JASPAR2026_CORE_plants_non-redundant_pfms.meme"

res_file = "results/11_differential/SIM6_vs_CIM7.results.csv"
summits_file = "/data/classes/seg_epig2026/data/atac_session/00_full_data/all_cutsites.summits.bed"
results_dir = "results/13_motifs"
logs_dir ="logs/13_motifs"


def get_summits_fasta(summits, peaks, outfile):
    (
        BedTool(summits)
        .filter(lambda x: x.name in peaks)
        .slop(b=100, g=genome)
        .sequence(fi=fasta, fo=outfile, name=True)
    )


def run_command(cmd):
    p = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    return p.wait()


def meme_sea_cmd(target, control, motifdb, outdir, log):
    return (
        f"{meme_path}/meme/bin/sea "
        f"--p {target} "
        f"--m {motifdb} "
        f"--n {control} "
        f"--oc {outdir} "
        "--qvalue "
        f"> {log} 2>&1"
    )


def main():
    Path(results_dir).mkdir(exist_ok=True)
    Path(logs_dir).mkdir(exist_ok=True)
    contrast = Path(res_file).name.removesuffix(".results.csv")

    logging.info("Reading dars")
    res = pd.read_csv(res_file, usecols=["peak", "dar", "padj"]).sort_values("padj").drop_duplicates()

    # 1. Get sequences of DAR summits

    # Random 500 non res (control set)
    logging.info("Selecting 500 random non DARs as control")
    control = set(res.loc[res.dar == "NS", "peak"].sample(500, random_state=42))
    logging.info("Extracting sequences of control regions")
    get_summits_fasta(summits_file, control, f"{results_dir}/{contrast}.control.fasta")

    # Top 500 up/down res
    for reg in ["up", "down"]:
        names = set(res.loc[res.dar == reg, "peak"].iloc[:500])
        logging.info(f"Extracting sequences of top 500 {reg} DARs")
        get_summits_fasta(summits_file, names, f"{results_dir}/{contrast}.{reg}.fasta")

    # 2. Peform motif enrichment (call MEME-SEA)
    for reg in ["up", "down"]:
        logging.info(f"Running MEME SEA over {reg} DARs sequences")
        cmd = meme_sea_cmd(
            target=f"{results_dir}/{contrast}.{reg}.fasta",
            control=f"{results_dir}/{contrast}.control.fasta",
            motifdb=motif_db,
            outdir=f"{results_dir}/{contrast}.{reg}",
            log=f"{logs_dir}/{contrast}.{reg}.meme_sea.log"
        )
        run_command(cmd)


if __name__ == "__main__":
    main()
