# ============================================================
# Simulation for EC with the NB model (plan: article/experiments_plan.txt, C8).
# 5 Mb region of chr7: the reference is the hg38 region; the normal consists of the two corresponding
# haplotypes of the HG008-N verkko 2.2 (hap1 + hap2). Known somatic variants on hap1, inside
# the region (no translocations and no events extending beyond the 5 Mb), in two clones:
#   hap1_T1 = hap1 + clonal variants; hap1_T2 = hap1 + clonal + subclonal.
# Tumor at depth C and purity rho, subclone in 50% of the tumor cells:
#   hap2 C/2; hap1 C/2 (1-rho); hap1_T1 C/2 rho/2; hap1_T2 C/2 rho/2
#   -> clonal VAF rho/2, subclonal rho/4.
# Normal: hap1 + hap2, 17.5x each, simulated separately (reads independent of those of the tumor).
# Grid (2026-09-29): C = 90, 60; rho = 1.0, 0.6, 0.3; EC of the tumor with the NB model
# (--tumor-purity true rho, --clonality 0.5, --somatic-fnr 0.05 / 0.10 / 0.20) and classic hifiasm.
# HiFi reads: PBSIM3 (QSHMM model) + ccs. Release binaries v1.0.0-somatic; the provenance dump
# (cluster and candidate rows) with the binary of the call-provenance branch, as for tab:res-bp-coverage.
# ============================================================
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
A=$HG008
SRC=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")   # script folder (differs from D with SIM_DIR)
D=${SIM_DIR:-$A/sim_ec_v2_2026-09-29}   # v2 (2026-09-29): SNV/indels visible after HPC, INS > 500 bp copies only, HP tags fixed, attributed measure. v1: sim_ec_2026-09-29 (scripts here)                # SIM_DIR: working folder (a test run goes in its own folder)
SV=$SV_BIN                                  # minimap2, samtools, bcftools, python with pysam, k8
SIMENV=$ENVS/simhifi                                     # pbsim3 + ccs (00_env.sh)
TAGBAM=$LIB/tag_bam_from_list.py
S2G=$LIB/sfs_to_genome.py
BPTAB=$SRC/bp_point_table_prov.py
PREP=$SRC/prep_vcf.py
HG38=$A/refs/hg38.primary.fa.gz
VERKKO=$A/refs/HG008N_verkko2.2_20250218/diploid.fa.gz    # hap1_* and hap2_*, bgzip + fai + gzi
REGION=${REGION:-chr7:110000001-115000000}            # hg38 region (1-based, ends included); chosen with region_scan/ (2026-09-29): hap1 and hap2 each on a single contig at 100%, 0.61 SNP/kb between hap1 and hap2, no N
CTG=sim                                               # name of the simulated reference sequence
SEED=${SEED:-20260929}
# somatic variants (split in half between clonal and subclonal)
N_SNV=${N_SNV:-1000}; N_IND=${N_IND:-300}; N_HP=${N_HP:-0}   # short non-homopolymer indels (2-20 bp): 300 (2026-09-29); homopolymer indels: 0 (2026-09-29: at cap 2 they are invisible by construction and would mix the effect of HPC with that of EC; the loss is measured on HG008, C2)
N_DEL=${N_DEL:-12}; N_INS=${N_INS:-12}; N_INV=${N_INV:-8}; N_DUP=${N_DUP:-8}                  # SV: DEL 50 bp-5 kb, INS 50 bp-2 kb, INV and DUP 0.5-10 kb
# samples
DEPTHS="${DEPTHS:-90 60}"; PURITIES="${PURITIES:-1.0 0.6 0.3}"; NORMAL_DEPTH=35
EC_SETS="${EC_SETS:-fnr0.05 fnr0.10 fnr0.20 classic raw}"   # raw: tumor without EC (2026-09-30)
CLONALITY=0.5                                         # subclone fraction: the lowest VAF level to protect
NORMAL_HETCOV=4                                       # as in the production EC of the normal
READ_LEN_MEAN=17500; READ_LEN_SD=6000; PASSES=10      # HG008-T: 245 Gb in 14 M reads = 17.5 kb on average
YIELD=${YIELD:-1.03}                                  # depth requested from PBSIM3 = target x YIELD (ccs keeps ~97% with rq >= 0.99, test 2026-09-29)
UNTAG_FRAC=0.15                                       # fraction of reads without HP (real data: ~85% tagged)
WEIGHTS="3 5 8"                                       # m_min of the call with the normal (at 60x and rho 0.3 the subclonal variants have ~4.5 reads)
THREADS=${THREADS:-16}; PAR=${PAR:-3}
export PATH=$SV:$PATH
log(){ echo "[$(date '+%F %T')] $*"; }
samples(){ echo normal; for c in $DEPTHS; do for r in $PURITIES; do echo "t${c}_p${r}"; done; done; }
