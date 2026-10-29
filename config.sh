# Paths and tools of every experiment: the only file to edit.
# Every script of reproduce/ sources it; all variables are exported, so the Python scripts read them too.
set -a

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LIB=$REPO/lib                 # helper scripts shared by the experiments
REPRO=$REPO/reproduce

# Data: inputs and outputs. Each experiment writes under the directory of its dataset.
DATA=/projects6/luca-nobackup/somatic
REFS=$DATA/refs               # GRCh38-GIABv3, centromere/telomere BED, tandem-repeat BED
HG008=$DATA/hg008/analysis_dcm
COLO829=$DATA/colo829
HG002=$DATA/hg002

# Our tools
VESPA=$DATA/tools/vespa/v1.0.1/vespa_linux_x86-64      # VESPA v1.0.1, static build, md5 1aa531662f2afccb587442a40aa0f615
HIFIASM=$DATA/tools/hifiasm/v1.0.0-somatic/hifiasm     # v1.0.0-somatic (hifiasm 0.25.0-r726), sha256 ad32741d...399d2

# Conda environments (envs/)
CONDA=/home/dmonti/miniconda3/bin/conda
CONDA_SH=/home/dmonti/miniconda3/etc/profile.d/conda.sh
ENVS=/home/dmonti/miniconda3/envs
SV_BIN=$ENVS/sv_somatic/bin   # samtools 1.23.1, bcftools, tabix, bgzip, minimap2 2.30-r1287, whatshap 2.8, truvari 5.4.0,
                              # python with pysam 0.24.0, PyVCF, numpy 2.4.6, scipy 1.17.1; Severus 1.7
PYTHON=$SV_BIN/python
SAMTOOLS=$SV_BIN/samtools
BCFTOOLS=$SV_BIN/bcftools
TABIX=$SV_BIN/tabix
MINIMAP2=$SV_BIN/minimap2
WHATSHAP=$SV_BIN/whatshap
SEVERUS=$SV_BIN/severus                         # 1.7
SAVANA=$ENVS/savana138/bin/savana               # 1.3.8
NANOMONSV=$ENVS/nanomonsv090/bin/nanomonsv      # 0.9.0
SNIFFLES=$ENVS/sniffles281/bin/sniffles         # 2.8.1
MINDA=$DATA/tools/minda/minda.py                # KolmogorovLab/minda, commit 6e13ee5
MINDA_PYTHON=$ENVS/minda/bin/python

set +a
