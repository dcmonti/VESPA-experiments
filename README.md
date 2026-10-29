# VESPA experiments

Scripts that produce the callsets and the evaluations of the VESPA paper.
VESPA itself is at https://github.com/dcmonti/VESPA; every experiment uses the release v1.0.1.

## Setup

1. Edit `config.sh`: `DATA` (root of inputs and outputs), the tool paths and the conda environments.
2. Create the conda environments in `envs/` (`conda env create -f envs/<name>.yml`).
3. Tools outside conda:
   - VESPA v1.0.1, static binary `vespa_linux_x86-64` of the GitHub release;
   - hifiasm v1.0.0-somatic, https://github.com/ldenti/hifiasm (tag `v1.0.0-somatic`);
   - Minda, https://github.com/KolmogorovLab/minda, commit `6e13ee5`;
   - DeepVariant 1.9.0, run with Singularity from `docker://google/deepvariant:1.9.0`.
4. Inputs: `bash download.sh` downloads every input to the path the scripts read, and checks it against its md5
   (`bash download.sh refs hg008 colo829 sim` for a subset). HG002 (experiment 10) downloads its own inputs.

## Published call sets

The benchmarks also score the published call sets below. They are not downloaded by `download.sh`; a benchmark
skips any of them that is missing.

| call set | source | path |
|---|---|---|
| HG008, Minda ensemble, HiFi: Severus, SAVANA, nanomonsv, Sniffles2 | GIAB FTP, `data_somatic/HG008/Liss_lab/analysis/NIH-NCI_minda-ensemble_20240710/vcf/{severus,savana,nanomonsv,sniffles}_hifi_vaf.vcf` | `$HG008/minda_ensemble/vcf/` |
| COLO829, Qin et al., v0.2 | Zenodo 10.5281/zenodo.17317401, `COLO829_hifi1.tar.gz` | extracted in `$COLO829/external_qinli_2026-09-24/v02/` |
| COLO829, Qin et al., minisv | as above, `COLO829_hifi1/minisv/msv_ltgs_COLO829_hifi1_somatic_generation5_filterasm.vcf`, converted with `python lib/msv2vcf.py <in> <dir>/minisv_generation5_filterasm.converted.vcf` | `$COLO829/external_qinli_2026-09-24/v02/COLO829_hifi1/minisv/` |
| COLO829, Qin et al., v4, SAVANA l+s | Zenodo 10.5281/zenodo.21852833, `ensembl__COLO829_hifi1__savana_l+s_3_filtered.vcf` | `$COLO829/external_qinli_v4_2026-10-08/` |
| COLO829, Keskus et al., Severus 1.2 | Zenodo 10.5281/zenodo.14541057, `variant_calls_and_benchmarks.tar.gz` (15 GB), file `variant_calls_and_benchmarks/COLO829/original_vcfs/severus_pb.vcf.gz` | extracted in `$COLO829/external_severus_2026-10-08/extract/` |

## Experiments

| | experiment | needs | script |
|---|---|---|---|
| 1 | HG008, comparison of the callers, Minda and truvari, clonal and clonal + subclonal regions | — | `run_exp_1.sh` |
| 2 | COLO829, comparison of the callers, minimap2 default (`-s200`) | — | `run_exp_2.sh` |
| 3 | COLO829, comparison of the callers, minimap2 `-s50` | 2 | `run_exp_3.sh` |
| 4 | effect of each germline filter | 1, 3 | `run_exp_4.sh` |
| 5 | HG008 without homopolymer compression and without the SFS length cap | 1 | `run_exp_5.sh` |
| 6 | HG008, somatic signal discovery with SFS: breakpoint coverage, germline cascade, SFS classes | 1 | `run_exp_6.sh` |
| 7 | VAF of the ground-truth events (svviz2) | 1, 3 | `run_exp_7.sh` |
| 8 | error correction on simulated tumor/normal pairs | — | `run_exp_8.sh` |
| 9 | over-dispersion of the read depth of the normal | 1, 2 | `run_exp_9.sh` |
| 10 | HG002, per-read detection rate q of LOWPOWER | — | `run_exp_10.sh` |

`run_exp_all.sh` runs them all in this order.
Each step skips the outputs it already finds, so a failed run resumes from the failed step when relaunched.
Outputs go under the directory of the dataset (`$HG008/reads_ec_align`, `$COLO829/reads_ec_align{,_s50}`, `$HG002/q_visibility`).

## Layout

- `config.sh`: paths and tools, the only file to edit.
- `download.sh`: the inputs; `data/`: the centromere and telomere BED, which has no public source.
- `run_exp_N.sh`, `run_exp_all.sh`: one script per experiment, and all of them.
- `reproduce/<experiment>/`: the steps of each experiment; the header of each script states its inputs, steps and outputs.
- `reproduce/bench/`: VCF preparation, truvari event counting and the excluded-truth rule, shared by every benchmark.
- `lib/`: helpers shared by the experiments (SFS deduplication, haplotag transfer, SFS projection on the genome).
- `envs/`: conda environments.
