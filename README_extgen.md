# NuWro request mode for MicroBooNE ExtGen (`bin/nuwro_extgen`)

Branch `v08_00_00_br`: NuWro 25.11.1 (tag `NuWro_25_11_1`) with the Σ⁻ PDG
fix of NuWro PR #61 and a request-mode program for MicroBooNE's ExtGen
(`ubsim/ExtGen`). Every change is marked `[uboone patch]`.

## What ExtGen needs from a generator

ExtGen lets GENIE choose each neutrino interaction's flavour, energy, target
nucleus and current (CC/NC) from the real flux files and detector geometry
(splines made from the external generator's own cross sections), and asks the
external generator for the final state of each of those interactions, called
*requests*. The final states are merged back into GENIE's events (vertex,
time, POT stay GENIE's) and overlaid on data as usual.

## The program

`bin/nuwro_extgen` is a separate executable. Stock `bin/nuwro` and NuWro's
physics are unchanged: the program calls NuWro's own `makeevent` and
`finishevent` (FSI) through a `friend` declaration in `src/nuwro.h`.

One process serves one *group*: a target nucleus, a neutrino flavour and a
current, fixed by the card and `-p` overrides (`nucleus_p`, `nucleus_n`,
`beam_particle`, the `dyn_*` flags of that current). Only the energy changes
from request to request (a fixed-energy beam along +z).

```
nuwro_extgen -m serve -r requests.txt -s served_events.txt -t trial_stats.txt \
             -i card.txt -o PREFIX -p random_seed=N [-p param=value ...]
nuwro_extgen -m scan -n TRIALS -r energies.txt -t scan.txt \
             -i card.txt -o PREFIX -p random_seed=N [-p param=value ...]
```

`-i`, `-o`, `-p` are NuWro's own options (`PREFIX.par` lists all parameters).
`random_seed` must be > 1; Pythia6 (DIS hadronization), which stock NuWro
never seeds, is seeded from it as well.

### Serve mode

Request file:
```
# extgen-nuwro-requests v1 channels=0,2,4,6,8,10
req_id E_GeV M_1 ... M_K
```
`channels` must equal the process's active channels (NuWro `dyn` codes:
0/1 QEL, 2/3 RES, 4/5 DIS, 6/7 COH, 8/9 MEC, 10 HYP; even CC, odd NC); `M_k`
is the largest weight allowed in channel k at that energy (1e-38 cm² per
nucleon; from the scan campaign with a margin).

For each request: draw channel k with probability M_k / ΣM, make a NuWro event
of channel k (weight w), accept it with probability w / M_k, else start again
(channel included). The accepted events follow the summed cross section
exactly, channel fractions included, as long as w ≤ M_k; a larger weight
stops the program with exit code 3 (`OVERWEIGHT ...`), so that the caller
retries with larger Max values. A channel with M_k = 0 is never generated.
Unlike stock NuWro's `real_events` there is no adaptive maximum
(`chooser::accept`), no rounding of per-channel event counts, and no DIS
energy bias.

Outputs (ExtGen formats, `ubsim/ExtGen/Modules/ExtGenIO.h`):

* `served_events.txt`: per request
  `EVENT req_id E_GeV gen_proc simb_mode simb_inttype nparts`, then
  `P status pdg px py pz E` (GeV, target rest frame, neutrino along +z):
  * status 11, the struck nucleon: NuWro's `in[1]` (QEL, RES, DIS, HYP);
    MEC 2p2h: the initial pair as one cluster, PDG 2000000200 + number of
    protons (as GENIE); none for COH and for MEC 3p3h;
  * status 1: every particle leaving the nucleus (`post`), the primary
    lepton (charged lepton for CC, the neutrino for NC) first;
  * `gen_proc` = `dyn`, + 100 for MEC 3p3h; `simb_mode`: QEL and HYP 0 (QE),
    RES 1, DIS 2, COH 3, MEC 10; `simb_inttype` as GENIE2ART (1001/1002 QE,
    1091/1092 DIS, 1097/1096 COH, 1000 RES and MEC).
* `trial_stats.txt`: `req_id tgt_pdg nu_pdg ccnc E_GeV n_trials sum_u sum_u2
  max_u n_dropped`, with the trial weight u = w·ΣM/M_k per nucleus (1e-38 cm²):
  its mean is the group's total cross section, which ExtGen's in-job check
  (V4) compares with the spline table. `n_dropped`: accepted events without
  their primary lepton after FSI (none seen so far).

NaN weights (numerical failures inside NuWro's models; seen in the hybrid RES
model for antineutrinos above ~1.5 GeV, ~0.3% of those trials) are discarded
as stock NuWro's chooser does: such a trial counts in no statistic, so cross
sections are means over valid trials, as NuWro reports them. Their number is
reported.

### Scan mode

Request lines `req_id E_GeV` (further columns ignored). For each request and
each active channel, TRIALS trials without FSI; output lines
`req_id E_GeV dyn n_trials sum_w sum_w2 max_w n_nan` (per nucleon, 1e-38
cm²; `n_trials` valid trials, `n_nan` discarded NaN trials on top): σ_k(E) for
the spline tables and the largest weights for the Max tables.

## Card requirements

* `target_type = 0` (one nucleus), `beam_test_only = 0`, `dyn_lep = 0`,
  electron channels off, the channels of one current only.
* **Spectral functions:** QEL with `sf_method > 0` exists for C12, O16, Ar40
  and Fe56 only; for every other nucleus `sfevent` returns weight 0 without a
  warning. `nuwro_extgen` refuses that combination: other nuclei (and free
  nucleons) need `sf_method = 0` (local Fermi gas). On argon, LFG gives a QEL
  cross section ~20 % above the spectral function (1 GeV ν_μ CC: 0.57 vs 0.47
  ×1e-38 cm² per nucleon), so a sample mixes models by target.

## Build

As stock NuWro (ROOT + ROOTEGPythia6): `make` builds `bin/nuwro_extgen` with
the other targets (or `make bin/nuwro_extgen`); CMake: target `nuwro_extgen`.

## Tests (local, ROOT 6.32 + ROOTEGPythia6, 2026-10-25)

* Scan vs stock `nuwro` test events (Ar40 ν_μ CC, 1 GeV, MicroBooNE card):
  QEL 0.470 ± 0.003 vs 0.465, RES 0.326 ± 0.003 vs 0.336, MEC 0.1188 vs
  0.1193, COH 0.0061 vs 0.0061 (×1e-38 cm² per nucleon).
* Serve, 20000 requests at 1 GeV: channel fractions QEL 0.509, RES 0.356,
  COH 0.006, MEC 0.130 against the scan's 0.510, 0.354, 0.007, 0.129; ~35
  trials per event, ~12 ms per event.
* Groups: ν_μ NC, ν̄_μ CC (hyperons, MEC clusters np/pp, 3p3h), ν_e CC on Ar40;
  C12, O16, Fe56 (SF); Na23, Ti48, Ni59, free p and n (LFG); overweight stop.
