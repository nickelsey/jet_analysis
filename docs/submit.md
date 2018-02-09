# Submit scripts, file lists, and data quality lists

## Overview
This subproject contains python submit scripts that are used to submit my analyses on the WSU grid. I use python instead
of a bash script because python, with the root library, allows me to check if the job was successful in producing output. This
became importing for run 14 data, where we have thousands of data files to run over.

## Data lists

### pp_list
a list of all our y6 ppHT data files, both on rhic121, and on the WSU grid (in pp_list/grid)

### y14_mb_file_list.txt
This list is used for embedding when running AJ on pp.

### tower lists
Bad tower lists for y6 pp, y7 AuAu, y14 AuAu, and combinations of y6 & y7 or y14. The combination lists should be used when
trying to do comparisons between pp & AuAu (like for AJ) so that the barrel calorimeter has essentially the same physical acceptance
for both pp and AuAu.
y6_bad_tower.txt
y7_y6_bad_tower.txt
y14_bad_tower.txt
y14_y6_bad_tower.txt
empty_list.txt (for testing, or when running data qa)

### bad run lists
These are lists generated from data QA of bad runs - either the detector had problems, there was a problem with RHIC, etc. QA is done for
vertex position, refmult, tracks, towers, etc, and runs that have values significantly different than the average are rejected.
y14_bad_run.txt

### qwrap.sh
This file is used during job submission to pbs, and is handled by the pbs submit scripts

### run ID files
My data QA needs the run IDs of all runs before analysis, so these lists must be pre-computed. These are trees of run IDs for each run period.
y14_runids.root
y7_runids.root

## Submit scripts
These are python scripts used to submit & monitor jobs on the WSU grid.

### pbs_submit_basic_dijet_imbalance.py
This is kept around mainly as legacy - I use the next routine instead.

### pbs_submit_dijet_imbalance.py
This file submits jobs to run the dijet imbalance measurement. It can submit Au+Au, p+p, and p+p embedded in Au+Au. Efficiency corrections
can be turned on as well, either for intra-system (i.e. Au+Au central vs peripheral) or inter-system (Au+Au vs p+p). The following are some of
the most important command line options.
```
--output              the output directory to store output root files
--triggerEfficiency   turn on relative efficiency corrections for triggered data
                      (Au+Au wrt centrality, or pp wrt Au+Au, etc)
--embedEfficiency     turn on efficiency corrections for embedding data
                      (Au+Au wrt centrality)
--pp                  specify trigger data is p+p
--no_pp               specify trigger data is not p+p (it is Au+Au)
--auau                specify trigger is Au+Au (it is not p+p)
--forcePPEffCent      if p+p is embedded, forces the efficiency corrections to be wrt
                      this centrality bin. If not, then it is selected from embedding
                      centrality, and if no centrality is selected, it is chosen randomly
--badRuns             CSV of runs to reject in both triggered and embedding data
--badTowers           CSV of tower IDs to reject in both triggered and embedding data
--triggers            string specifying which triggers to use for triggered data (e.g. y14ht)
--embedTriggers       same but for embedding data
--embed               text file with all embedding data files listed inside. It will select one
                      randomly for each trigger data file
--reuseTrigger        if using p+p, tell the analysis to use the same event multiple times, to
                      increase statistics
--readerSetting       can specify a text file with reader settings if you do not want to use the
                      defaults. Example is given in the test directory
--embedReaderSetting  same for the embedding data

ANALYSIS VARIABLES
--constEta            list of maximum eta values to use for constituents
--leadConstPt         list of leading jet constituent pt cuts
--leadConstPtMatch    list of leading matched jet constituent pt cuts
--subConstPt          list of subleading jet constituent pt cuts
--subConstPtMatch     list of subleading matched jet consttituent pt cuts
--leadR               list of leading hard jet radii
--leadJetPt           list of leading hard jet pt cuts
--subR                list of subleading hard jet radii
--subJetPt            list of subleading hard jet pt cuts

```
### pbs_submit_check_refmult.py
This script is used to generate refmult distributions with different track quality cuts. It was used mainly to make sure we knew
how refmult was calculated originally in the muDST files. Not in current use. Important command line options given here:
```
--output              the output directory to store output root files
--badRuns             CSV of runs to reject in both triggered and embedding data
--badTowers           CSV of tower IDs to reject in both triggered and embedding data
--triggers            string specifying which triggers to use for triggered data (e.g. y14ht)
--beginRun            accept only runs > beginRun
--endRun              accept only runs < endRun
```

### pbs_submit_comprehensive_dq.py
This script submits jobs that do relatively comprehensive data quality analysis on the input data. It produces event level
histograms such as refmult, grefmult, nglobal, zdc rate, vz, vx, vy, etc. Also analyzes tracks and towers.
REQUIREMENT: needs a Tfile with a ttree of all run IDs present in the data set PRIOR to running comprehensive_dq
```
--output              the output directory to store output root files
--badRuns             CSV of runs to reject in both triggered and embedding data
--badTowers           CSV of tower IDs to reject in both triggered and embedding data
--triggers            string specifying which triggers to use for triggered data (e.g. y14ht)
--runIDs              a TFile with a TTree named runid, with a branch of unsigned ints named runid
                      that contains all the runIDs present in the data set of interest
--trigger             which triggers to accept: discard all other events
```

