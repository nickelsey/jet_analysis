# Data quality analysis & helper routines

## overview
This subproject is used primarily as a common repo for some of the routines I have used when analyzing a new data set.


### check_refmult
This routine was used only once for run 14, to make sure that the reported refmult in the header matched with our expected definition of refmult -
specifically, |eta| < 0.5, nHitsFit >= 10, DCA < 3.0. The results confirmed that refmult was calculated properly. The reason for the recalculation was
that the refmult distribution for run 14 looked different from year 7, but there were changes in tracking, TPC efficiency, and luminosity, which more than
likely fully explain the difference.

This routine can be submitted via the following command. Non-default settings can be set using command line flags, they are listed in the python script.
```
python submit/pbs_submit_check_refmult.py /path/to/data/files/*

```

### comprehensive_data_quality
Relatively comprehensive data quality of the input data. Requires a pre-made list of all runIDs of interest in the run period in a TTree.
Makes histograms of vertex position, zdc rate, refmult, track kinematics, tower kinematics, and others. TH2Ds are created with many variables
plotted as a function of time (runid). Many variables that are expected to be correlated are also saved in TH2Ds.
Actual quality analysis is not done in this script - it is a little misleading. Here, I only do the histogram filling. The actual analysis is done in routines
such as print_data_quality.cc, for instance

This routine can be submitted via a python script
```
python submit/pbs_submit_comprehensive_dq.py /path/to/data/files/*
```

### print_data_quality
This script does the analysis of the histograms produced in comprehensive_data_quality. It takes relevant projections and profiles. However,
this routine is mainly to compare y14 to y7. A more generic printing & analysis routine will need to be made.

This routine can be run like this:
```
./bin/data_qa/print_data_quality --output=output/location --input=hadded/root/file 
```
