# Dijet Imbalance

## Overview
Implementation of the Dijet Imbalance measurement for Au+Au or p+p STAR data, along with helper functionality.

### basic_dijet_imbalance
Depricated

### differential_aj
This is currently the active AJ analysis code. It implements event, track and tower cuts, then converts TStarJetPicoEvents
into fastjet::PseudoJets, then passes those PseudoJets to a DijetWorker, which implements the actual clustering & matching.
The DijetWorker is set up to handle multiple dijet definitions, and these definitions are constructed from command line arguments.
To see how to submit this on the grid, look in docs/submit.md, and look for pbs_submit_dijet_imbalance.py

### print_auau_pp_results
Printing routines to beautify and overlay the results from Au+Au & p+p if applicable. Can specify input root files vias --auau &
--pp, and an output directory for printing via  --outputDir
