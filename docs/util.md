# Utility functions and helper classes

## Overview
These are helper functions and quality of life improvements.

### arg_helper
Implements a set of functions that help set up command line option parsing. Look at dijet imbalance analysis for example
of how it works.

### dijet_key
The keys in the dictionary used by the dijet_worker have to be quite long to avoid clashes. This function handles creating the
keys from the input dijet_definition.

### histogram_routines
Simple routines to calculate mean & standard deviation in histograms, given constraints along one axis. (For instance, calculating
the mean of the X axis, but ignoring bins without a certain number of entries)

### pt_distribution
An implemention of our normal pt*exp(-pt/T) probability distribution of pT of the thermal background. This
works with <random> generators and conforms to <random>'s distribution interface. However, it does have to do an MC integration
so its not as fast as the analytical functions that <random> has.

### reader_util
Functions to make initializing TStarJetPicoReader more convenient and less error prone.

### root_print_routines
Histogram and profile printing routines used for overlaying series of histograms, or just beautifying the ugly root defaults

### selector_compare
It is impossible to get a fastjet::selector's numerical cuts without some hacking, so I parse them from the selector's info string.
This allows me to compare selectors numerically (and check that they implement similar cuts)

### trigger_lookup
A lookup table of all triggers we would be interested in, including AuAu, pp and some dAu, from y6, y7, y10, y11,  y12, y14, etc.
Includes high tower, jet patch, min bias, etc.

### vector_conversion
Functions to convert from TStarJetPico or Pythia to fastjet::pseudojet.
