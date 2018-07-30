# DijetWorker

## Overview
The DijetWorker class implements the dijet imbalance measurement, given a DijetMatrix & a set of fastjet::PseudoJets.
The DijetMatrix class defines, given a set of constituent & jet level kinematics, all possible permutations of those sets.
(However it discards permutations that don't make sense. For instance, if the set has subleading jet pt > leading jet pt).

### JetDef
Contains constituent cuts, R, jet algorithm, background subtraction options, etc. Essentially, just a more complete fastjet::JetDefinition.

### MatchDef
Two JetDefs make a MatchDef: one "hard" jet definition and one matching jet definition. However JetDef/MatchDef combo is
robust enough that the initial JetDef does not need to be "hard" and matching can be disabled, if not needed.

### DijetDefinition
Two MatchDefs make a DijetDefinition. STAR's AJ dijet definition involves a two steps, One finding a hard dijet pair, and one step matching
full jets to the hard jets. A DijetDefinition represents a more generic version of two-step di-jet analyses, where, for instance, the initial phase
does not have to be "hard", or the matching can be disabled, etc.

### DijetMatrix
The DijetMatrix is the class that, given a set of constituent and jet cuts, creates all valid permutations of those cuts, and creates a
DijetDefinition for each.

### DijetWorker
The DijetWorker inherits from the DijetMatrix. It would be more apt to call it the DijetAjWorker, as it specifically implements the AJ analysis.
Given a set of fastjet::PseudoJets, the DijetWorker will run the AJ analysis, once for each DijetDefinition that the DijetMatrix has defined.
The results of all clusterings are stored in a custom data structure that is returned in a dictionary, and easily accessed by the
user via unique identifier strings.
