# Table of Content
1. [Introduction](#introduction)
2. [Setup](#setup)



# Introduction <a name="introduction"></a>
The scripts used for QmisID rate estimation can be found under fold `Scripts/`. 
- To get the MC Truth-matched rate, run `QmisID/Scripts/ChargeFlipTruth.C`
- To get the Likelihood rate, run `QmisID/Scripts/Likelihood.C`

To get some inspirations, see examples commands under fold `run/`
- `QmisID/run/run_commands_ChargeFlipTruth.sh`
- `QmisID/run/run_commands_Likelihood.sh`


*To be added*: 
1. scripts for closure test and systematics

# Setup <a name="setup"></a>
For the first time:
```
setupATLAS
lsetup git
git clone ssh://git@gitlab.cern.ch:7999/shuhui/QmisID.git
cd QmisID/
source setup.sh
```

When starting a fresh terminal, you need to run again:
```
source setup.sh
```
