# PxlSkimmer
## Quick Installation

```bash
# Setup and source CMSSW
export SCRAM_ARCH=slc6_amd64_gcc491
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_7_4_6_patch2
cd CMSSW_7_4_6_patch2/src
cmsenv

# Miscellaneous setup
git cms-addpkg PhysicsTools/PatAlgos
git clone git@github.com:Aachen-3A/PxlSkimmer.git
./PxlSkimmer/setup.sh

# Compile the code
scram build -j 6
```

## Installation
Before installing the PxlSkimmer, one has to setup CMSSW.

When CMSSW is set up and sourced, one can either clone the PxlSkimmer repository
into `${CMSSW_BASE}/src` or add a CMS package to set up the CMSSW git
infrastructure. As a recommended package to install the CMSSW git
infrastructure, the `PatAlgos` can be added by calling `git cms-addpkg
PhysicsTools/PatAlgos` in the `${CMSSW_BASE}/src` directory. To clone the
PxlSkimmer one has to call `git clone git@github.com:Aachen-3A/PxlSkimmer.git`

After cloning the PxlSkimmer, one has to execute the `setup.sh` script inside the
PxlSkimmer directory, which will take care of the remaining additions required
to run the skimmer.

## Using the Skimmer on the GRID
The submission of Jobs to the GRID can be done using the script music_crab3, see https://github.com/Aachen-3A/tools3a/wiki/music_crab3/
# Coding conventions
TODO
