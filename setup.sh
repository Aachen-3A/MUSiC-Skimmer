# !/usr/bin/bash

if [ -z "${CMSSW_BASE+x}" ]; then
    echo "CMSSW_BASE is not set!"
    echo "Source CMSSW before setting up the skimmer."
    exit 1
fi

bold=$(tput bold)
normal=$(tput sgr0)

echo "CMSSW_BASE is set to '${CMSSW_BASE}'"
echo ""
echo "${bold}Setting up the skimmer ..."
echo ""
cd ${CMSSW_BASE}/src

# Add the recommended electron idenfitication source files
git cms-merge-topic 9003 #this is the version that is in CMSSW_7_4_X
rm -rf RecoEgamma/ElectronIdentification/data
git clone https://github.com/cms-data/RecoEgamma-ElectronIdentification.git RecoEgamma/ElectronIdentification/data
rm -rf RecoEgamma/PhotonIdentification/data
git clone https://github.com/cms-data/RecoEgamma-PhotonIdentification.git RecoEgamma/PhotonIdentification/data

# Pxl installation
hg clone https://forge.physik.rwth-aachen.de/hg/cmssw-modules/Pxl
cd Pxl/
hg up pxl-3.5.1
cd ..

echo ""
echo "${bold}Done."
