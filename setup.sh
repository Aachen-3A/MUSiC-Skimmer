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
svn export https://github.com/ikrav/cmssw/branches/egm_id_74X_v1/RecoEgamma/ElectronIdentification RecoEgamma/ElectronIdentification
svn export https://github.com/ikrav/cmssw/branches/egm_id_74X_v1/RecoEgamma/EgammaTools RecoEgamma/EgammaTools
svn export https://github.com/ikrav/cmssw/branches/egm_id_74X_v1/PhysicsTools/SelectorUtils PhysicsTools/SelectorUtils

# Pxl installation
hg clone https://forge.physik.rwth-aachen.de/hg/cmssw-modules/Pxl
cd Pxl/
hg up pxl-3.5.1
cd ..

echo ""
echo "${bold}Done."
