## This is a bootstrap script for the MUSiC Skimmer
# Please source it before using the Skimmer, for him
# to work properly
 
 
# first get th dir where the bootstrap script is placed
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

    
if [ -z ${CMSSW_BASE+x} ]; then 
    echo "CMSSW_BASE is unset!"; 
    echo "Please source CMSSW before bootstraping the skimmer"; 
else 
    echo "CMSSW_BASE is set to '$CMSSW_BASE'"; 
    ln -fs  $DIR $CMSSW_BASE/src/MUSiC-Skimmer
fi
