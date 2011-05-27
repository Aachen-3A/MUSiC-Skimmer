#!/bin/bash

SCRIPTNAME=$(basename $0)

# Redirect stdout of echo to stderr
#
err() {
   echo "$SCRIPTNAME: ERROR: $@" 1>&2
}

abort() {
   echo "$SCRIPTNAME: Aborting..."
   exit $1
}

cleanUp() {
   killNC $1
   rm $FIFO
}

killNC() {
   # get the process id of nc by hand
   PID=$(pgrep -f "nc -l $1")
   if [[ -n "$PID" ]]; then
      kill $PID
   fi
}

wrongArg() {
   err "Malformatted argument '$*'"
   abort 2
}

# check the environment, we need CMSSW
#
if [[ -z "$CMSSW_BASE" ]]; then
   err "Environment 'CMSSW_BASE' not set"
   abort 1
fi

if [[ -z "$CMSSW_VERSION" ]]; then
   err "Environment 'CMSSW_VERSION' not set"
   abort 1
fi

if [[ -z "$SCRAM_ARCH" ]]; then
   err "Environment 'SCRAM_ARCH' not set"
   abort 1
fi

# parse the arguments
#
TEMP=`getopt -n $0 -o "f:kKl:p:r:s:" -l "firstRun:,kinit,perKeylastRun:,port:,rrurl:,server:" -- "$@"`

if [ $? != 0 ]; then
   err "Argument parsing failed!"
   abort 2
fi

eval set -- "$TEMP"

FIRST=1
LAST=999999
KINIT=false
KEY=false
SERVER="pccmsdqm04.cern.ch"
RRURL="http://$SERVER/runregistry/xmlrpc"
PORT=9999

while true; do
   case "$1" in
      -f|--firstRun)
         case "$2" in
            -*) wrongArg "$2";;
            *) FIRST="$2"; shift 2;;
         esac;;
      -k|--kinit) KINIT=true; shift;;
      -K|--perKey) KEY=true; shift;;
      -l|--lastRun)
         case "$2" in
            -*) wrongArg "$2";;
            *) LAST="$2"; shift 2;;
         esac;;
      -r|--rrurl)
         case "$2" in
            -*) wrongArg "$2";;
            *) RRURL="$2"; shift 2;;
         esac;;
      -s|--server)
         case "$2" in
            -*) wrongArg "$2";;
            *)  SERVER="$2"; shift 2;;
         esac;;
      -p|--port)
         case "$2" in
            -*) wrongArg "$2";;
            *)  PORT="$2"; shift 2;;
         esac;;
      --) shift;
         if [[ -n "$1" ]]; then
            err "Unknown option: '$1'"
            abort 3
         else
            break
         fi;;
      *) err "Unknown option: '$1'"
         abort 3 ;;
   esac
done

# SERVER must be a part of RRURL
#
if [[ ! "$RRURL" =~ "$SERVER" ]]; then
   err "Server '$SERVER' and rrurl '$RRURL' do not fit together!"
   abort 4
fi

################################################################################

# check if we have a ticket
# if the kinit option is set a new ticket is forced
#
if ! klist -5 -s || $KINIT; then
   kdestroy
   kinit -f $LOGNAME@CERN.CH -l 100d
fi

#exit

# create a temporary fifo
# it is needed to take care of the output sent to the server
#
FIFO=/tmp/fifo.hlt.$$
mkfifo $FIFO

echo "$SCRIPTNAME: Created fifo at '$FIFO'."

# some plumbing here:
# this starts a nc that listens to the specified port, pipes it to sed where
# the server name is correced, pipes the output through ssh to lxplus as the
# server is only accesible from there and sends it back to the fifo that puts
# the output back to the tcp port it came from
#
killNC $PORT
echo "$SCRIPTNAME: Using port $PORT for redirection."
< $FIFO nc -l $PORT | sed -ue "s/localhost:$PORT/$SERVER/" | ssh -x $LOGNAME@lxplus.cern.ch "nc $SERVER 80" | sed -ue 's%</methodResponse>%</methodResponse>\n%' > $FIFO &

GETHLT="$CMSSW_BASE/src/HLTrigger/Tools/python/getHLTkey.py"
if [[ ! -f $GETHLT ]]; then
   err "$GETHLT not found! Please make sure you have it in the right place!"
   cleanUp $PORT
   abort 5
fi
if $KEY; then
   $GETHLT --rrurl http://localhost:$PORT/runregistry/xmlrpc --perKey --firstRun $FIRST --lastRun $LAST
else
   $GETHLT --rrurl http://localhost:$PORT/runregistry/xmlrpc --firstRun $FIRST --lastRun $LAST
fi

cleanUp $PORT
