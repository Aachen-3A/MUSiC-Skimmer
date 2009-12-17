#!/usr/bin/env python

import sys
import subprocess
import imp
import pickle
import optparse
import datetime


parser = optparse.OptionParser( description='Submit MUSiCSkimmer jobs, using CMSSW config CFG_FILE, on all samples listed in DATASET_FILE',  usage='usage: %prog [options] CFG_FILE DATASET_FILE' )
parser.add_option( '-n', '--name', metavar='NAME', help='Output will be written in store/user/{your_name}/NAME/{short_dataset_name} [default: MUSiC/{current_date}]' )
parser.add_option( '-r', '--runs', metavar='RUNS', help='Only analyze the given runs' )

(options, args ) = parser.parse_args()

if len( args ) < 2:
    parser.error( 'CFG_FILE and DATASET_FILE required' )

del parser

pset = args[0]
samples = args[1]

if options.name:
    outname = options.name
    while outname.startswith( '/' ):
        outname = outname[ 1: ]
    while outname.endswith( '/' ):
        outname = outname[ :-1 ]
else:
    outname = 'MUSiC/'+datetime.date.today().isoformat()
outname += '/'

if options.runs:
    run_line = 'runselection = '+options.runs+'\n'
else:
    run_line = ''


print 'Reading config', pset
file = open( pset )
cfo = imp.load_source("pycfg", pset, file )
del file
process = cfo.process
del cfo


for line in open( samples ):
    line = line[:-1]
    if not line or line.startswith( '#' ): continue
    line = line.split( ':' )
    name = line[0]
    sample = line[1]
    print name, ': Generating CRAB cfg...',
    cfg = ( '[CRAB]\n'+
            'jobtype = cmssw\n'+
            'scheduler = glite\n'+
            '[CMSSW]\n'+
            'datasetpath = '+sample+'\n'+
            'pset = '+name+'_cfg.py\n'+
            'total_number_of_events=-1\n'+
            'events_per_job = 50000\n'+
            run_line+
            'output_file = '+name+'.pxlio\n'+
            '[USER]\n'+
            'return_data = 0\n'+
            'copy_data = 1\n'+
            'storage_element = T2_DE_RWTH\n'+
            'user_remote_dir = '+outname+name+'\n'+
            '[GRID]\n'+
            'rb = CERN\n'+
            'group = dcms\n'+
            'se_black_list = T0,T1\n'+
            '[CONDORG]\n'
            )
    cfg_file = open(name+'.cfg', 'w')
    cfg_file.write(cfg)
    cfg_file.close()

    print 'Generating CMSSW config...',
    process.Skimmer.FileName = name+'.pxlio'
    process.Skimmer.Process = name

    pset_file = open( name+'_cfg.py', 'w' )
    pset_file.write( "import FWCore.ParameterSet.Config as cms\n" )
    pset_file.write( "import pickle\n" )
    pset_file.write( "pickledCfg=\"\"\"%s\"\"\"\n" % pickle.dumps( process ) )
    pset_file.write("process = pickle.loads(pickledCfg)\n")
    pset_file.close()

    print 'done and submitting...'
    subprocess.call( [ 'crab', '-create', '-submit', '-cfg', name+'.cfg' ] )
