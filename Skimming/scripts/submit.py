#!/usr/bin/env python

import sys
import subprocess
import imp
import pickle
import optparse
import datetime
import ConfigParser


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
    config = ConfigParser.RawConfigParser()
    config.add_section( 'CRAB' )
    config.set( 'CRAB', 'jobtype', 'cmssw' )
    config.set( 'CRAB', 'scheduler', 'glite' )
    config.add_section( 'CMSSW' )
    config.set( 'CMSSW', 'datasetpath', sample )
    config.set( 'CMSSW', 'pset', name+'_cfg.py' )
    config.set( 'CMSSW', 'total_number_of_events', '-1' )
    config.set( 'CMSSW', 'events_per_job', '50000' )
    if options.runs:
        config.set( 'CMSSW', 'runselection', options.runs )
    config.set( 'CMSSW', 'output_file', name+'.pxlio' )
    config.add_section( 'USER' )
    config.set( 'USER', 'return_data', '0' )
    config.set( 'USER', 'copy_data', '1' )
    config.set( 'USER', 'storage_element', 'T2_DE_RWTH' )
    config.set( 'USER', 'user_remote_dir', outname+name )
    config.add_section( 'GRID' )
    config.set( 'GRID', 'rb', 'CERN' )
    config.set( 'GRID', 'group', 'dcms' )
    config.set( 'GRID', 'se_black_list', 'T0,T1' )
    config.set( 'GRID', 'additional_jdl_parameters', 'rank=-other.GlueCEStateEstimatedResponseTime+(RegExp("rwth-aachen.de",other.GlueCEUniqueID)?100000:0)+(RegExp("desy.de",other.GlueCEUniqueID)?100000:0)' )

    cfg_file = open(name+'.cfg', 'wb')
    config.write( cfg_file )

 
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
