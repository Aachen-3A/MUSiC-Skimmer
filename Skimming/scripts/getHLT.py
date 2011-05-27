#!/usr/bin/env python

import os
import sys
major, minor, micro, releaselevel, serial = sys.version_info
if major < 2 or minor < 4:
   print 'FATAL: Requires Python >= 2.4...Activate a CMSSW to get that.'
   sys.exit(1)
cmssw_version = os.getenv( "CMSSW_VERSION" )
if cmssw_version == None:
   print 'FATAL: Environment CMSSW_VERSION not set'
   sys.exit(1)
cmssw_base = os.getenv( "CMSSW_BASE" )
if cmssw_base == None:
   print 'FATAL: Environment CMSSW_BASE not set'
   sys.exit(1)
scram_arch = os.getenv( "SCRAM_ARCH" )
if scram_arch == None:
   print 'FATAL: Environment SCRAM_ARCH not set'
   sys.exit(1)

import optparse
import simplejson as json
import subprocess

usage = """Usage: %prog [options] DATASETS
           DATASETS should be a space separated list of datasets you are
           interested in.
           If DATASETS is omitted, the default datasets are used:
           Jet, METBTag, Photon, SingleElectron, SingleMu"""
parser = optparse.OptionParser( usage = usage, version = '%prog version 1' )
parser.add_option( '-v', '--verbose', action = 'store_const', const = 1, dest = 'verbose', help = 'Tell the script to be verbose, i.e. to dump warnings.' )
parser.add_option( '-d', '--dcs-file', metavar = 'DSC', default = os.path.join( cmssw_base, 'src/MUSiCProject/Skimming/test/lumi/json_DCSONLY.txt' ),
                   help = 'Define the dcs file. [default: $CMSSW_BASE/src/MUSiCProject/Skimming/test/lumi/json_DCSONLY.txt]' )
parser.add_option( '-f', '--firstRun', metavar = 'FIRST', default = 1, help = 'Define the first run to be analysed. [default: %default]' )
parser.add_option( '-l', '--lastRun', metavar = 'LAST', default = 999999, help = 'Define the last run to be analysed. [default: %default]' )
parser.add_option( '-F', '--filter', metavar = 'FILTER', action = 'append', default = None,
                   help = """Define filter words. Trigger bit names must contain FILTER to be considered.
                             If you use ~FILTER the filter is negated and the trigger bit name must not contain FILTER to be considered. [default: %default]""" )
parser.add_option(       '--debug', action = 'store_true', default = False, help = 'Activate additional output for debugging. [default: %default]' )
parser.add_option( '-o', '--output', metavar = 'OUTPUT', help = 'Define the output file for the trigger ranges. [default: HLT_<firstRun>-<lastRun>.txt]' )

( options, args ) = parser.parse_args()

# take care of dynamic output file name
#
if not options.output:
   options.output = 'HLT_' + options.firstRun + '-' + options.lastRun + '.txt'

options.neg_filters = None
if options.filter:
   options.neg_filters = list()
   for filter in options.filter[:]:
      if filter.startswith( '~' ):
         options.neg_filters.append( filter[1:] )
         options.filter.remove( filter )


if options.debug:
   print 'Using filters:', options.filter
   print 'Using negated filters:', options.neg_filters
   print
   print 'Using DCS JSON file:', options.dcs_file

options.fullCall = ' '.join( sys.argv )

options.datasets = [ 'Jet', 'MET', 'METBTag', 'Photon', 'SingleElectron', 'SingleMu' ]

if args:
   # remove duplicates by converting into a set and back again
   #
   options.datasets = list( set( args ) )


class Run:
   num = None
   config = None

   def __str__( self ):
      return ( 'Run: ' +
               'num = '    + str( self.num )    + '; ' +
               'config = ' + str( self.config )
               )

   def __repr__( self ):
      return self.__str__()

   def __cmp__( self, alien ):
      if self.num < alien.num:
         return -1
      elif self.num > alien.num:
         return 1
      else:
         return 0


class Range:
   start_run = None
   end_run = None
   unprescaled_triggers = None

   def __str__( self ):
      return ( 'Run range: ' + str( self.start_run ) + ' - ' + str( self.end_run ) + '; ' +
               'unprescaled triggers: ' + str( self.unprescaled_triggers )
               )

   def __repr__( self ):
      return self.__str__()


# get runs that are specified in the DCS only file
#
def getDCSRunList( dcs_file ):
   dcs = json.load( file( dcs_file ) )

   runs = [ int( run ) for run in dcs ]

   return runs


# get the runs in the specified run range sorted by their HLT config
#
def getHLTRunList( options ):
   getHLTKeys = os.path.expandvars( '$CMSSW_BASE/src/MUSiCProject/Skimming/scripts/HLTKeys.sh' )

   proc = subprocess.Popen( [ getHLTKeys, '-f', str( options.firstRun ), '-l', str( options.lastRun ), '-K' ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
   #for line in proc.stdout: #.readline():
      #print line
   output = proc.communicate()[0]
   
   if options.debug:
      print
      print "From '%s'" % getHLTKeys
      print output
      print "End from '%s'" % getHLTKeys
      print

   raw_lines = output.splitlines()

   # filter HI and non physics runs and lines that have nothing to do with HLT
   #
   hlt_runs = []
   for line in raw_lines:
      if 'cdaq' in line:
         if 'HI' not in line:
            if 'physics' in line:
               # the first part is the HLT config
               # the rest are run numbers with this config
               #
               config = line.split()[ 0 ]
               runs = line.split()[ 1: ]
               for run in runs:
                  r = Run()
                  r.num = int( run )
                  r.config = config
                  hlt_runs.append( r )

   return hlt_runs


# filter runs that are not in the dcs file
#
def getHLTDCSRunList( options, hlt_runs, dcs_runs ):
   hlt_dcs_runs = []

   if options.debug:
      print "hlt_runs:", hlt_runs
      print "dcs_runs:", dcs_runs

   for run in hlt_runs:
      if run.num in dcs_runs:
         hlt_dcs_runs.append( run )

   hlt_dcs_runs.sort()
   return hlt_dcs_runs


def getUnprescaledTriggers( run, options ):
   datasets = options.datasets[:]
   edmScript = 'edmConfigFromDB'
   config = subprocess.Popen( [ edmScript, '--format', 'python', '--orcoff', '--configName', run.config ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
   # The output of the edmConfigFromDB script delivered by 'communicate' is
   # python code. It is made accessible with help of 'exec'.
   #
   output = config.communicate()[0]
   #if options.debug:
      #print "From '%s'" %edmScript
      #print output
      #print "End from '%s'" %edmScript, '\n'
   exec output
   triggers = []
   for dataset in datasets[:]:
      try:
         triggers.append( [ dataset, getattr( process.datasets, dataset ) ] )
      except AttributeError:
         if options.verbose:
            print "WARNING: Caught '%s': %s" %( sys.exc_type, sys.exc_value )
            print "         This can happen when different HLT configs are used"
            print "         Igoring dataset: '%s' for config '%s'" %( dataset, run.config )
            print "         and proceeding..."
         datasets.remove( dataset )


   prescaled_triggers = []
   for pset in process.PrescaleService.prescaleTable:
      prescaled_triggers.append( pset.pathName.pythonValue().strip('\'' ) )  # ugly CMSSW stuff

   # initialise a dict with the datasets as keys and empty lists as values
   #
   #unprescaled_triggers = dict.fromkeys( datasets, list() )
   unprescaled_triggers = {}
   for dataset in datasets:
      unprescaled_triggers[ dataset ] = set()

   # fill the dict only with triggers that are not in the prescaled table
   # i.e. their prescale factor is 1
   # if the user has specified any filters for trigger paths they are also
   # removed from the trigger list
   #
   for dataset, trigList in triggers:
      for trig in trigList:
         if trig not in prescaled_triggers:
            addTrig = True
            if options.filter:
               addTrig = False
               for filter in options.filter:
                  if filter in trig:
                     addTrig = True
                     #unprescaled_triggers[ dataset ].add( trig )
            if options.neg_filters and addTrig:
               for neg_filter in options.neg_filters:
                  if neg_filter in trig:
                     addTrig = False
            if addTrig:
               unprescaled_triggers[ dataset ].add( trig )

   if options.debug:
      print 'prescaled_triggers:', '\n', prescaled_triggers, '\n'
      print 'unprescaled_triggers:', '\n', unprescaled_triggers, '\n'

   return unprescaled_triggers



# merge run ranges that have the same HLT config
#
def getRangesByConfig( hlt_dcs_runs ):

   def addNewRunRange( ranges, start_run, end_run ):
      range = Range()
      range.start_run = start_run.num
      range.end_run = end_run.num
      range.unprescaled_triggers = getUnprescaledTriggers( start_run, options )
      ranges.append( range )

   start_run = hlt_dcs_runs[0]
   end_run = hlt_dcs_runs[0]

   ranges = []
   for run in hlt_dcs_runs:
      if run.config == start_run.config:
         end_run = run
      else:
         addNewRunRange( ranges, start_run, end_run )

         start_run = run
         end_run = run
      if run == hlt_dcs_runs[-1]:
         addNewRunRange( ranges, start_run, end_run )

      print run

   return ranges


# merge run ranges that have the same unprescaled triggers
#
def getRangesByTrigger( run_ranges ):

   def addNewRange( new_ranges, start_range, end_range ):
      new_range = Range()
      new_range.start_run = start_range.start_run
      new_range.end_run = end_range.end_run
      new_range.unprescaled_triggers = start_range.unprescaled_triggers
      new_ranges.append( new_range )

   start_range = run_ranges[0]
   end_range = run_ranges[0]

   new_ranges = []
   for range in run_ranges:
      if range.unprescaled_triggers == start_range.unprescaled_triggers:
         end_range = range
      else:
         addNewRange( new_ranges, start_range, end_range )

         start_range = range
         end_range = range
      if range == run_ranges[-1]:
         addNewRange( new_ranges, start_range, end_range )


   return new_ranges


def printOutputFile( run_ranges_trig, options ):
   outFile = open( options.output, 'w' )
   print >> outFile, 'Script called with the parameters:'
   print >> outFile, options.fullCall
   print >> outFile
   print >> outFile, 'Number of ranges:', len( run_ranges_trig )
   print >> outFile

   all_unprescaled_triggers = set()

   for range in run_ranges_trig:
      print >> outFile
      rangesString = 'Range: %i-%i' % (range.start_run, range.end_run)
      print >> outFile, rangesString
      print >> outFile, '-' * len( rangesString )
      for dataset in range.unprescaled_triggers.keys():
         print >> outFile
         datasetString =  '%s:' % dataset
         print >> outFile, datasetString
         print >> outFile, '-' * len( datasetString )
         trigs = list( range.unprescaled_triggers[ dataset ] )
         trigs.sort()
         #print trigs
         for t in trigs:
            print >> outFile, t
            all_unprescaled_triggers.add( t )
      print

   triggerList = list( all_unprescaled_triggers )
   triggerList.sort()
   print >> outFile
   print >> outFile, 'All unprescaled triggers for all run ranges:'
   print >> outFile, '--------------------------------------------'
   for trig in triggerList:
      if trig == triggerList[ -1 ]:
         print >> outFile, "'%s'" %trig
      else:
         print >> outFile, "'%s'," %trig
   print >> outFile

   outFile.close()

   print "Output written to file: '%s'" %options.output



################################################################################
##################################### MAIN #####################################
################################################################################

dcs_runs = getDCSRunList( options.dcs_file )
hlt_runs = getHLTRunList( options )
hlt_dcs_runs = getHLTDCSRunList( options, hlt_runs, dcs_runs )

run_ranges = getRangesByConfig( hlt_dcs_runs )

run_ranges_trig = getRangesByTrigger( run_ranges )

printOutputFile( run_ranges_trig, options )
