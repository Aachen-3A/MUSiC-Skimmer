#!/usr/bin/env python

import subprocess
import sys
import os
import shutil
import optparse
import ConfigParser
import xml.dom.minidom
import time
import MUSiCProject.Skimming.outputCleaner as outputCleaner
import logging


class JobStatusError( Exception ):
    def __init__( self, job ):
        self.job = job
    def __str__( self ):
        return 'Job with a weird status. Please report this with the following information:\n'+str(self.job)

class defaultdict( dict ):
    def __init__( self, factory, *args ):
        self.factory = factory
        dict.__init__( self, args )

    def __getitem__( self, key ):
        try:
            return dict.__getitem__( self, key )
        except KeyError:
            self[ key ] = self.factory()
            return self[ key ]


class Job:
    num = None
    state = None
    host = None
    grid = None
    exe = None

    def __str__( self ):
        return ( 'Job: '+
                 'num='+str( self.num )+'; '+
                 'state='+str( self.state )+'; '+
                 'host='+str( self.host )+'; '+
                 'grid_exit='+str( self.grid )+'; '+
                 'exe_exit='+str( self.exe )
                 )


def format_job_list( jobs ):
    if jobs == None: return ''
    if jobs == 'all': return jobs
    job_list = [ x.num for x in jobs ]
    output = str( job_list )
    output = output[1:-1] #remove []
    output = output.replace( ' ','' ) #remove space after comma
    return output


def call_crab( command, jobs, dir, stdout=False, additionals=None ):
    jobs = format_job_list( jobs )
    cmndlst = [ 'crab', command, jobs, '-c', dir ]
    if additionals:
        cmndlst += additionals
    if stdout:
        retcode = subprocess.call( cmndlst )
        if retcode != 0:
            print
            print 'Failed command:', cmndlst
            print 'Returncode:', retcode
    else:
        proc = subprocess.Popen( cmndlst, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
        output = proc.communicate()[0]
        if not proc.returncode == 0:
            print 'Failed command:', cmndlst
            print 'Returncode:', proc.returncode
            print '=====================Output:'
            print output
        else:
            return output.splitlines()


def make_state( job, use_server ):
    if job.state == 'Retrieved' or job.state == 'Cleared' or (not use_server and job.state == 'Done'):
        if job.grid == 0 and job.exe == 0:
            return 'Success'
        elif job.exe != None and job.exe != 0:
            return 'App-Fail'
        elif job.grid != None and job.grid != 0 and ( job.exe == None or job.exe == 0 ):
            return 'Grid-Fail'
        elif job.grid == None and job.exe == None and job.state == 'Done':
            return job.state
        else:
            raise JobStatusError( job )
    else:
        return job.state




def parse_output( output, use_server ):
    found_start = False
    jobs = defaultdict( list )

    for line in output:
        if 'WARNING: Error getting the status from server. Please issue crab -status again' in line:
            #CRAB failed to get the status, so what we got is outdated
            return None

        if not found_start:
            if 'STATUS' in line:
                found_start = True
            else:
                continue
    
        if '------' in line:
            continue

        if 'crab' in line:
            break

        split_line = line.split()
        if len( split_line ) == 0:
            continue

        job = Job()

        job.num = split_line[0]
        if job.num.isdigit():
            job.num = int( job.num )
        else:
            continue

        job.state = split_line[1]
        if job.state == 'Cancelled':
            del split_line[2:4]

        if len( split_line ) == 2:
            job.host = None
            job.exe = None
            job.grid = None
        elif len( split_line ) == 3:
            if split_line[2].isdigit():
                job.host = None
                job.exe = None
                job.grid = int( split_line[2] )
            else:
                job.host = split_line[2]
                job.exe = None
                job.grid = None
        elif len( split_line ) == 4:
            if split_line[2].isdigit():
                job.host = None
                job.exe = int( split_line[2] )
                job.grid = int( split_line[3] )
            else:
                job.host = split_line[2]
                job.exe = None
                job.grid = int( split_line[3] )
        elif len( split_line ) == 5:
            job.host = split_line[2]
            job.exe = int( split_line[3] )
            job.grid = int( split_line[4] )
        else:
            raise Exception( 'Unexpected line in crab output: %s' % line )

        job.state = make_state( job, use_server )

        jobs[ job.state ].append( job )


    return jobs



def get_status( dir, use_server ):
    for trial in xrange( 3 ):
        print 'Retrieving status (trial %s)...' % (trial+1)
        output = call_crab( '-status', None, dir )
        states = parse_output( output, use_server )
        if states: break
    else:
        print 'Failed 3 times to get the status!'
        return None

    for state,jobs in states.iteritems():
        print state, format_job_list( jobs )

    print

    return states



def print_exit_codes( jobs ):
    grid_codes = defaultdict( list )
    exe_codes = defaultdict( list )
    for job in jobs:
        if job.exe and job.exe != 0:
            exe_codes[ job.exe ].append( job.num )
        else:
            grid_codes[ job.grid ].append( job.num )

    if exe_codes:
        print 'Application failures:'
        for code,jobs in exe_codes.iteritems():
            print code, ':', jobs
        print
    if grid_codes:
        print 'Grid failures:'
        for code,jobs in grid_codes.iteritems():
            print code, ':', jobs
        print



class statistics:
    def __init__( self ):
        self.states = defaultdict( int )
        self.resub_states = defaultdict( int )

        def def_dict_fact():
            return defaultdict( int )
        self.sites = defaultdict( def_dict_fact )
        self.resub_sites = defaultdict( def_dict_fact )
    def add_jobs( self, states ):
        for jobs in states.values():
            for job in jobs:
                self.states[ job.state ] += 1
                self.sites[ job.host ][ job.state ] += 1
    def add_resubmitted( self, jobs):
        for job in jobs:
            self.resub_states[ job.state ] += 1
            self.resub_sites[ job.host ][ job.state ] += 1

    def print_statistics( self ):
        print '\n\n\n++++++ Global statistics ++++++'
        
        print '\n++++ States ++++'
        for state, num in self.states.items():
            print state, num,
            if state in self.resub_states:
                print '(Resubmitted', self.resub_states[ state ], ')'
            else:
                print

        print '\n++++ States by sites ++++'
        for site, states in self.sites.items():
            for state, num in states.items():
                print site, state, num,
                if site in self.resub_sites and state in self.resub_sites[ site ]:
                    print '(Resubmitted', self.resub_sites[ site ][ state ], ')'
                else:
                    print

        print '\n\n\n'


def move_dir( dir, target ):
    if not os.access( target, os.F_OK ):
        os.mkdir( target )
    if not os.access( target, os.F_OK | os.W_OK ):
        print "Can't move, something is wrong with the target directory:", target
        return
    shutil.move( dir, os.path.join( target, dir ) )



       



def resubmit( dir, jobs, state, options ):
    print 'Resubmitting jobs in state:', state
    if options.blacklist:
        bl_cmd = [ '-GRID.ce_black_list', options.blacklist ]
    else:
        bl_cmd = None
    call_crab( '-resubmit', jobs, dir, stdout=True, additionals=bl_cmd )
                  
    stat.add_resubmitted( jobs )
                  
    print





def clean_output( dir, parser, jobs ):
    if parser.getboolean( 'USER', 'copy_data' ) and not options.noclean:
        print 'Cleaning output'
        remote_dir = parser.get( 'USER', 'user_remote_dir' )
        outputCleaner.clean( dir, [ job.num for job in jobs ], remote_dir, clear_all=True )
        print 'Done\n\n'

    moved_dirs.append( dir )
    move_dir( dir, 'done' )





parser = optparse.OptionParser( usage='usage: %prog [options] crab_dirs...' )
parser.add_option( '-a', '--resubmit-aborted', action='store_true', default=False, help='Resubmit aborted jobs' )
parser.add_option( '-g', '--resubmit-grid-failed', action='store_true', default=False, help='Resubmit jobs with grid failures' )
parser.add_option( '-f', '--resubmit-app-failed', action='store_true', default=False, help='Resubmit jobs with application failures' )
parser.add_option( '-r', '--resubmit-running', action='store_true', default=False, help='Kill running jobs and resubmit them' )
parser.add_option( '-s', '--resubmit-scheduled', action='store_true', default=False, help='Kill scheduled jobs and resubmit them' )
parser.add_option( '--kill-all', action='store_true', default=False, help='Issue -kill all for all tasks' )
parser.add_option( '-u', '--user', help='Set the grid user name, in case it is not the same as the login name' )
parser.add_option( '-b', '--blacklist', metavar='STRING', help='Blacklist STRING during resubmission.' )
parser.add_option( '-n', '--noclean', action='store_true', default=False, help='Do not check for and delete files left over by resubmission' )
parser.add_option( '--debug', metavar='LEVEL', default='INFO', help='Set the debug level. Allowed values: ERROR, WARNING, INFO, DEBUG. [default: %default]' )
(options, crab_dirs) = parser.parse_args()

#set log level
logging.basicConfig( level=logging._levelNames[ options.debug ] )


if not options.user:
    options.user = os.getenv( 'LOGNAME' )
    if not options.user:
        print 'Cannot get user name, please provide one on the command line.'
        sys.exit(1)

if options.kill_all and ( options.resubmit_aborted or options.resubmit_grid_failed or options.resubmit_app_failed or options.kill_resubmit ):
    print '--kill-all must be the only command'
    sys.exit(1)


stat = statistics()
moved_dirs = []

for dir in crab_dirs:
    print '=====', dir, '====='

    #if we just need to kill, then do so
    if options.kill_all:
        call_crab( '-kill', 'all', dir, stdout=True )
        #and continue with the next task
        continue

    parser = ConfigParser.SafeConfigParser()
    parser.read( os.path.join( dir, 'share/crab.cfg' ) )
    try:
        use_server = parser.getboolean( 'CRAB', 'use_server' )
    except ConfigParser.NoOptionError:
        use_server = False

    states = get_status( dir, use_server )
    if states == None:
        print 'Skipping', dir
        print
        continue

    if 'Done' in states:
        print 'Jobs done, getting output...'
        call_crab( '-get', states[ 'Done' ], dir, stdout=True )

        print
        if use_server:
            print 'Got output, sleeping 10 seconds to allow server to settle.'
            time.sleep( 10 )
            print 'Getting new states:'
        else:
            print 'Got output, new states:'

        states = get_status( dir, use_server )
        if states == None:
            print 'Skipping', dir
            continue


    stat.add_jobs( states )


    if 'Grid-Fail' in states:
        print_exit_codes( states[ 'Grid-Fail' ] )
    if 'App-Fail' in states:
        print_exit_codes( states[ 'App-Fail' ] )

    if len( states ) == 1 and 'Success' in states:
        clean_output( dir, parser, states[ 'Success' ] )

    if options.resubmit_aborted and 'Aborted' in states:
        resubmit( dir, states[ 'Aborted' ], 'Aborted', options )

    if options.resubmit_grid_failed and 'Grid-Fail' in states:
        resubmit( dir, states[ 'Grid-Fail' ], 'Grid-Fail', options )

    if options.resubmit_app_failed and 'App-Fail' in states:
        resubmit( dir, states[ 'App-Fail' ], 'App-Fail', options )

    to_kill = []
    if options.resubmit_running:
        to_kill += states[ 'Running' ]
    if options.resubmit_scheduled:
        to_kill += states[ 'Scheduled' ]
    if to_kill:
        call_crab( '-kill', to_kill, dir, stdout=True )
        resubmit( dir, to_kill, 'Killed', options )

stat.print_statistics()

if moved_dirs:
    print 'Dirs moved to done:', moved_dirs
