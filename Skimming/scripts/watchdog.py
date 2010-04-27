#!/usr/bin/env python

import subprocess
import sys
import os
import shutil
import optparse
import ConfigParser


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
    job_list = [ x.num for x in jobs ]
    output = str( job_list )
    output = output[1:-1] #remove []
    output = output.replace( ' ','' ) #remove space after comma
    return output


def call_crab( command, jobs, dir, stdout=False ):
    jobs = format_job_list( jobs )
    if stdout:
        retcode = subprocess.call( [ 'crab', command, jobs, '-c', dir ] )
        if retcode != 0:
            print
            print 'Failed command:', 'crab', command, jobs, '-c', dir
            print 'Returncode:', retcode
    else:
        proc = subprocess.Popen( [ 'crab', command, jobs, '-c', dir ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
        output = proc.communicate()[0]
        if not proc.returncode == 0:
            print 'Failed command:', 'crab', command, jobs, '-c', dir
            print 'Returncode:', proc.returncode
            print '=====================Output:'
            print output
        else:
            return output.splitlines()


def make_state( job ):
    if job.state == 'Retrieved' or job.state == 'Cleared' or job.state == 'Done':
        if job.grid == 0 and job.exe == 0:
            return 'Success'
        elif job.exe != None and job.exe != 0:
            return 'App-Fail'
        elif job.grid != None and job.grid != 0 and ( not job.exe or job.exe == 0 ):
            return 'Grid-Fail'
        elif job.grid == None and job.exe == None and job.state == 'Done':
            return job.state
        else:
            raise JobStatusError( job )
    else:
        return job.state




def parse_output( output ):
    found_start = False
    jobs = defaultdict( list )

    for line in output:
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

        if len( split_line ) > 2:
            job.host = split_line[2]
        else:
            job.host = None

        if len( split_line ) == 4:
            job.grid = int( split_line[3] )
            job.exe = None
        elif len( split_line ) > 4:
            job.exe = int( split_line[3] )
            job.grid = int( split_line[4] )
        else:
            job.exe = None
            job.grid = None

        job.state = make_state( job )

        jobs[ job.state ].append( job )


    return jobs



def get_status( dir ):
    output = call_crab( '-status', None, dir )
    states = parse_output( output )

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
    call_crab( '-resubmit', jobs, dir, stdout=True )
                  
    stat.add_resubmitted( jobs )
                  
    print






parser = optparse.OptionParser( usage='usage: %prog [options] crab_dirs...' )
parser.add_option( '-a', '--resubmit-aborted', action='store_true', default=False, help='Resubmit aborted jobs' )
parser.add_option( '-g', '--resubmit-grid-failed', action='store_true', default=False, help='Resubmit jobs with grid failures' )
parser.add_option( '-f', '--resubmit-app-failed', action='store_true', default=False, help='Resubmit jobs with application failures' )
parser.add_option( '-u', '--user', help='Set the grid user name, in case it is not the same as the login name' )
(options, crab_dirs) = parser.parse_args()

if not options.user:
    options.user = os.getenv( 'LOGNAME' )
    if not options.user:
        print 'Cannot get user name, please provide one on the command line.'
        sys.exit(1)



stat = statistics()
moved_dirs = []

for dir in crab_dirs:
    print '=====', dir, '====='

    states = get_status( dir )

    if 'Done' in states:
        print 'Jobs done, getting output...'
        call_crab( '-get', states[ 'Done' ], dir, stdout=True )

        print
        print 'Got output, new states:'

        states = get_status( dir )


    stat.add_jobs( states )


    if 'Grid-Fail' in states:
        print_exit_codes( states[ 'Grid-Fail' ] )
    if 'App-Fail' in states:
        print_exit_codes( states[ 'App-Fail' ] )

    if len( states ) == 1 and 'Success' in states:
        moved_dirs.append( dir )
        move_dir( dir, 'done' )

    if options.resubmit_aborted and 'Aborted' in states:
        resubmit( dir, states[ 'Aborted' ], 'Aborted', options )

    if options.resubmit_grid_failed and 'Grid-Fail' in states:
        resubmit( dir, states[ 'Grid-Fail' ], 'Grid-Fail', options )

    if options.resubmit_app_failed and 'App-Fail' in states:
        resubmit( dir, states[ 'App-Fail' ], 'App-Fail', options )
        


stat.print_statistics()

if moved_dirs:
    print 'Dirs moved to done:', moved_dirs
