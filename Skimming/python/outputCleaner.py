#!/usr/bin/env python

import os.path
import xml.dom.minidom
import subprocess
import logging
import ConfigParser
import optparse
import glob
import collections
import sys


log = logging.getLogger( "outputCleaner" )

def clean( crab_dir, jobs, r_dir, clear_all=False, simulate=False ):
    if r_dir.endswith( '/' ):
        remote_dir = r_dir[ :-1 ]
    else:
        remote_dir = r_dir

    log.info( 'Parsing framework job reports' )
    valid_files = collections.defaultdict( set )
    for job in jobs:
        fjr_file_name = 'crab_fjr_'+str( job )+'.xml'
        log.debug( 'Parsing: %s', fjr_file_name )
        fjr_file = os.path.join( crab_dir, 'res', fjr_file_name )
        dom = xml.dom.minidom.parse( fjr_file )
        for analysisFile in dom.getElementsByTagName( 'AnalysisFile' ):
            PFN = analysisFile.getElementsByTagName( 'PFN' )[0].getAttribute( 'Value' )
            log.debug( 'Found valid file: '+PFN )
            dir_name, file_name = os.path.split( PFN )
            if dir_name.endswith( remote_dir ):
                valid_files[ dir_name ].add( file_name )
            else:
                log.warning( 'PFN %s in %s does not match configured user_remote_dir %s and will be ignored.', PFN, fjr_file_name, remote_dir )
    
    log.info( 'Building list of files on storage element' )
    existing_files = {}
    for dir_name in valid_files.keys():
        log.debug( 'srmls '+dir_name )
        proc = subprocess.Popen( [ 'srmls', dir_name ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
        output = proc.communicate()[0]
        files = set()
        if proc.returncode != 0:
            if 'does not exist' in output:
                #the directory doesn't exist, so there can't be any relic files
                log.warning( 'Directory not found, probably no output has been written.' )
                log.warning( 'Missing directory: %s', dir_name )
            else:
                log.error( 'Something went wrong with getting the stageout directory content: %s', dir_name )
                log.error( 'Output of srmls:' )
                log.error( output )
                sys.exit(1)
        else:
            for line in output.splitlines():
                if 'WARNING: SRM_PATH' in line: continue
                if '  0' in line: continue
                if not line: continue
                files.add( os.path.basename( line.split()[1] ) )
                log.debug( 'Found file: '+os.path.basename( line.split()[1] ) )
        existing_files[ dir_name ] = files
    
    log.info( 'Comparing directories' )
    relic_files = {}
    relic_files_present = False
    missing_files = {}
    missing_files_present = False
    unknown_files = {}
    unknown_files_present = False
    for dir_name in valid_files.keys():
        log.info( 'Comparing for directory %s', dir_name )

        v_f = valid_files[ dir_name ]
        e_f = existing_files[ dir_name ]

        m_f = v_f - e_f
        missing_files[ dir_name ] = m_f
        if len( m_f ) > 0:
            missing_files_present = True
        log.info( '%i files missing', len( m_f ) )

        if clear_all:
            r_f = e_f - v_f
        else:
            r_f = set()
            for file_name in v_f:
                full_name = os.path.join( dir_name, file_name )
                log.debug( 'Checking file: %s', full_name )
                if file_name in m_f:
                    log.warning( 'File %s is in a framework job report, but missing on the storage element. Relic files of this job will be ignored.', full_name )
                else:
                    log.debug( 'Valid file exists: %s', full_name )
                    root, ext = os.path.splitext( file_name )
                    root_part, num = root.rsplit( '_', 1 )
                    num = int( num )
                    for i in xrange( num ):
                        check_name = root_part + '_' + str( i ) + ext
                        log.debug( 'Checking potential relic file %s', os.path.join( dir_name, check_name ) )
                        if check_name in e_f:
                            log.info( 'Relic file: '+os.path.join( dir_name, check_name ) )
                            r_f.add( check_name )
        relic_files[ dir_name ] = r_f
        if len( r_f ) > 0:
            relic_files_present = True
        log.info( '%i relic files', len( r_f ) )

        u_f = e_f - v_f - r_f
        unknown_files[ dir_name ] = u_f
        if len( u_f ) > 0:
            unknown_files_present = True
        log.info( '%i unknown files', len( r_f ) )


    if relic_files_present:
        log.info( 'Deleting relic files' )
        for dir_name, files in relic_files.items():
            for file_name in files:
                full_name = os.path.join( dir_name, file_name )
                if simulate:
                    print 'srmrm', full_name
                else:
                    log.info( 'srmrm %s', full_name )
                    proc = subprocess.Popen( [ 'srmrm', full_name ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
                    output = proc.communicate()[0]
                    if proc.returncode != 0:
                        log.error( 'Failed to delete %s', full_name )
                        log.error( 'srmrm output:' )
                        log.error( output )
                        sys.exit(1)

    if missing_files_present:
        log.warning( 'There are files in the framework job reports that are not in the output directory:' )
        for dir_name, files in missing_files.items():
            for file_name in files:
                log.warning( os.path.join( dir_name, file_name ) )
        log.warning( 'Relic versions of those missing files have not been deleted!' )

    if unknown_files_present:
        log.warning( 'The following files have not been deleted but are not in any analyzed framework job report:' )
        for dir_name, files in unknown_files.items():
            for file_name in files:
                log.warning( os.path.join( dir_name, file_name ) )
        log.warning( 'Those files have not been deleted!' )




def build_job_list_from_crab( crab_dir, use_server ):
    log.info( 'Calling crab -status -c %s', crab_dir )
    proc = subprocess.Popen( [ 'crab', '-status', '-c', crab_dir ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
    output = proc.communicate()[0]
    if not proc.returncode == 0:
        log.error( 'Failed command: crab -status -c %s', crab_dir )
        log.error( 'Returncode: %s', proc.returncode )
        log.error( '=====================Output:' )
        log.error( output )
        sys.exit(1)
    
    log.info( 'Parsing output' )
    found_start = False
    jobs = []
    for line in output.splitlines():
        log.debug( line )
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

        num = split_line[0]
        if num.isdigit():
            num = int( num )
        else:
            continue
        log.debug( 'Found job with num %i', num )

        if len( split_line ) > 4:
            state = split_line[1]
            exe = int( split_line[3] )
            grid = int( split_line[4] )
        else:
            continue

        if state == 'Retrieved' or state == 'Cleared' or (not use_server and state == 'Done'):
            if grid == 0 and exe == 0:
                log.debug( 'Adding job %i to list', num )
                jobs.append( num )

    return jobs


def build_job_list_from_glob( crab_dir ):
    log.info( 'Looking for framework job reports in %s', crab_dir )
    jobs = []
    for path in glob.iglob( crab_dir+'/res/crab_fjr_*.xml' ):
        log.debug( 'Found framework job report file %s', path )
        #get the file name without path and extention
        file_name = os.path.splitext( os.path.basename( path ) )[0]
        num = int( file_name.rsplit( '_', 1 )[1] )
        log.debug( 'Adding job %i to list', num )
        jobs.append( num )
    return jobs



def build_job_list_from_list( jobs ):
    job_list = []
    for par in jobs.split( ',' ):
        if '-' in par:
            start, end = par.split( '-' )
            log.debug( 'Adding job range %s-%s to list', start, end )
            job_list += range( int(start), int(end)+1 )
        else:
            log.debug( 'Adding job %s to list', par )
            job_list.append( int( par ) )
    return job_list



def clean_dirs( dirs, options ):
    for crab_dir in dirs:
        log.info( 'Cleaning output of %s', crab_dir )
        parser = ConfigParser.SafeConfigParser()
        parser.read( os.path.join( crab_dir, 'share/crab.cfg' ) )
        if parser.getboolean( 'USER', 'copy_data' ):
            remote_dir = parser.get( 'USER', 'user_remote_dir' )
            if options.jobs == None:
                log.info( 'Cleaning all successful jobs' )
                jobs = build_job_list_from_crab( crab_dir, parser.getboolean( 'CRAB', 'use_server' ) )
                clean( crab_dir, jobs, remote_dir, clear_all=options.clear, simulate=options.dry_run )
            elif options.jobs == 'all':
                log.info( 'Cleaning all jobs with returned framework job reports' )
                jobs = build_job_list_from_glob( crab_dir )
                clean( crab_dir, jobs, remote_dir, clear_all=options.clear, simulate=options.dry_run )
            else:
                log.info( 'Cleaning jobs specified in %s', options.jobs )
                jobs = build_job_list_from_list( options.jobs )
                clean( crab_dir, jobs, remote_dir, clear_all=options.clear, simulate=options.dry_run )
        else:
            log.warning( 'The task in'+crab_dir+'has not been configured to return data, so no cleaning possible.' )



def check_proxy():
    log.info( 'Checking proxy' )
    proc = subprocess.Popen( ['voms-proxy-info' ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    output = proc.communicate()[0]
    log.debug( output )
    return proc.returncode == 0


def main():
    parser = optparse.OptionParser( usage='usage: %prog [options] crab_dirs...',
                                    description='Delete left-over output files of the jobs in crab_dirs.' )
    parser.add_option( '-j', '--jobs', metavar='JOBLIST', help='CRAB like list of jobs' )
    parser.add_option( '-c', '--clear', action='store_true', default=False, help='Delete all files that are not in a framework job report, but are in a directory that contains at least one file that is.' )
    parser.add_option( '-n', '--dry-run', action='store_true', default=False, help='Do not actually delete anything, just report what would be deleted. [default: %default]' )
    parser.add_option( '--debug', metavar='LEVEL', default='WARNING', help='Set the debug level. Allowed values: ERROR, WARNING, INFO, DEBUG. [default: %default]' )
    options, crab_dirs = parser.parse_args()

    #set log level
    logging.basicConfig( level=logging._levelNames[ options.debug ], stream=sys.stdout )

    if options.dry_run:
        logging.warning( 'Dry run: Nothing will be deleted!' )

    if check_proxy():
        log.info( 'Proxy valid' )
    else:
        log.error( 'Proxy not valid. Make a new one!' )
        sys.exit(1)

    clean_dirs( crab_dirs, options )


if __name__ == '__main__':
    main()
