#!/usr/bin/env python

import optparse
import sys

# for json support
try: # FUTURE: Python 2.6, prior to 2.6 requires simplejson
    import json
except:
    try:
        import simplejson as json
    except:
        print "Please use lxplus or set an environment (for example crab) with json lib available"
        sys.exit(1)


def createNewJSONFile( runmin, runmax, jsonfilename, prefix = 'DCS' ):
    new_file_name = prefix + '-' + str( runmin ) + '-' + str( runmax ) + '.json'

    json_dict = readJSONfile( jsonfilename )
    json_filtered = selectRunsFromJSON( json_dict, runmin, runmax )
    writeNewJSON( new_file_name, json_filtered )


def writeNewJSON( newfilename, jsonfiltered ):
    new_file = open( newfilename, 'w' )
    json.dump( jsonfiltered, new_file )
    new_file.close()

    print 'New filtered json file produced in:', newfilename


def selectRunsFromJSON( jsondict, runmin, runmax ):
    json_filtered = {}

    for run in jsondict.keys():
        print 'Reading run:', run
        if int( run ) >= runmin and int( run ) <= runmax:
            json_filtered[ run ] = jsondict[ run ]
            print "-----------> accepted"

    return json_filtered


def readJSONfile( jsonfilename ):
    json_dict = {}
    json_file = file( jsonfilename, 'r' )
    json_dict = json.load( json_file )

    return json_dict


################################################################################
################################# MAIN #########################################
################################################################################


def main():
    description = "Select a given run range from the given JSON file and create a new file."
    usage = "Useage: %prog <runmin> <runmax> <oldjson>"

    parser = optparse.OptionParser( description = description,  usage = usage, version = '%prog version 1' )
    parser.add_option( '-p', '--prefix', metavar = 'PREFIX', default = 'DCS',
                       help = 'Set the prefix for the output file name. The file name is always PREFIX-<runmin>-<runmax>. [default: %default]' )

    ( options, args ) = parser.parse_args()
    del parser

    run_min_str = args[0]
    run_max_str = args[1]
    JSON = args[2]

    if len( args ) < 3:
        parser.error( "ERROR: At least three arguments needed!" )

    if not run_min_str.isdigit():
        parser.error( "ERROR: <runmin> must be an integer!" )

    if not run_max_str.isdigit():
        parser.error( "ERROR: <runmax> must be an integer!" )

    run_min = int( run_min_str )
    run_max = int( run_max_str )

    if run_min > run_max_str:
        parser.error( "ERROR: <runmax> must be bigger than <runmin>!" )

    createNewJSONFile( run_min, run_max, JSON, options.prefix )

if __name__ == '__main__':
    main()
