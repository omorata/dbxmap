#!/usr/bin/env python3
##
##  configure.py
##
##  O. Morata2020
##
##   Module to read and process dbxmap.py configuration files
##
##
import argparse
import os

import yaml

import framework as fw


##-- Functions ---------------------------------------------------------

def read_configuration_file(cfgfile):
    """ Read the YAML configuration file
    """
    print("  + reading configuration file:", cfgfile, "...")

    with open(cfgfile, 'r') as ymlfile:
        try:
            cnfg = yaml.safe_load(ymlfile)
        except yaml.YAMLError as exc:
            print(exc)
        
    return cnfg



def dump_config(dfile, cnfg) :
    """ dump the configuration in the logfile
    """ 

    with open(dfile, 'a+') as outfile:
        yaml.dump(cnfg, outfile, default_flow_style=False)


        
def process_config(ifiles, dirs) :
    """ process the configuration file config
    """

    cfg = read_configuration_file(ifiles['config'])

    if 'log' in ifiles:
        logfile = os.path.join(dirs['outdir'], ifiles['log'])
        print("  + dumping configuration file in lofgfile:", logfile)
        dump_config(logfile, cfg)
        
    G = fw.Figure(cfg, dirs)

    return G


def read_command_line() :
    """ read command line arguments
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', action='store',
                        help='configuration file', default='plot_cfg.yml')
    parser.add_argument('-l', '--log', help='log_file', default=None)
    parser.add_argument('-o', '--output_dir', help='output directory',
                        metavar='DIR', default='.')
    parser.add_argument('-w', '--work_dir', help='working directory',
                        default='.')
    return parser.parse_args()


##-- End of functions --------------------------------------------------

