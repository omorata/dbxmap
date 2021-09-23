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

import numpy as np
import copy

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
        
    G = set_pages(fw.Figure(cfg, dirs))

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



def set_pages(G):

    pyx = G.frame.yx
    outf = os.path.splitext(G.outfile)
    sz = pyx[0] * pyx[1]
    
    K = []
    
    p = G.frame.panels

    pages = int(np.ceil(len(p) / sz))
    
    if pages > 1:
        dig = int(np.log10(pages)+1)
        outf_fmt = '{0:s}-page_{1:0'+str(dig)+'d}{2:s}'
        
     
        for jp in range(pages) :
            K.append(copy.deepcopy(G))
            K[jp].frame.panels = []
            K[jp].outfile = outf_fmt.format(outf[0], jp, outf[1])

        for i in p:
            page = int((i.position[2] - 1) / sz)
            i.position = (pyx[0], pyx[1], i.position[2] - page * sz)
            K[page].frame.panels.append(i)

    else:
        K.append(G)
        
    return K


##-- End of functions --------------------------------------------------

