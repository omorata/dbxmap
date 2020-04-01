#!/usr/bin/env python3
##
##  dbxmap.py
##
##  O. Morata 2019-2020
##
##   Tool to plot maps using Aplpy
##
##

import matplotlib

import configure as cfg

# allow use of LaTex in the texts
#
matplotlib.rcParams['text.usetex'] = True


##-- Main --------------------------------------------------------------

if __name__ == "__main__" :
    ## defaults
    #
    #config_file = 'plot_cfg.yml'

    print(" Starting...")

    matplotlib.use('Agg')
    
    args = cfg.read_command_line()

    dirs = {'wkdir' : args.work_dir, 'outdir' : args.output_dir}
    files = {'config' : args.config}

    if args.log is not None:
        files['log'] = args.log

    # process the configuration
    #
    fig = cfg.process_config(files, dirs)

    # plot the figure
    #
    fig.create()

    fig.f_print()
    
    fig.end()
    
##-- End of main -------------------------------------------------------
