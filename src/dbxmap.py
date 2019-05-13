#!/usr/bin/env python3
##
## dbxmap.py
## O. Morata 2019
##
## Tool to plot maps
##

import argparse

import matplotlib
import matplotlib.pyplot as plt
import aplpy
import numpy as np

import yaml

import os
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


print("a", os.getcwd())

    
## defaults
#
config_file = 'plot_cfg.yml'
configFilePath = (os.path.join(os.getcwd(), config_file))
print("b", configFilePath)
#matplotlib.use('Agg')

print(" Starting...")


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config', action='store',
                    help='configuration file', default='plot_cfg.yml')
parser.add_argument('-l', '--log', help='log_file', default='plot.log')
parser.add_argument('-o', '--output_dir', help='output directory',
                    metavar='DIR', default='.')
parser.add_argument('-w', '--work_dir', help='working directory', default='.')
args = parser.parse_args()


config_file = args.config
wkdir = args.work_dir
outdir = args.output_dir
logfile = args.log

print("  + reading configuration file...")

with open(config_file, 'r') as ymlfile:
    try:
        cfg = yaml.safe_load(ymlfile)
    except yaml.YAMLError as exc:
        print(exc)



for section in cfg:
    print(section)
    print(cfg[section])
    for sub in cfg[section]:
        print(sub)
        print(cfg[section][sub])


print(cfg['figure']['outfile'])
