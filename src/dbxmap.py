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

from astropy import wcs
from astropy.io import fits

import yaml

import os
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


    
## defaults
#
config_file = 'plot_cfg.yml'
#configFilePath = (os.path.join(os.getcwd(), config_file))
#rint("b", configFilePath)
matplotlib.use('Agg')

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



#for section in cfg:
#    print(section)
#    print(cfg[section])
#    for sub in cfg[section]:
#        print(sub)
#        print(cfg[section][sub])


#print(cfg['figure']['outfile'])


hdulist = fits.open('test_map2.fits')
print(hdulist.info())
w = wcs.WCS(hdulist[0].header)

center = (w.wcs.crval[0], w.wcs.crval[1])

gc = []

# figure definition
#
fig = plt.figure(figsize=(8.25, 8.25))

gc.append(aplpy.FITSFigure('test_map2.fits', figure=fig, subplot=(1,1,1),
                       dimensions=[0,1]))


gc[0].recenter(center[0], center[1], radius=5./3600.)

gc[0].show_colorscale(vmin=-7e-6, vmax=7e-5, stretch='linear', cmap='inferno',
                      aspect='auto')

gc[0].show_contour(levels=[1e-5,2e-5, 3e-5, 4e-5, 5e-5,6e-5], colors="white")
gc[0].add_beam()
gc[0].save('out.pdf', dpi=100)

