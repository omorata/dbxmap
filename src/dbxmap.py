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

#import os
#try:
#    from yaml import CLoader as Loader, CDumper as Dumper
#except ImportError:
#    from yaml import Loader, Dumper


class Label(object):

    def __init__(self, cfg):
        if 'text' in cfg :
            self.text = cfg['text']
        else :
            self.text =""

        if 'text' in cfg:
            self.relative = cfg['relative']
        else :
            self.relative = True

        if 'position' in cfg:
            self.position = cfg['position']
        else :
            print(" ERROR")

        if 'color' in cfg :    
            self.color = cfg['color']
        else :
            self.color = 'black'

        if 'size' in cfg:
            self.size = cfg['size']
        else :
            self.size = 12

        if 'style' in cfg:
            self.style = cfg['style']
        else :
            self.style = 'normal'


    def add_label(self, pp, idx) :
        pp[idx].add_label(self.position[0], self.position[1],
                          self.text,
                          relative=self.relative,
                          color=self.color,
                          size=self.size,
                          style=self.style)

        return pp

        
##-- Functions ---------------------------------------------------------

def add_panel(xx, pstn, xcfg, idx) :
    """ Add a new panel
    """

    ini_panel = False
    pxflg = False
    ctrflg = False
    
    dataset_str = [k for k in xcfg if 'dataset' in k]
    for dset in dataset_str:
        dset_cfg =  xcfg[dset]

        ctrset = []
        
        if dset_cfg['type'] == 'pixel' :
            pixset = dset

            if not pxflg :
                pxflg = True
            else :
                print(" ERROR: two pixel datasets cannot be defined"
                      "simultaneously")
                
        elif dset_cfg['type'] == 'cntr':
            ctrset.append(dset)
            ctrflg = True

        else :
            print(" **ERROR: plotting option unknown for ", dset)

            
    if pxflg :
        xx, ini_panel = add_dataset(xx, pstn, xcfg, idx, pixset, ini_panel)

    if ctrflg:
        for ctr in ctrset :
            xx, ini_panel = add_dataset(xx, pstn, xcfg, idx, ctr, ini_panel)


    if 'labels' in xcfg:
        it_lab = [k for k in xcfg['labels'] if 'lab' in k]

        for lab in it_lab:
            l = Label(xcfg['labels'][lab])
            xx = l.add_label(xx, idx)

            
    return xx



def add_dataset(xx, pos, dcfg, idx, data, ini_flag) :
    """ Add a new dataset to a panel
    """
    fname = dcfg[data]['filename']
    dims = dcfg[data]['dims']
    pltype = dcfg[data]['type']

    if not ini_flag :
        hdulist = fits.open(fname)
        #print(hdulist.info())
        w = wcs.WCS(hdulist[0].header)

        center = (w.wcs.crval[0], w.wcs.crval[1])

    
        xx.append(aplpy.FITSFigure(fname,
                                   figure=fig,
                                   subplot=pos,
                                   dimensions=dims))

        if 'view' in dcfg :
            xx[idx].recenter(center[0], center[1],
                             dcfg['view']['radius'])

        ini_flag = True


    if pltype == 'pixel':
        colmap = 'inferno'
        xx[idx].show_colorscale(vmin=-7e-6, vmax=7e-5,
                                stretch='linear',
                                cmap=colmap,
                                aspect='auto')

    elif pltype == 'cntr' :
        linecolor = 'white'
        linestyle = None
        linewidth = 1
        if 'cntr_props' in dcfg[data] :
            if 'color' in dcfg[data]['cntr_props']:
                linecolor = dcfg[data]['cntr_props']['color']
            if 'style' in dcfg[data]['cntr_props']:
                linestyle = dcfg[data]['cntr_props']['style']
            if 'width' in dcfg[data]['cntr_props']:
                linewidth = dcfg[data]['cntr_props']['width']

                
        xx[idx].show_contour(fname, levels=[-2e-5, 1e-5,2e-5, 3e-5, 4e-5, 5e-5,
                                     6e-5],
                             colors=linecolor,
                             linestyles=linestyle,
                             linewidths=linewidth)
            
        # TODO
        # check if in every panel
        # do not draw it if it is the same
        # how to plot several beams
        #
    xx[idx].add_beam()

    return (xx, ini_flag)


##-- End of functions --------------------------------------------------


    
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

print("  + reading configuration file:", config_file, "...")

with open(config_file, 'r') as ymlfile:
    try:
        cfg = yaml.safe_load(ymlfile)
    except yaml.YAMLError as exc:
        print(exc)


        
#for section in cfg:
#    if section == 'figure':
#        print(cfg[section].keys())
#    elif section == 'panels':
#        print(cfg[section].keys())
#    print(section)
#    print(cfg[section])
#    for sub in cfg[section]:
#        print(sub)
#        print(cfg[section][sub])

cfg['global'] = cfg['figure']
del cfg['figure']

#process cfg


gc = []

panel_idx = 0


# figure definition
#
fig = plt.figure(figsize=(cfg['global']['size']))


panel_str = [k for k in cfg['panels'] if 'panel' in k]

#defined_panels = np.shape(panel_str)[0]

# what to do if defined_panels is different than numpanels

for panel in panel_str:

    print("  + adding panel:", panel, "...")
    
    position = (cfg['global']['xy'][0], cfg['global']['xy'][1], panel_idx+1)

    panel_cfg = cfg['panels'][panel]

    gc = add_panel(gc, position, panel_cfg, panel_idx)

    panel_idx += 1


ofile = cfg['global']['outfile']
print("  + plotting output file:", ofile, "...")
gc[0].save(ofile, dpi=cfg['global']['dpi'])

