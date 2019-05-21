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

class General(object):
    """ define object General
    """ 
    def __init__(self, cnfg):
        if 'outfile' in cnfg:
            self.outfile = cnfg['outfile']
        else :
            self.outfile = 'outfile.pdf'

        if 'dpi' in cnfg :
            self.dpi = cnfg['dpi']
        else:
            self.dpi = 100

        if 'xy' in cnfg:
            self.xy = cnfg['xy']
        else:
            self.xy = [1,1]

        if 'size' in cnfg:
            self.size = cnfg['size']
        else :
            self.size = [10,10]


    def create(self) :
        f = plt.figure(figsize=(self.size))
        return f

    
    def f_print(self, gc):
        print("  + plotting output file:", self.outfile, "...")
        gc[0].save(self.outfile, dpi=self.dpi)

    

        
class DbxFig(object):
    """ Class defining the elements of the figure
    """
    def __init__(self,cnfg, xy):

        # self.view
        # self.labels
        # self.contours
        # self.pixrange


        if 'view' in cnfg:
            self.view = View(cnfg['view'], self)

        #if 'labels' in cnfg:

        if 'pixrange' in cnfg:
            self.pixrange = Pixrange(cnfg['pixrange'],self)

        if 'contours' in cnfg:
            self.contour = Contour(cnfg['contours'], self)


            
        panel_list = []
    
        panel_str = [k for k in cnfg if 'panel' in k]

        if panel_str :

            p_idx = 0
            
            for panel in panel_str:

                print("  + adding panel:", panel, "...")
                panel_list.append(Panel(cnfg[panel], panel, xy, p_idx, self))
                p_idx += 1

            self.panels = panel_list
            

        else:
            self.panels = []

            
    def add_panels(self, fig, gc):

        p_idx = 0
        for p in self.panels :
            print("  + plotting panel:", p.name)
            p.add_panel(fig, gc, p_idx)
            p_idx += 1
            
        return 1


    
            
class Panel(object) :
    """ Create a panel
    """ 
    def __init__(self, cnfg, name, xy, idx, parent):

        # self.contours
        # self.pixrange
        
        self.name = name

        
        if 'position' in cnfg:
            cpos = cnfg['position']
            self.position = (cpos[0], cpos[1], cpos[2])
        else:
            self.position = (xy[0], xy[1], idx+1)

            
        if 'view' in cnfg:
            self.view = View(cnfg['view'], parent)

        elif hasattr(parent, 'view'):
            self.view = View(None, parent)
            
            
        if 'labels' in cnfg:

            label_list = []
            for l in cnfg['labels'] :
                label_list.append(Label(cnfg['labels'][l]))

            self.labels = label_list
        else:
            self.labels = []

            
        dataset_list = []

        dataset_str = [k for k in cnfg if 'dataset' in k]

        if dataset_str :
            for d in dataset_str:
                dataset_list.append(Dataset(cnfg[d], d, parent))
                
            self.datasets = dataset_list
        else:
            self.datasets = []


            
    def add_panel(self, fig, gc, idx):
        """ add a panel to the figure
        """
        self.set_pixel_first()

        vw = self.view

        cid = 0
        for d in self.datasets :
            if cid == 0 :
                if vw.center == None:
                    vw.center = d.get_reference()

                gc.append(aplpy.FITSFigure(d.filename,
                                   figure=fig,
                                   subplot=self.position,
                                   dimensions=d.dims))

                vw.set_view(gc, idx)

            d.show(gc, idx)

            cid += 1

        gc[idx].add_beam()

        for lb in self.labels :
            gc = lb.add_label(gc, idx)

            
            
    def set_pixel_first(self):
        """ put the dataset with a pixel range in the first position
        """
        idx = 0
        dp = self.datasets

        for d in dp:
            if d.dtype == 'pixel':
                if idx != 0 :
                    tmp = dp[0]
                    dp[0] = d
                    dp[idx] = tmp
                    self.datasets = dp
                    break

            idx += 1

    

            
class View(object) :
    """ Create a view
    """
    def __init__(self, cnfg, parent) :
        if cnfg != None and 'type' in cnfg:
            self.vtype = cnfg['type']
        else :
            try:
                self.vtype = parent.view.vtype
                
            except AttributeError:
                self.vtype = None

        if self.vtype == 'radius' :
            if cnfg != None and 'radius' in cnfg:
                self.radius = cnfg['radius']

            else:
                try:
                    self.radius = parent.view.radius
                except AttributeError:
                    self.radius = None

        elif self.vtype == 'box' :
            if cnfg != None  and 'box' in cnfg:
                self.box = cnfg['box']
            else :
                try:
                    self.box = parent.view.box
                except AttributeError:
                    self.box = None

        
            

        if cnfg != None and 'center' in cnfg:
            self.center = cnfg['center']
        else :
            try:
                self.center = parent.view.center
            except: 
                self.center = None

            

    def set_view(self, gp, idx):
        """set the view of the panel
        """
        if self.vtype == 'radius' :
            gp[idx].recenter(self.center[0], self.center[1],
                             radius=self.radius)

        elif sef.vtype == 'box' :
            gp[idx].recenter(self.center[0], self.center[1],
                             height=self.box[0], width=self.box[1])



            
class Dataset(object) :
    """ create a dataset object
    """
    def __init__(self, cnfg, name, parent):

        self.name = name

        
        if 'filename' in cnfg:
            self.filename = cnfg['filename']

        if 'dims' in cnfg:
            self.dims = cnfg['dims']
        else :
            self.dims = [0,1]

        if 'type' in cnfg:
            self.dtype = cnfg['type']

            if self.dtype == 'pixel' :
                if 'pixrange' in cnfg:
                    self.pixrange = Pixrange(cnfg['pixrange'], parent)
                elif hasattr(parent, 'pixrange'):
                    self.pixrange = Pixrange(None, parent)
                else :
                    raise Exception('Error: pixrange is not defined anywhere')
                    
            elif self.dtype == 'cntr' or self.dtype == 'contour':
                if 'contours' in cnfg:
                    self.contour = Contour(cnfg['contours'], parent)
                elif hasattr(parent, 'contour'):
                    self.contour = Contour(None, parent)
                else :
                    raise Exception('Error: contour is not defined anywhere')

                

    def get_reference(self):
        """ read the reference position of the datasetfrom the FITS
            header
        """
        hdulist = fits.open(self.filename)
        w = wcs.WCS(hdulist[0].header)

        center = (w.wcs.crval[0], w.wcs.crval[1])
        return center

    

    def show(self, gc, idx) :
        """ show the dataset
            It takes into account whether it is a pixel map or a contour
            map
        """
        if self.dtype == 'pixel' :
            self.show_colorscale(gc, idx)
        elif self.dtype == 'cntr' :
            self.show_contour(gc, idx)


            
    def show_colorscale(self, gc, idx) :
        """ show a colorscale of the dataset
        """
        gc[idx].show_colorscale(vmin=self.pixrange.range[0],
                                vmax=self.pixrange.range[1],
                                stretch=self.pixrange.stretch,
                                cmap=self.pixrange.colormap,
                                aspect='auto')


        
    def show_contour(self, gc, idx) :
        """ show a contour map of the dataset
        """
        gc[idx].show_contour(self.filename,
                             levels=self.contour.levels,
                             colors=self.contour.colors)
                             #linestyles=linestyle,
                             #linewidths=linewidth)
        

                    
class Pixrange(object):
    """ create a pixrange
    """
    def __init__(self, cnfg, parent):

        if cnfg != None and 'base' in cnfg:
            self.base = np.float(cnfg['base'])
        else :
            try:
                self.base = parent.pixrange.base

            except AttributeError:
                self.base = 1.

        if cnfg != None and 'colormap' in cnfg:
            self.colormap = cnfg['colormap']
            
        else:
            try:
                self.colormap = parent.pixrange.colormap

            except AttributeError:
                self.colormap = 'inferno'

        if cnfg != None and 'range' in cnfg:
            cfgrange = cnfg['range']

            self.range = [(float(i) * self.base) for i in cfgrange[0:2]]

            if np.shape(cfgrange)[0] == 3 :
                self.stretch = cfgrange[2]

            else :
                self.stretch = 'linear'

        else :
            try:
                self.range = parent.pixrange.range

            except AttributeError:
                pass

            try:
                self.stretch = parent.pixrange.stretch

            except AttributeError:
                pass


            
class Contour(object):
    """ create contours
    """
    def __init__(self,cnfg, parent) :

        if cnfg != None and 'base' in cnfg:
            self.base = float(cnfg['base'])
        else :
            try:
                self.base = parent.contour.base

            except AttributeError:
                self.base = 1.


        if cnfg != None and 'colors' in cnfg:
            self.colors = cnfg['colors']
        else :
            try:
                self.colors = parent.contour.colors

            except AttributeError :
                self.colors = 'black'

        if cnfg != None and 'levels' in cnfg:
            sc_levels = [(float(i) * self.base) for i in cnfg['levels']]
            
            if not hasattr(self, 'levels') :
                self.levels = sc_levels

            else:
                prev = self.levels
                prev.append(sc_levels)
                prev.sort()
                self.levels = prev
        else:
            try:
                self.levels = parent.contour.levels
                if self.base != 1 :
                    self.levels = [(i * self.base) for i in self.levels]
                
            except:
                pass

        if cnfg != None and 'gen_levels' in cnfg:
            lev_range = [(float(i) * base) for i in cnfg['gen_levels']]

            lv = lev_range[0]
            inc = lev_range[2]
            
            if not hasattr(self, 'levels'):
                levs = []
                while np.sign(inc) * (lev_range[1] - lv) > -1e-10 :
                    levs.append(lv)
                    lv += inc

                self.levels = levs

            else :
                prev = self.levels
                while np.sign(inc) * (lev_range[1] - lv) > -1e-10 :
                    prev.append(lv)
                    lv += inc

                prev.sort()
                self.levels = prev




class Label(object):
    """ create a label
    """
    def __init__(self, cfg):
        if 'text' in cfg :
            self.text = cfg['text']
        else :
            self.text =""

        if 'relative' in cfg:
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

        if 'props' in cfg:
            self.props = cfg['props']
                             


    def add_label(self, pp, idx) :
        pp[idx].add_label(self.position[0], self.position[1],
                          self.text,
                          relative=self.relative,
                          color=self.color,
                          size=self.size,
                          style=self.style)

        return pp

        
##-- Functions ---------------------------------------------------------
##

def read_configuration_file(file):
    """ Read the YAML configuration_file file
    """
    print("  + reading configuration file:", file, "...")

    with open(file, 'r') as ymlfile:
        try:
            cnfg = yaml.safe_load(ymlfile)
        except yaml.YAMLError as exc:
            print(exc)
        
    # modify figure keyword
    #
    if 'figure' in cnfg :
        cnfg['global'] = cnfg['figure']
        del cnfg['figure']
    else :
        cnfg['global'] = {}

    return cnfg



def process_config(config) :
    """ process the configuration file config
    """

    cfg = read_configuration_file(config)

    if 'global' in cfg:
        G = General(cfg['global'])
    else :
        G = General()

    if 'panels' in cfg:
        P = DbxFig(cfg['panels'], G.xy)

    return G, P



##-- End of functions --------------------------------------------------


    
## defaults
#
config_file = 'plot_cfg.yml'

print(" Starting...")

matplotlib.use('Agg')



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


General, Panels = process_config(config_file)


gc = []


# figure definition
#
fig = General.create()

a = Panels.add_panels(fig, gc)

General.f_print(gc)


