#!/usr/bin/env python3
##
## framework.py
##
## O. Morata 2020
##
##  module containing the classes that build the basic structure of
##  the map
##

import os
import sys

import aplpy
import matplotlib.pyplot as plt
import numpy as np

import display as dsp
import markers as mrk


class Figure(object):
    """Define object Figure."""
    
    def __init__(self, cnfg, dirs):

        if 'page' in cnfg :
            fig_cfg= cnfg['page']
            
            if 'outfile' in fig_cfg:
                self.outfile = os.path.join(dirs['outdir'], fig_cfg['outfile'])
            else :
                self.outfile = os.path.join(dirs['outdir'], 'outfile.pdf')

            if 'dpi' in fig_cfg :
                self.dpi = fig_cfg['dpi']
            else:
                self.dpi = 100

            if 'size' in fig_cfg:
                self.size = fig_cfg['size']
            else :
                self.size = [10,10]
                    
            if 'name' in fig_cfg:
                self.name = fig_cfg['name']
            else :
                self.name = "NO_NAME"

        if 'panels' in cnfg:
            self.frame = Frame(cnfg['panels'], dirs)
        else:
            print("ERROR: no Frame definition")
            sys.exit(1)


            
    def create(self) :
        """Create figure."""

        print("  + creating figure [", self.name,"]")
        self.f = plt.figure(figsize=(self.size))
        self.frame.add_panels(self.f)
        return self.f



    def f_print(self):
        """Print figure. """
        
        print("  + plotting output file:", self.outfile)
        self.f.savefig(self.outfile, dpi=self.dpi)



    def end(self):
        """Ending comment."""
        
        print("\n  ... Figure [", self.name, "] done!\n")
    


        
class Frame(object):
    """ Class defining the elements of the figure
    """
    def __init__(self, cnfg, dirs):

        # self.labels

        self.wkdir = dirs['wkdir']
        
        if 'yx' in cnfg:
            self.yx = cnfg['yx']
        else:
            self.yx = [1,1]

        if 'view' in cnfg:
            self.view = View(cnfg['view'], self)

        if 'markers' in cnfg:
            self.markers = mrk.Marker(cnfg['markers'], self)

        if 'labels' in cnfg:
            self.labels = mrk.Label(cnfg['labels'], self)
        else :
            self.labels = None

        if 'pixrange' in cnfg:
            self.pixrange = dsp.Pixrange(cnfg['pixrange'], self)

        if 'contours' in cnfg:
            self.contour = dsp. Contour(cnfg['contours'], self)

            
        self.wd = dirs['wkdir']

        panel_str = [k for k in cnfg if 'panel' in k]

        if panel_str :
            panel_list = []
            p_idx = 0
            
            for panel in panel_str:
                print("    + adding panel:", panel, "...")
                panel_list.append(Panel(cnfg[panel], panel, p_idx, self))
                p_idx += 1

            self.panels = panel_list
        
        else:
            self.panels = []



    def add_panels(self, fig):

        p_idx = 0
        gc = []
        
        for p in self.panels :
            print("    + plotting panel:", p.name, "...")
            p.add_panel(fig, gc, p_idx)
            p_idx += 1


            
            
class Panel(object) :
    """ Create a panel
    """ 
    def __init__(self, cnfg, name, idx, parent):

        self.name = name

        
        if 'position' in cnfg:
            cpos = cnfg['position']
            self.position = (cpos[0], cpos[1], cpos[2])
        else:
            self.position = (parent.yx[0], parent.yx[1], idx+1)

            
        if 'view' in cnfg:
            self.view = View(cnfg['view'], parent)

        elif hasattr(parent, 'view'):
            self.view = View(None, parent)

        if 'axes' in cnfg:
            self.axes = Axes(cnfg['axes'], parent)
        elif hasattr(parent, 'axes'):
            self.axes = Axes(None, parent)
        else:
            self.axes = None
            
        if hasattr(parent, 'labels') and parent.labels != None:
            self.labels = parent.labels
            self.label_props = parent.labels.label_props.copy()
        else :
            self.labels = None

        if 'labels' in cnfg:
            self.labels = mrk.Label(cnfg['labels'], self)

            
        if 'markers' in cnfg:
            self.markers = mrk.Marker(cnfg['markers'], parent)
        elif hasattr(parent, 'markers'):
            self.markers = mrk.Marker(None, parent)

        if 'colorbar' in cnfg:
            self.colorbar = dsp.Colorbar(cnfg['colorbar'], self)
            
        dataset_list = []

        dataset_str = [k for k in cnfg if 'dataset' in k]

        if dataset_str :
            for d in dataset_str:
                dataset_list.append(dsp.Dataset(cnfg[d], d, parent))
                
            self.datasets = dataset_list
        else:
            self.datasets = []


            
    def add_panel(self, fig, gc, idx):
        """Adds a panel to the figure.

        Arguments:
            fig - figure to add the panel to
            gc - list of panels
            idx - index of panel
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

            try:
                # it has to be tested yet
                if hasattr(d, 'beam_args') :
                    gc[idx].add_beam(**d.beam_args)
                else :
                    gc[idx].add_beam()

            except KeyError:
                if hasattr(d, 'beam_args') :
                    gc[idx].add_beam(**d.beam_shape, **d.beam_args)
                    
                pass

        if self.axes != None:

            self.set_axes(gc, idx)



        if self.labels != None :
            for lb in self.labels.label_list :
                gc = lb.add_label(gc, idx)

        if self.markers != None:
            for mk in self.markers.marklist :
                mk.add_markers(gc, idx)

        if hasattr(self, 'colorbar'):
            self.set_colorbar(gc[idx])

            

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



    def set_colorbar(self, pnl) :
        pnl.add_colorbar()
        pnl.colorbar.set_width(self.colorbar.width)
        pnl.colorbar.set_location(self.colorbar.location)
        pnl.colorbar.set_axis_label_text(self.colorbar.text)


        
    def set_axes(self, gc, ix):
        """Sets axes properties."""

        ax = self.axes

        if hasattr(ax, 'xspacing'):
            gc[ix].ticks.set_xspacing(ax.xspacing)
        if hasattr(ax, 'yspacing'):
            gc[ix].ticks.set_yspacing(ax.yspacing)

        if hasattr(ax, 'xminor_freq'):
            if hasattr(ax, 'yminor_freq'):
                gc[ix].ticks.set_minor_frequency(ax.xminor_freq,ax.yminor_freq)
            else :
                gc[ix].ticks.set_minor_frequency(ax.xminor_freq)


        if hasattr(ax, 'tick_color'):
            gc[ix].ticks.set_color(ax.tick_color)

        if hasattr(ax, 'tick_length'):
            gc[ix].ticks.set_length(ax.tick_length)

        if hasattr(ax, 'tick_direction'):
            print(ax.tick_direction)
            gc[ix].ticks.set_tick_direction(ax.tick_direction)


        if hasattr(ax, 'xhide'):
            if ax.xhide :
                gc[ix].ticks.hide_x()

        if hasattr(ax, 'yhide'):
            if ax.yhide :
                gc[ix].ticks.hide_y()


                
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

        elif self.vtype == 'box' :
            gp[idx].recenter(self.center[0], self.center[1],
                             height=self.box[0], width=self.box[1])


            
class Axes(object):

    def __init__(self, cfg, parent):

        if cfg != None and 'minor_freq' in cfg :
            self.xminor_freq = cfg['minor_freq']

        if cfg != None and 'tick_color' in cfg :
            self.tick_color = cfg['tick_color']

        if cfg != None and 'tick_length' in cfg :
            self.tick_length = cfg['tick_length']

        if cfg != None and 'tick_direction' in cfg :
            self.tick_direction = cfg['tick_direction']

        if cfg != None and 'ticks_hide' in cfg:
            if cfg['ticks_hide']:
                self.xhide = True
                self.yhide = True
                

        if cfg != None and 'ticks_x' in cfg:
            tx = cfg['ticks_x']

            if 'spacing' in tx :
                self.xspacing = tx['spacing']

            if 'minor_freq' in tx:
                self.xminor_freq = tx['minor_freq']

            if 'hide' in tx:
                self.xhide = tx['hide']

        if cfg != None and 'ticks_y' in cfg:
            ty = cfg['ticks_y']

            if 'spacing' in ty :
                self.yspacing = ty['spacing']

            if 'minor_freq' in ty:
                self.yminor_freq = ty['minor_freq']

            if 'hide' in ty:
                self.yhide = ty['hide']
