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
            ax = self.axes
            ax.set_axes(gc, idx)

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

        #if cfg != None and 'axes_labels' in cfg
        #    self.read_axes_labels(cfg['axes_labels'])

        #if cfg != None and 'tick_labels' in cfg:
        #    self.read_tick_labels(cfg['tick_labels'])

        if cfg != None and 'ticks' in cfg :
            self.read_ticks(cfg['ticks'])



    def read_ticks (self, tx):

            if 'xspacing' in tx :
                self.xspacing = tx['xspacing']

            if 'xminor_freq' in tx:
                self.xminor_freq = tx['xminor_freq']

            if 'yspacing' in tx :
                self.yspacing = tx['yspacing']

            if 'yminor_freq' in tx:
                self.yminor_freq = tx['yminor_freq']

            if 'color' in tx:
                self.tick_color = tx['color']

            if 'length' in tx:
                self.tick_length = tx['length']

            if 'linewidth' in tx:
                self.tick_linewidth = tx['linewidth']

            if 'direction' in tx:
                self.tick_direction = tx['direction']

            if 'hide' in tx:
                self.tick_xhide = tx['hide']
                self.tick_yhide = tx['hide']

            if 'xhide' in tx:
                self.tick_xhide = tx['xhide']

            if 'yhide' in tx:
                self.tick_yhide = tx['yhide']



    def set_axes(self, gc, ix):
        """Sets axes properties."""

        #self.set_axes_labels(gc, ix)
        #self.set_tick_labels(gc, ix)
        self.set_ticks(gc,ix)



    def set_ticks(self, gc, ix):
        """Set tick properties."""

        if hasattr(self, 'xspacing'):
            gc[ix].ticks.set_xspacing(self.xspacing)
        if hasattr(self, 'yspacing'):
            gc[ix].ticks.set_yspacing(self.yspacing)

        if hasattr(self, 'xminor_freq'):
            if hasattr(self, 'yminor_freq'):
                gc[ix].ticks.set_minor_frequency(self.xminor_freq,
                                                 self.yminor_freq)
            else :
                gc[ix].ticks.set_minor_frequency(self.xminor_freq)


        if hasattr(self, 'tick_color'):
            gc[ix].ticks.set_color(self.tick_color)

        if hasattr(self, 'tick_length'):
            gc[ix].ticks.set_length(self.tick_length)

        if hasattr(self, 'tick_linewidth'):
            gc[ix].ticks.set_linewidth(self.tick_linewidth)

        if hasattr(self, 'tick_direction'):
            gc[ix].ticks.set_tick_direction(self.tick_direction)


        if hasattr(self, 'tick_xhide'):
            if self.tick_xhide :
                gc[ix].ticks.hide_x()

        if hasattr(self, 'tick_yhide'):
            if self.tick_yhide :
                gc[ix].ticks.hide_y()
