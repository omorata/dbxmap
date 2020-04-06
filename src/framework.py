#!/usr/bin/env python3
##
## framework.py
##
## O. Morata 2020
##
##  module containing the classes that build the basic structure of
##  the map
##

import copy
import os
import sys

import aplpy
import matplotlib.pyplot as plt
import numpy as np

import display as dsp
import markers as mrk
import astropy.units as u

import re

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
    """ Class defining the elements of the figure."""

    def __init__(self, cnfg, dirs):

        self.wkdir = dirs['wkdir']

        if 'font' in cnfg:
            self.fonts = self.read_font(cnfg['font'])
        else :
            self.fonts = None

        if 'yx' in cnfg:
            self.yx = cnfg['yx']
        else:
            self.yx = [1,1]

        if 'view' in cnfg:
            self.view = View(cnfg['view'], self)

        if 'markers' in cnfg:
            self.markers = mrk.Marker(cnfg['markers'], self, fonts=self.fonts)

        if 'labels' in cnfg:
            self.labels = mrk.Label(cnfg['labels'], self, fonts=self.fonts)
        else :
            self.labels = None

        if 'pixrange' in cnfg:
            self.pixrange = dsp.Pixrange(cnfg['pixrange'], self)

        if 'contours' in cnfg:
            self.contour = dsp. Contour(cnfg['contours'], self)

        if 'axes' in cnfg:
            self.axes = Axes(cnfg['axes'], self, fonts=self.fonts)

        if 'colorbar' in cnfg:
            self.colorbar = dsp.Colorbar(cnfg['colorbar'], self,
                                         fonts=self.fonts)

            
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



    @staticmethod
    def read_font(ff):

        dct = {}
        for prop in ['family', 'style', 'size', 'variant', 'stretch', 'weight']:
            if prop in ff :
                dct[prop] = ff[prop]
        return dct



    def add_panels(self, fig):
        """Add panels to the figure.

        Argument:
            fig: handle of a FITSFigure instance
        """

        p_idx = 0
        gc = []
        
        for p in self.panels :
            print("    + plotting panel:", p.name, "...")
            p.add_panel(fig, gc, p_idx)
            p_idx += 1


            
            
class Panel(object) :
    """ Create a panel."""
     
    def __init__(self, cnfg, name, idx, parent):

        self.name = name

        if 'font' in cnfg:
            self.fonts = self.read_font(cnfg['font'], parent)
        elif hasattr(parent, 'fonts'):
            self.fonts = self.read_font(None, parent)
        else:
            self.fonts = None

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
            self.axes = Axes(cnfg['axes'], parent, fonts=self.fonts)
        elif hasattr(parent, 'axes'):
            self.axes = Axes(None, parent, fonts=self.fonts)
        else:
            self.axes = None
            
        if 'labels' in cnfg:
            self.labels = mrk.Label(cnfg['labels'], parent, fonts=self.fonts)
        elif hasattr(parent, 'markers'):
            self.labels = mrk.Label(None, parent, fonts=self.fonts)

        if 'markers' in cnfg:
            self.markers = mrk.Marker(cnfg['markers'], parent, fonts=self.fonts)
        elif hasattr(parent, 'markers'):
            self.markers = mrk.Marker(None, parent, fonts=self.fonts)

        if 'colorbar' in cnfg:
            self.colorbar = dsp.Colorbar(cnfg['colorbar'], parent,
                                         fonts=self.fonts)
        elif hasattr(parent, 'colorbar'):
            self.colorbar = dsp.Colorbar(None, parent, fonts=self.fonts)
            
        dataset_list = []

        dataset_str = [k for k in cnfg if 'dataset' in k]

        if dataset_str :
            for d in dataset_str:
                dataset_list.append(dsp.Dataset(cnfg[d], d, parent))
                
            self.datasets = dataset_list
        else:
            self.datasets = []


            
    @staticmethod
    def read_font(ff, parent):

        if hasattr(parent, 'fonts'):
            dct = copy.deepcopy(parent.fonts)
        else :
            dct = {}

        if ff != None:
            for prop in ['family', 'style', 'size', 'variant', 'stretch',
                         'weight']:
                if prop in ff :
                    dct[prop] = ff[prop]
        return dct



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

                vw.set_view(gc[idx])

            d.show(gc[idx])

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
            self.axes.set_axes(gc[idx])

        if self.labels != None :
            for lb in self.labels.label_list :
                lb.add_label(gc[idx])

        if self.markers != None:
            for mk in self.markers.marklist :
                mk.add_markers(gc[idx])

        if hasattr(self, 'colorbar'):
            self.colorbar.set_colorbar(gc[idx])

            

    def set_pixel_first(self):
        """Put the dataset with a pixel range in the first position."""

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
    """Create a view."""

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
                self.radius = self.read_units(cnfg['radius'] )

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

            

    def set_view(self, g):
        """set the view of the panel."""

        if self.vtype == 'radius' :
            g.recenter(self.center[0], self.center[1], radius=self.radius)

        elif self.vtype == 'box' :
            g.recenter(self.center[0], self.center[1],
                       height=self.box[0], width=self.box[1])



    @staticmethod
    def read_units(val):
        a = str(val)
        if re.match('^[0-9\.]*$', a) or re.match('^[0-9]*$',a) :
            return float(val)
        
        elif re.match('.*arcmin$', a):
            new = a.split('arcmin')[0] * u.arcmin
            return (new.to(u.degree))
        elif re.match('.*arcsec$', a):
            new = a.split('arcsec')[0] * u.arcsec
            return (new.to(u.degree))
        elif re.match('.*deg$', a):
            return a.split('deg')[0]
        else:
            print("ERROR: unrecognized units")
            sys.exit(1)




class Axes(object):
    """Class that contains the definition of the plot axes."""
    
    def __init__(self, cfg, parent, fonts=None):

        if cfg != None and 'axes_labels' in cfg:
            if hasattr(parent, 'axes'):
                self.read_axes_labels(cfg['axes_labels'], parent.axes,
                                      fonts=fonts)
            else :
                self.read_axes_labels(cfg['axes_labels'], None, fonts=fonts)
        else:
            try:
                self.read_axes_labels(None, parent.axes, fonts=fonts)
            except AttributeError:
                pass

        if cfg != None and 'tick_labels' in cfg:
            if hasattr(parent, 'axes'):
                self.read_tick_labels(cfg['tick_labels'], parent.axes,
                                      fonts=fonts)
            else :
                self.read_tick_labels(cfg['tick_labels'], None, fonts=fonts)
        else:
            try:
                self.read_tick_labels(None, parent.axes, fonts=fonts)
            except AttributeError:
                pass

        if cfg != None and 'ticks' in cfg :
            if hasattr(parent, 'axes'):
                self.read_ticks(cfg['ticks'], parent.axes)
            else :
                self.read_ticks(cfg['ticks'], None)
        else:
            try:
                self.read_ticks(None, parent.axes)
            except AttributeError:
                pass



    def read_axes_labels(self, ax, parent, fonts=None):
        """Read the configuration labels for the axes labels.""" 

        for at in ['xposition', 'yposition', 'xpad', 'ypad', 'xtext', 'ytext']:
            if ax != None and at in ax:
                setattr(self, at, ax[at])
            elif hasattr(parent, at):
                setattr(self, at, getattr(parent, at))

        if hasattr(parent, 'axis_font'):
            self.axis_font = copy.deepcopy(parent.axis_font)
        elif fonts != None :
            self.axis_font = copy.deepcopy(fonts)
        else :
            self.axis_font = {}

        if ax != None and 'font' in ax:
            ff = ax['font']

            for prop in ['family', 'style', 'size', 'variant', 'stretch',
                         'weight']:
                if prop in ff :
                    self.axis_font[prop] = ff[prop]

        if ax != None and 'hide' in ax:
            self.axis_xhide = ax['hide']
            self.axis_yhide = ax['hide']

        for at in ['xhide', 'yhide']:
            if ax != None and at in ax :
                setattr(self, 'axis_'+at, ax[at])
            elif hasattr(parent, 'axis_'+at):
                setattr(self, 'axis_'+at, getattr(parent, 'axis_'+at))



    def read_ticks (self, tx, parent):
        """Read the configuration for ticks."""

        for at in ['xspacing', 'xminor_freq', 'yspacing', 'yminor_freq']:
            if tx != None and at in tx:
                setattr(self, at, tx[at])
            elif hasattr(parent, at):
                setattr(self, at, getattr(parent, at))
                
        if tx != None and 'hide' in tx:
            self.tick_xhide = tx['hide']
            self.tick_yhide = tx['hide']

        for at in ['color', 'length', 'linewidth', 'direction', 'xhide',
                   'yhide']:
            if tx != None and at in tx:
                setattr(self, 'tick_'+at, tx[at])
            elif hasattr(parent, 'tick_'+at):
                setattr(self, 'tick_'+at, getattr(parent, 'tick_'+at))
                

            
    def read_tick_labels(self, tl, parent, fonts=None):
        """Read the configuration for the tick labels.""" 

        if tl != None and 'hide' in tl:
            self.ticklabel_xhide = tl['hide']
            self.ticklabel_yhide = tl['hide']

        for at in ['xposition', 'yposition', 'xformat', 'yformat', 'style',
                   'xhide', 'yhide']:
            if tl != None and at in tl:
                setattr(self, 'ticklabel_'+at, tl[at])
            elif hasattr(parent, 'ticklabel_'+at):
                setattr(self, 'ticklabel_'+at, getattr(parent, 'ticklabel_'+at))

        if hasattr(parent, 'ticklabel_font'):
            self.ticklabel_font = copy.deepcopy(parent.ticklabel_font)
        elif fonts != None:
            self.ticklabel_font = copy.deepcopy(fonts)
        else :
            self.ticklabel_font = {}

        if tl != None and 'font' in tl:
            ff = tl['font']

            for prop in ['family', 'style', 'size', 'variant', 'stretch',
                         'weight']:
                if prop in ff :
                    self.ticklabel_font[prop] = ff[prop]



    def set_axes(self, g):
        """Sets axes properties."""

        self.set_axes_labels(g.axis_labels)
        self.set_tick_labels(g.tick_labels)
        self.set_ticks(g.ticks)



    def set_axes_labels(self, gp) :
        """Set axes labels."""

        if hasattr(self, 'axis_font'):
            gp.set_font(**self.axis_font)

        if hasattr(self, 'xposition'):
            gp.set_xposition(self.xposition)

        if hasattr(self, 'yposition'):
            gp.set_yposition(self.yposition)

        if hasattr(self, 'xpad'):
            gp.set_xpad(self.xpad)

        if hasattr(self, 'ypad'):
            gp.set_ypad(self.ypad)

        if hasattr(self, 'xtext'):
            gp.set_xtext(self.xtext)

        if hasattr(self, 'ytext'):
            gp.set_ytext(self.ytext)

        if hasattr(self, 'axis_xhide'):
            if self.axis_xhide:
                gp.hide_x()

        if hasattr(self, 'axis_yhide'):
            if self.axis_yhide:
                gp.hide_y()



    def set_ticks(self, gp):
        """Set tick properties."""

        if hasattr(self, 'xspacing'):
            gp.set_xspacing(self.xspacing)
        if hasattr(self, 'yspacing'):
            gp.set_yspacing(self.yspacing)

        if hasattr(self, 'xminor_freq'):
            if hasattr(self, 'yminor_freq'):
                gp.set_minor_frequency(self.xminor_freq,
                                                 self.yminor_freq)
            else :
                gp.set_minor_frequency(self.xminor_freq)

        if hasattr(self, 'tick_color'):
            gp.set_color(self.tick_color)

        if hasattr(self, 'tick_length'):
            gp.set_length(self.tick_length)

        if hasattr(self, 'tick_linewidth'):
            gp.set_linewidth(self.tick_linewidth)

        if hasattr(self, 'tick_direction'):
            gp.set_tick_direction(self.tick_direction)

        if hasattr(self, 'tick_xhide'):
            if self.tick_xhide :
                gp.hide_x()

        if hasattr(self, 'tick_yhide'):
            if self.tick_yhide :
                gp.hide_y()



    def set_tick_labels(self, gp) :
        """Set the tick labels."""

        if hasattr(self, 'ticklabel_font'):
            gp.set_font(**self.ticklabel_font)

        if hasattr(self, 'ticklabel_xposition'):
            gp.set_xposition(self.ticklabel_xposition)

        if hasattr(self, 'ticklabel_yposition'):
            gp.set_yposition(self.ticklabel_yposition)

        if hasattr(self, 'ticklabel_xformat'):
            gp.set_xformat(self.ticklabel_xformat)

        if hasattr(self, 'ticklabel_yformat'):
            gp.set_yformat(self.ticklabel_yformat)

        if hasattr(self, 'ticklabel_style'):
            gp.set_style(self.ticklabel_style)

        if hasattr(self, 'ticklabel_xhide'):
            if self.ticklabel_xhide:
                gp.hide_x()

        if hasattr(self, 'ticklabel_yhide'):
            if self.ticklabel_yhide:
                gp.hide_y()

            


