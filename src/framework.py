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
import re
import sys

import aplpy
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import display as dsp
import markers as mrk

from astropy.io import fits
from astropy import wcs



class Figure(object):
    """Define object Figure."""
    
    def __init__(self, cnfg, dirs):

        if 'figure' in cnfg :
            fig_cfg= cnfg['figure']
            
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

        if 'frame' in cnfg:
            self.frame = Frame(cnfg['frame'], dirs)
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
        self.calc = None

        if 'font' in cnfg:
            self.fonts = self.read_font(cnfg['font'])
        else :
            self.fonts = None

        if 'yx' in cnfg:
            self.yx = cnfg['yx']
        else:
            self.yx = [1,1]

        self.gridsize = self.yx[0] * self.yx[1]
        self.gridpanels = {}
        for gr in range(self.gridsize) :
            nwstr = {'pstr' :  '' , 'data' : [] }
            self.gridpanels[gr] = nwstr

            
        if 'pad' in cnfg:
            self.pad = cnfg['pad']
        else :
            self.pad = None

        if 'margins' in cnfg:
            self.margins = cnfg['margins']
        else :
            self.margins = None
            
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
            self.contour = dsp.Contour(cnfg['contours'], self)

        if 'axes' in cnfg:
            self.axes = Axes(cnfg['axes'], self, fonts=self.fonts)

        if 'colorbar' in cnfg:
            self.colorbar = dsp.Colorbar(cnfg['colorbar'], self,
                                         fonts=self.fonts)

            
        self.wd = dirs['wkdir']


        data_str = [ d for d in cnfg if 'dataset' in d]

        if data_str : 
            self.datasets = {}
            self.load_generic_datasets(data_str, cnfg)


        panel_str = [k for k in cnfg if 'panel' in k]

        if panel_str :
            self.load_panels(panel_str,cnfg)


        panel_list = []

        for pnl in range(self.gridsize):

            if self.gridpanels[pnl]['pstr']:
                print("    + adding panel:", pnl, "...")
                panel = self.gridpanels[pnl]['pstr']
                panel_list.append(Panel(cnfg[panel], panel, pnl, self))

            elif self.gridpanels[pnl]['data']:
                print("    + adding panel:", pnl, "...")
                panel_list.append(Panel(None, 'panel_'+'g'+str(pnl), pnl, self))


        self.panels = panel_list
        
                

    def load_generic_datasets(self, d_str, cfg):
        """Reads the values of a generic dataset"""
        
        for g in d_str:

            if 'panels' in cfg[g]:
                p = cfg[g]['panels']

                lst_pan = []
                for it in p:
                    sit = str(it).split("-")
                    if len(sit) == 1:
                        pn = int(sit[0])

                        if pn > self.gridsize-1:
                            print("ERROR: dataset panels definition is out of",
                                  "bounds")
                            sys.exit(1)
                                
                        lst_pan.append(pn)

                    elif len(sit) == 2:
                        if not sit[0]:
                            print("ERROR: wrong dataset panels definition")
                            sys.exit(1)
                                
                        ini_rg = int(sit[0])
                        end_rg = int(sit[1])
                        if ini_rg > end_rg:
                            swap = ini_rg
                            ini_rg = end_rg
                            end_rg = swap

                        if end_rg > self.gridsize-1:
                            end_rg = self.gridsize - 1
                            
                        b = list( range(ini_rg, end_rg+1) )
                            
                        lst_pan.extend(b)

                    else:
                        print("ERROR: wrong format in dataset panels",
                              "definition")
                        sys.exit(1)
                    
                panel_set = list(set(lst_pan))

            else:
                panel_set = list( range(0, self.gridsize) )

            if 'slices' in cfg[g]:
                slc = cfg[g]['slices']
                
                sequences = {}

                for ix, sl in enumerate(slc):
                    rglist = []
                    
                    tp = tuple(sl)

                    for it in tp:

                        slrg = str(it).split("-")

                        if len(slrg) == 1 :
                            rglist.append(int(slrg[0]))

                        elif len(slrg) == 2:
                            if not slrg[0]:
                                print("ERROR: wrong dataset slices definition")
                                sys.exit(1)
                                
                            ini_sl = int(slrg[0])
                            end_sl = int(slrg[1])
                            if ini_sl > end_sl:
                                swap = ini_sl
                                ini_sl = end_sl
                                end_sl = swap

                            rglist.extend(list(range(ini_sl, end_sl+1)))

                        else:
                            print("ERROR: wrong format in dataset slices",
                                  "definition")

                    while len(rglist) < len(panel_set):
                        rglist.append(rglist[-1])

                    sequences[ix] = rglist
                        
                    
            arr_slc = np.transpose(np.array([sequences[0], sequences[1]]))

            for ct, panel in enumerate(panel_set):
                ds_str = g+'__'+str(panel)
            
                self.datasets[ds_str] = { k: cfg[g][k] for k in
                                          cfg[g].keys() - {'panels'} -
                                          {'slices'}}
                self.datasets[ds_str]['slices'] = list(arr_slc[ct])

                self.gridpanels[panel]['data'] = self.add_to_list(
                    self.gridpanels[panel]['data'], ds_str)



    def load_panels(self, p_str, cfg):
        """Load panels definition"""
        
        p_idx = 0

        p_order = self.read_panel_order(p_str, cfg)
        self.npanels = len(p_order)

        for ct, panel in enumerate(p_str):
            p_ord, p_idx = self.set_panel_order(p_order, ct, p_idx)

            if p_ord in self.gridpanels:
                self.gridpanels[p_ord]['pstr'] = panel



    @staticmethod
    def add_to_list(plist, new):

        nstr = plist
        nstr.append(new)
        return nstr

    

    @staticmethod
    def read_font(ff):

        dct = {}
        for prop in ['family', 'style', 'size', 'variant', 'stretch', 'weight']:
            if prop in ff :
                dct[prop] = ff[prop]
        return dct



    def read_panel_order(self, plist, cfg):
        """look for order key in panels config and add value to order list"""
        
        order_list = []
        
        for pnl in plist:
            if 'order' in cfg[pnl] :
                ord_val = cfg[pnl]['order']
                if ord_val < 0 or ord_val > self.gridsize-1:
                    print("ERROR: wrong panel order value")
                    sys.exit(1)
                order_list.append(ord_val)
            else:
                order_list.append(None)

        return order_list



    @staticmethod
    def set_panel_order(ord_id, i, idx):
        """get the order id of the panel"""
        
        if ord_id[i] != None :
            order = ord_id[i]
        else:
            while idx in ord_id:
                idx += 1
            order = idx
            idx += 1

        return order, idx

    

    def add_panels(self, fig):
        """Add panels to the figure.

        Argument:
            fig: handle of a FITSFigure instance
        """

        p_idx = 0
        gc = []


        if self.pad :
            fig.subplots_adjust(wspace=self.pad[0], hspace=self.pad[1])

        if self.margins:
            fig.subplots_adjust(left=self.margins[0], right=self.margins[1],
                                bottom=self.margins[2], top=self.margins[3])

        for p in self.panels :
            print("    + plotting panel:", p.name, "...")
            p.add_panel(fig, gc, p_idx)
            p_idx += 1

            
            
class Panel(object) :
    """ Create a panel."""
     
    def __init__(self, cnfg, name, idx, parent):

        self.name = name

        if cnfg != None and 'font' in cnfg:
            self.fonts = self.read_font(cnfg['font'], parent)
        elif hasattr(parent, 'fonts'):
            self.fonts = self.read_font(None, parent)
        else:
            self.fonts = None

        if cnfg != None and 'position' in cnfg:
            cpos = cnfg['position']
            if np.shape(cpos)[0] == 3:
                self.position = (cpos[0], cpos[1], cpos[2])
            elif np.shape(cpos)[0] == 4:
                self.position = [cpos[0], cpos[1], cpos[2], cpos[3]]
            else :
                print("  ERROR: wrong definition of panel position")
                sys.exit(1)
                
        else:
            self.position = (parent.yx[0], parent.yx[1], idx+1)

            
        if cnfg != None and 'order' in cnfg:
            self.order = cnfg['order']
        else:
            self.order = None
            
            
        if cnfg != None and 'view' in cnfg:
            self.view = View(cnfg['view'], parent)
        elif hasattr(parent, 'view'):
            self.view = View(None, parent)

        if cnfg != None and 'axes' in cnfg:
            self.axes = Axes(cnfg['axes'], parent, fonts=self.fonts, p_id=idx,
                             tc=(parent.npanels, parent.yx[1]))
        elif hasattr(parent, 'axes'):
            self.axes = Axes(None, parent, fonts=self.fonts, p_id=idx,
                             tc=(parent.npanels, parent.yx[1]))
        else:
            self.axes = Axes(None, None, fonts=self.fonts, p_id=idx,
                             tc=(parent.npanels, parent.yx[1]))
            
        if cnfg != None and 'labels' in cnfg:
            self.labels = mrk.Label(cnfg['labels'], parent, fonts=self.fonts)
        elif hasattr(parent, 'labels'):
            self.labels = mrk.Label(None, parent, fonts=self.fonts)
        else:
            self.labels = None

        if cnfg != None and 'markers' in cnfg:
            self.markers = mrk.Marker(cnfg['markers'], parent, fonts=self.fonts)
        elif hasattr(parent, 'markers'):
            self.markers = mrk.Marker(None, parent, fonts=self.fonts)
        else:
            self.markers = None

        if cnfg != None and 'colorbar' in cnfg:
            self.colorbar = dsp.Colorbar(cnfg['colorbar'], parent,
                                         fonts=self.fonts)
        elif hasattr(parent, 'colorbar'):
            self.colorbar = dsp.Colorbar(None, parent, fonts=self.fonts)
        else:
            self.colorbar = None
            
        dataset_list = []

        plist = parent.gridpanels[idx]['data']
        if plist:
            for d in plist:
                dataset_list.append(dsp.Dataset(parent.datasets[d], d, parent))
                

        if cnfg != None:
            dataset_str = [k for k in cnfg if 'dataset' in k]

            if dataset_str :
                for d in dataset_str:
                    dataset_list.append(dsp.Dataset(cnfg[d], d, parent))

                
        self.datasets = dataset_list



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

                # TODO:
                # check a different center
                # if offsets, calc new center
                #  below has to be separated
                if vw.center == None:
                    vw.center = d.get_reference()
                    #vw.center = d.center

                refpos = vw.center
                if vw.coords == 'off' :
                    hdu, ra, dec = d.to_offsets(vw.center)
                    vw.center = (ra, dec)
                else :
                    hdu = d.filename


                if d.slices :
                    hdu = d.get_slice(hdu=hdu)

                # hdu ==> d.hdu after doing all conversions

                gc.append(aplpy.FITSFigure(hdu, figure=fig,
                                           subplot=self.position,
                                           dimensions=d.dims))
                

                vw.set_view(gc[idx])

            d.show(gc[idx], coords=vw.coords, ref=refpos)


            # adds channel label
            try:
                if d.lblchan:
                    #d.ds_label.label_props['text'] = d.lbltxt
                    self.labels.label_list.append(d.ds_label)
                    
            except AttributeError:
                pass
            
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
            if self.colorbar != None :
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
                self.box = (self.read_units(cnfg['box'][0]),
                            self.read_units(cnfg['box'][1]))
                      
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

                
        if cnfg != None and 'coords' in cnfg:
            self.coords = cnfg['coords']
        else :
            try:
                self.coords = parent.view.coords
            except: 
                self.coords = 'abs'


                

            

    def set_view(self, g):
        """set the view of the panel."""

        if self.vtype == 'radius' :
            if self.coords == 'off':
                ang = coord.Angle(self.radius * u.deg)
                self.radius = (ang.to(u.arcsec)).arcsec
            g.recenter(self.center[0], self.center[1], radius=self.radius)

        elif self.vtype == 'box' :
            #TODO
            # change and test box for offsets
            #
            g.recenter(self.center[0], self.center[1],
                       height=self.box[0], width=self.box[1])



    @staticmethod
    def read_units(val):
        """Proces a configuration file value containing an angle unit 
           string.
        """
        
        strval = str(val)

        if re.match("^.*[a-z]+.*$", val):
            try:
                return (coord.Angle(val).to(u.deg)).degree
            except ValueError: 
                print("ERROR: unrecognized definition of units:", val)
                sys.exit(1)
        else:
            return (coord.Angle(val+'deg').to(u.deg)).degree



class Axes(object):
    """Class that contains the definition of the plot axes."""
    
    def __init__(self, cfg, parent, fonts=None, p_id=-1, tc=(1,1)):

        
        if p_id > -1 :
            self.set_axes_in_grid(p_id, tc, parent)

            
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



    def set_axes_in_grid(self, pid, t, parent):
        """Hides axis in grid when no padding"""
        
        if (pid % t[1]) != 0 and parent.pad[0] < 0.16:
            self.axis_yhide = True
            self.ticklabel_yhide = True 

        if (t[0] - pid) > t[1] and parent.pad[1] < 0.12:
            self.axis_xhide = True
            self.ticklabel_xhide = True



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

        
