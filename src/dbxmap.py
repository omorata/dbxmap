#!/usr/bin/env python3
##
##  dbxmap.py
##
##  O. Morata 2019-2020
##
##   Tool to plot maps using Aplpy
##
##
import argparse

import matplotlib.pyplot as plt
import aplpy

import numpy as np
from astropy import wcs
from astropy.io import ascii, fits
import astropy.coordinates as coord
from astropy import units as u

import yaml

import os
import sys

import matplotlib


# allow use of LaTex in the texts
#
matplotlib.rcParams['text.usetex'] = True

##-- Class definitions -------------------------------------------------

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
            self.markers = Marker(cnfg['markers'], self)

        if 'labels' in cnfg:
            self.labels = Label(cnfg['labels'], self)
        else :
            self.labels = None

        if 'pixrange' in cnfg:
            self.pixrange = Pixrange(cnfg['pixrange'], self)

        if 'contours' in cnfg:
            self.contour = Contour(cnfg['contours'], self)

            
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
            
            
        if hasattr(parent, 'labels') and parent.labels != None:
            self.labels = parent.labels
            self.label_props = parent.labels.label_props.copy()
        else :
            self.labels = None

        if 'labels' in cnfg:
            self.labels = Label(cnfg['labels'], self)

            
        if 'markers' in cnfg:
            self.markers = Marker(cnfg['markers'], parent)
        elif hasattr(parent, 'markers'):
            self.markers = Marker(None, parent)

        if 'colorbar' in cnfg:
            self.colorbar = Colorbar(cnfg['colorbar'], self)
            
        dataset_list = []

        dataset_str = [k for k in cnfg if 'dataset' in k]

        if dataset_str :
            for d in dataset_str:
                dataset_list.append(Dataset(cnfg[d], d, parent))
                
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

        if self.labels != None :
            for lb in self.labels.label_list :
                gc = lb.add_label(gc, idx)

        if hasattr(self, 'markers') :
            a = self.markers
            a.add_markers(gc, idx)

            
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


            
class Dataset(object) :
    """Create a dataset object."""
    
    def __init__(self, cnfg, name, parent):

        self.name = name

        
        if 'filename' in cnfg:
            self.filename = os.path.join(parent.wd, cnfg['filename'])

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

                if not hasattr(self.pixrange, 'range'):
                    if hasattr(self.pixrange, 'scale') :
                        self.pixrange.range = self.get_range_from_scale()
                    else :
                        raise Exception('Error: undefined pixel range')
                        sys.exit(1)
                    
            elif self.dtype == 'cntr' or self.dtype == 'contour':
                if 'contours' in cnfg:
                    self.contour = Contour(cnfg['contours'], parent,
                                           name=self.name)
                elif hasattr(parent, 'contour'):
                    self.contour = Contour(None, parent, name=self.name)
                else :
                    raise Exception('Error: contour is not defined anywhere')
                
        if 'beam' in cnfg:
            self.read_beam_parameters(dict(cnfg['beam']))


            
    def get_range_from_scale(self):
        """Calculates the pixel range from the given percentile scale."""
        
        hdulist = fits.open(self.filename)
        dt = hdulist[0].data

        diff = (100-self.pixrange.scale)*0.5
        rg = [np.nanpercentile(dt, diff), np.nanpercentile(dt, 100-diff)]

        return rg


    
    def get_reference(self):
        """ read the reference position of the datasetfrom the FITS
            header
        """

        hdulist = fits.open(self.filename)
        w = wcs.WCS(hdulist[0].header)

        center = (w.wcs.crval[0], w.wcs.crval[1])

        return center
    

    
    def read_beam_parameters(self, battr):
        """Reads the beam parameters from the configuration file

        It populates self.beam_shape and self.beam_args

        TODO:
            proper unit treatment
        """
        
        if 'bmaj' in battr :
            bmaj = battr['bmaj'] / 3600
            self.beam_shape = {'major' : bmaj}
            
        if 'major' in self.beam_shape :

            if 'bmin' in battr :
                self.beam_shape['minor'] = battr['bmin'] / 3600
            else :
                self.beam_shape['minor'] = bmaj
            
            if 'bpa' in battr :
                self.beam_shape['angle'] = battr['bpa']
            else :
                self.beam_shape['angle'] = 0.
                               
        if 'linewidth' in battr :
            linewidth = float(battr['linewidth'])
        else :
            linewidth = 0.5

        if 'linestyle' in battr :
            linestyle = battr['linestyle']
        else :
            linestyle = 'solid'

        if 'edgecolor' in battr :
            edgecolor = battr['edgecolor']
        else :
            edgecolor = 'steelblue'

        if 'facecolor' in battr :
            facecolor = battr['facecolor']
        else :
            facecolor = 'steelblue'

        if 'alpha' in battr :
            alpha = float(battr['alpha'])
        else :
            alpha = 1

        if 'frame' in battr:
            frame = battr['frame']
        else:
            frame = False

        self.beam_args = {'linewidth' : linewidth, 'linestyle' : linestyle,
                          'edgecolor' : edgecolor, 'alpha' : alpha,
                          'facecolor' : facecolor, 'frame' : frame}

        if 'corner' in battr:
            self.beam_args['corner'] = battr['corner']

        if 'borderpad' in battr:
            self.beam_args['borderpad'] = battr['borderpad']
                
        if 'pad' in battr:
            self.beam_args['pad'] = battr['pad']
                

            
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
                                aspect='equal')


        
    def show_contour(self, gc, idx) :
        """ show a contour map of the dataset
        """
        gc[idx].show_contour(self.filename,
                             levels=self.contour.levels,
                             colors=self.contour.colors,
                             linewidths=self.contour.linewidth,
                             linestyles=self.contour.linestyle)
        


        
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
            r_dim = np.shape(cfgrange)[0]

            if '%' in str(cfgrange[0]):
                self.scale = float(cnfg['range'][0].split('%')[0])

                if self.scale > 100 or self.scale < 1 :
                    print("ERROR: wrong scale for the pixel range")
                    sys.exit(1)

                if r_dim == 2 :
                    self.stretch = cfgrange[1]
                else :
                    self.stretch = 'linear'

            else :
                if r_dim > 1 :
                    self.range = [(float(i) * self.base) for i in cfgrange[0:2]]
                    
                    if r_dim == 3 :
                        self.stretch = cfgrange[2]
                    else :
                        self.stretch = 'linear'
                    
                else :
                    print("ERROR: wrong number of parameters in pixrange",
                          "definition")
                    sys.exit(1)

        else :
            try:
                self.range = parent.pixrange.range

            except AttributeError:
                try :
                    self.scale = parent.pixrange.scale

                except AttributeError:
                    print("We are in trouble")
                    pass

            try:
                self.stretch = parent.pixrange.stretch
            except AttributeError:
                pass

            

            
class Colorbar(object):
    """creat colorbar"""

    def __init__(self, cnfg, parent) :
        """Initialization of a colorbar object."""
        
        self.add = 'y'
        
        if cnfg != None and 'width' in cnfg:
            self.width = float(cnfg['width'])
        else :
            try:
                self.width = parent.colorbar.width

            except AttributeError:
                self.width = 0.3
            
        if cnfg != None and 'location' in cnfg:
            self.location = cnfg['location']
        else :
            try:
                self.location = parent.colorbar.location

            except AttributeError:
                self.location = 'right'

            
        if cnfg != None and 'text' in cnfg:
            self.text = cnfg['text']
        else :
            try:
                self.text = parent.colorbar.text

            except AttributeError:
                self.text = ''


                
class Contour(object):
    """ Create contours.
    """

    tol = 1.e-10
    
    def __init__(self, cnfg, parent, name='') :

        self.append_contours = False
        
        if cnfg !=None and 'append' in cnfg:
            if cnfg['append'] == True :
                self.append_contours = True
                try :
                    self.levels = parent.contour.levels.copy()
                except:
                    pass

                
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
            self.add_levels(self.scale_values(self.base, cnfg['levels']))

            
        if cnfg != None and 'gen_levels' in cnfg:
            self.generate_levels(self.scale_values(self.base,
                                                   cnfg['gen_levels']))
            
                
        if cnfg != None and 'pbase' in cnfg:
            self.pbase = float(cnfg['pbase'])
        else :
            try:
                self.pbase = parent.contour.pbase

            except AttributeError:
                pass

            
        if cnfg != None and 'plevels' in cnfg:
            try :
                self.add_levels(self.scale_values(self.pbase*0.01,
                                                  cnfg['plevels']))

            except AttributeError :
                print("WARNING: no base defined for percentage levels")
                pass

            
        if cnfg != None and 'gen_plevels' in cnfg:
            try :
                self.generate_levels(self.scale_values(self.pbase*0.01,
                                                       cnfg['gen_plevels']))
                
            except AttributeError :
                print("WARNING: no base defined for percentage levels")
                pass

        if cnfg != None and 'flevels' in cnfg:
            self.add_levels(cnfg['flevels'])

            
        if cnfg != None and 'gen_flevels' in cnfg:
            self.generate_levels(cnfg['gen_flevels'])
                
            
        if cnfg != None and 'linewidth' in cnfg:
            self.linewidth = float(cnfg['linewidth'])
        else :
            try:
                self.linewidth = parent.contour.linewidth

            except AttributeError:
                self.linewidth = 1.
                
        if cnfg != None and 'linestyle' in cnfg:
            self.linestyle = str(cnfg['linestyle'])
        else :
            try:
                self.linestyle = parent.contour.linestyle

            except AttributeError:
                self.linestyle = '-'


        if name:
            if not hasattr(self, 'levels'):
                try :
                    self.levels = parent.contour.levels.copy()
                except:
                    print("ERROR: no level definition anywhere for the dataset",
                             name)
                    sys.exit("Exiting")

                    
    @staticmethod
    def scale_values(factor, vlist):
        """Scales the values in vlist by factor."""

        return [(float(i) * factor) for i in vlist]

        
    def add_levels(self, new_levels) :
        """Adds new levels to the current list of levels.

        If the list of levels already existed, it sorts them.
        Args:
           new_levels: list of new levels to add
        """
        
        if not hasattr(self, 'levels') :
            self.levels = new_levels

        else:
            prev = self.levels
            prev.extend(new_levels)
            prev.sort()
            self.levels = prev


    def generate_levels(self, lev_list):
        """Automatically creates a set of levels and adds them to the
        list of levels

        Args:
            lev_list - three element list with the minimum and maximum
                       values of the range, plus the increment
        """

        if len(lev_list) != 3 :
            print("ERROR: wrong number of elements to generate levels")
            sys.exit(1)
            
        lv = lev_list[0]
            
        gen_levs = []
        while np.sign(lev_list[2]) * (lev_list[1] - lv) > -self.tol :
            gen_levs.append(lv)
            lv += lev_list[2]

        self.add_levels(gen_levs)
            

        
class Label(object):
    """Create a label."""

    def __init__(self, cfg, parent):

        if hasattr(parent, 'labels') and parent.labels != None :
            self.label_list = parent.labels.label_list.copy()
        else:
            self.label_list = []

        if hasattr(parent, 'label_props') :
            self.label_props = parent.label_props.copy()
        else :
            self.label_props = self.default_label_props()

        property_list = ['text', 'relative', 'position', 'color', 'size',
                         'style']

        if cfg != None:
            for prop in property_list :
                if prop in cfg:
                    self.label_props[prop] = cfg[prop]

            file_str = [k for k in cfg if 'label' in k]

            if file_str :
                for lab in file_str:
                    self.label_list.append(Label(cfg[lab], self))
        else :
            print("\n  +++ WARNING: naked 'label' field in the configuration",
                  "file +++\n")

            

    def add_label(self, pp, idx) :
        """Adds label to the panel."""

        pp[idx].add_label(self.label_props['position'][0],
                          self.label_props['position'][1],
                          **self.label_props)

        return pp



    def default_label_props(self):
        """Define default values for label properties."""

        props = { 'color' : 'black',
                  'relative' : False,
                  'position' : [0.1, 0.9],
                  'size' : 12,
                  'style' :  'normal',
                  'text' : "" }
        return props




class Marker(object) :
    """Class to define markers (including polygons)."""

    
    
    properties = ['type', 'edgecolor', 'linewidth', 'linestyle',
                  'facecolor', 'show_label', 'color', 'weight', 'size']
    defaults = ['', 'black', 1.0, 'solid', 'none', False, 'black',
                'normal', 12]

    
    def __init__ (self, cfg, parent):

        self.wkdir = parent.wkdir
        
        if hasattr(parent, 'markers'):
            self.marklist = parent.markers.marklist.copy()
            #list_markers = self.marklist
        else:
            #list_markers = []
            self.marklist = []

        if hasattr(parent, 'marker_props') :
            self.marker_props = parent.marker_props.copy()
        else :
            self.marker_props = self.default_marker_props()

            
        if cfg != None :
            for prop in self.properties :
                if prop in cfg:
                    self.marker_props[prop] = cfg[prop]


            if 'file' in cfg:
                fname = os.path.join(self.wkdir, cfg['file'])

                #self.marklist.extend(self.read_markers(fname))
                self.marker = self.read_markers(fname))

            marker_str = [k for k in cfg if 'marker' in k]

            if marker_str :
                #if not hasattr(self, 'marklist'):
                #    self.marklist = []
                    
                for mrk in marker_str :
                    #fname = os.path.join(self.wkdir, cfg[files])

                    #self.marklist.extend(self.read_markers(fname))
                    self.marklist.append(Marker(cfg[mrk], self))

            #else :
            #    self.marklist = None
                

                
    def default_marker_props(self):
        """Define default values for marker properties."""

        props = { key:value for key, value in
                  zip(self.properties, self.defaults))
        return props


    
    def add_markers(self, gc, i):

        lmarkers = self.marklist
        
        for mk in lmarkers:
            type = mk['type']

            if type == 'polygon' :
                gc[i].show_polygons(mk['corners'], **mk['style'])
                
                gc[i].add_label(mk['bcenter'][0], mk['bcenter'][1],
                                mk['id'], **mk['l_style'])
                
            elif type == 'cross' :
                gc[i].show_markers(mk[0]['x'].degree, mk[0]['y'].degree,
                                   edgecolor=mk[0]['color'],
                                   c=mk[0]['filled'],
                                   linewidths=mk[0]['linewidth'],
                                   s=mk[0]['size'],
                                   marker='$+$')



    def read_markers(self, fname):
        """Reads a markers file."""

        tbl = ascii.read(fname, delimiter=" ", format="basic")
        
        mark_list = []
        for marker in tbl:

            #if marker[0] == "cross" :
            #    props = Marker.check_cross(marker)
            #    if props != None :
            #        mark_list.append(props)

            if self.marker_props['type'] == 'polygon' :
                if marker['type'] == 'Polygon' :
                    properties = self.read_polygon(marker)
                    mark_list.append(properties)
                else :
                    print("\n  +++ ERROR: marker type polygon reading a",
                          "file without Polygons\n")
                    sys.exit(1)

        return mark_list

    
    
    def read_polygon(self, it):
        """Read out the definition of a polygon marker."""

        attrib = {'type' : 'polygon', 'id' : it['id']}
        
        coord_elements = it['corners'].split(" ")

        c_dim = np.shape(coord_elements)[0]
        if c_dim % 2 != 0 :
            print("ERROR: uneven coordinates for polygon", it['id'])
            sys.exit(1)

        array_corners = np.empty([int(c_dim *0.5), 2])

        lcorners = []
        for p in range(int(c_dim * 0.5)) :
            array_corners[p,:] = [float(coord_elements.pop(0)),
                                  float(coord_elements.pop(0))]

        lcorners.append(array_corners)
        attrib['corners'] = lcorners

        # baricenter
        #
        attrib['bcenter'] = [np.average(array_corners[:,0]),
                             np.average(array_corners[:,1])]

        style = {}
        for ss in ['edgecolor', 'linewidth', 'linestyle', 'facecolor']:
            style[ss] =  self.marker_props[ss]
        #if it['edgecolor'] != "" :
        #    edgecolor = it['edgecolor']
        #else :
        #    edgecolor = 'black'

        #if it['linewidth'] != "" :
        #    linewidth = it['linewidth']
        #else :
        #    linewidth = 1.0

        #if it['linestyle'] != "" :
        #    linestyle = it['linestyle']
        #else :
        #    linestyle = 'solid'

        #style = { 'edgecolor': edgecolor, 'linewidth' : linewidth,
        #          'linestyle' : linestyle}
        
        attrib['style'] = style

        #if it['labelcolor'] != "":
        #    labelcolor = it['labelcolor']
        #else :
        #    labelcolor = 'black'

        #if it['weight'] != "":
        #    weight = it['weight']
        #else :
        #    weight = 'normal'
            
        #if it['fontsize'] != "":
        #    size = it['fontsize']
        #else :
        #    size = 12
        labelstyle = {}
        for ls in ['color', 'weight', 'size']:
            labelstyle[ls] =  self.marker_props[ls]
            
        #labelstyle = {'color' : labelcolor, 'weight' : weight, 'size' : size}

        attrib['l_style'] = labelstyle
        
        return attrib
        
    
    @staticmethod
    def check_cross(it) :
        """ check configuration of a cross marker
        """
        if len(it) != 8 :
            print("Wrong number of elements in cross marker")
            return None
        else :
            attribs = {'sym' : (4,1,0),
                       'size': float(it[4]),
                       'linewidth' : float(it[5]),
                       'color' : it[6]
            }
            if it[1] == "abs" :
                attribs["x"] = coord.Angle(it[2], unit=u.hour)
                attribs["y"] = coord.Angle(it[3], unit=u.degree)

            if it[7] == "y" :
                attribs['filled'] = it[6]
            else :
                attribs['filled'] = 'none'
                
            return attribs

        
##-- End of class definitions ------------------------------------------
        
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
        
    G = Figure(cfg, dirs)

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


##-- End of functions --------------------------------------------------

##-- Main --------------------------------------------------------------

if __name__ == "__main__" :
    ## defaults
    #
    #config_file = 'plot_cfg.yml'

    print(" Starting...")

    matplotlib.use('Agg')
    
    args = read_command_line()

    dirs = {'wkdir' : args.work_dir, 'outdir' : args.output_dir}
    files = {'config' : args.config}

    if args.log is not None:
        files['log'] = args.log

    # process the configuration
    #
    fig = process_config(files, dirs)

    # plot the figure
    #
    fig.create()

    fig.f_print()
    
    fig.end()
    
##-- End of main -------------------------------------------------------
