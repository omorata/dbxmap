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

import matplotlib


# allow use of LaTex in the texts
#
matplotlib.rcParams['text.usetex'] = True

##-- Class definitions -------------------------------------------------

class Page(object):
    """ define object Page
    """ 
    def __init__(self, cnfg, dirs):
        if 'outfile' in cnfg:
            self.outfile = os.path.join(dirs['outdir'], cnfg['outfile'])
        else :
            self.outfile = os.path.join(dirs['outdir'], 'outfile.pdf')

        if 'dpi' in cnfg :
            self.dpi = cnfg['dpi']
        else:
            self.dpi = 100

        if 'size' in cnfg:
            self.size = cnfg['size']
        else :
            self.size = [10,10]

        if 'name' in cnfg:
            self.name = cnfg['name']


            
    def create(self) :
        """ create figure
        """
        self.f = plt.figure(figsize=(self.size))
        return self.f



    def f_print(self):
        """ print figure
        """
        print("  + plotting output file:", self.outfile, "...")
        self.f.savefig(self.outfile, dpi=self.dpi)



    def end(self):
        """ ending comment
        """
        if hasattr(self, 'name') :
            print("\n  ... figure <<", self.name, ">> done!\n")
        else :
            print("\n  ... done!\n")
    


        
class DbxFig(object):
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
            self.markers = Markers(cnfg['markers'], self)

        #if 'labels' in cnfg:
        #    self.labels = Label(cnfg['labels'], self)

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
                print("  + adding panel:", panel, "...")
                panel_list.append(Panel(cnfg[panel], panel, p_idx, self))
                p_idx += 1

            self.panels = panel_list
        
        else:
            self.panels = []



    def add_panels(self, fig):

        p_idx = 0
        gc = []
        
        for p in self.panels :
            print("  + plotting panel:", p.name, "...")
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
            
            
        if 'labels' in cnfg:
            label_list = []
            for l in cnfg['labels'] :
                label_list.append(Label(cnfg['labels'][l]))

            self.labels = label_list
        else:
            self.labels = []

            
        if 'markers' in cnfg:
            self.markers = Markers(cnfg['markers'], parent)
        elif hasattr(parent, 'markers'):
            self.markers = Markers(None, parent)

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

        #gc[idx].add_beam()

        for lb in self.labels :
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

        elif sef.vtype == 'box' :
            gp[idx].recenter(self.center[0], self.center[1],
                             height=self.box[0], width=self.box[1])



            
class Dataset(object) :
    """ create a dataset object
    """
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
        #print(center,"-----")
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
    """ Create contours."""

    def __init__(self, cnfg, parent) :

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


            
    def add_label(self, pp, idx) :
        pp[idx].add_label(self.position[0], self.position[1],
                          self.text,
                          relative=self.relative,
                          color=self.color,
                          size=self.size,
                          style=self.style)

        return pp



                                   
class Markers(object) :

    def __init__ (self, cnfg, parent):

        self.wkdir = parent.wkdir
        
        if hasattr(parent, 'markers'):
            self.marklist = parent.markers.marklist
            list_markers = self.marklist
        else:
            list_markers = []
                                                                      
        if cnfg != None :
            file_str = [k for k in cnfg if 'file' in k]

            if file_str :

                for files in file_str :
                    fname = os.path.join(self.wkdir, cnfg[files])
                    list_markers.append(self.read_markers(fname))

                self.marklist = list_markers
                    
            else :
                self.marklist = None
                
        #print("------>", self.marklist)



    def add_markers(self, gc, i):

        lmarkers = self.marklist

        for mk in lmarkers:
            
            gc[i].show_markers(mk[0]['x'].degree, mk[0]['y'].degree,
                               edgecolor=mk[0]['color'],
                               c=mk[0]['filled'],
                               linewidths=mk[0]['linewidth'],
                               s=mk[0]['size'],
                               marker='$+$')



            
    def read_markers(self, fname):
        """ reads a markers file
        """

        tbl = ascii.read(fname, delimiter=" ", format="basic")

        mark_list = []
        for marker in tbl:

            if marker[0] == "cross" :
                props = Markers.check_cross(marker)
                if props != None :
                    mark_list.append(props)

        return mark_list

                    
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
        
    if 'page' in cfg:
        G = Page(cfg['page'], dirs)
    else :
        G = Page()

    if 'panels' in cfg:
        P = DbxFig(cfg['panels'], dirs)

    return G, P



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
    Page, Panels = process_config(files, dirs)


    # make the figure
    #

    fig = Page.create()

    Panels.add_panels(fig)

    Page.f_print()
    Page.end()
    
##-- End of main -------------------------------------------------------
