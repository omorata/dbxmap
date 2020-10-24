#!/usr/bin/env python3
##
## display.py
##
## O. Morata 2020
##
## module containing classes that process and display the map data
##

import copy
import os
import sys

from astropy import wcs
from astropy.io import fits
import numpy as np


import astropy.units as u
from astropy.coordinates import SkyCoord


import markers as mrk
import framework as fw
import copy

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

                
        if 'slices' in cnfg:
            self.slices = cnfg['slices']
        else:
            self.slices = None
            
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
        


    def to_offsets(self, reference, frame='icrs'):
        """ transforms the coordinates of the fits file to a new frame"""

        hdu = fits.open(self.filename)[0]
        w = wcs.WCS(hdu.header)

        center = (w.wcs.crval[0], w.wcs.crval[1])

        ref = SkyCoord(ra=reference[0]*u.arcsec, dec=reference[1]*u.arcsec,
                      frame=frame)

        target = SkyCoord(ra=center[0]*u.arcsec, dec=center[1]*u.arcsec,
                          frame=frame) 

        new = target.transform_to(ref.skyoffset_frame())
        ra, dec = new.lon.arcsec, new.lat.arcsec
                    
        w.wcs.crval = ra, dec
        w.wcs.ctype = 'RA -- XOFFSET', 'DEC -- YOFFSET'
        w.wcs.cunit = 'arcsec', 'arcsec'
        w.wcs.cdelt = wcs.utils.proj_plane_pixel_scales(w)*3600

        hdu_linear = hdu.copy()
        hdu_linear.header = w.to_header()

        return hdu_linear, ra, dec

    
    
    def read_beam_parameters(self, battr):
        """Reads the beam parameters from the configuration file

        It populates self.beam_shape and self.beam_args

        TODO:
            proper unit treatment
        """

        if 'bmaj' in battr :
            bmaj = fw.View.read_units(battr['bmaj'])
            self.beam_shape = {'major' : bmaj}
            
        if 'major' in self.beam_shape :

            if 'bmin' in battr :
                self.beam_shape['minor'] = fw.View.read_units(battr['bmin'])
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
                

            
    def show(self, g, **kwargs) :
        """ show the dataset
            It takes into account whether it is a pixel map or a contour
            map
        """
        if self.dtype == 'pixel' :
            self.show_colorscale(g)

        elif self.dtype == 'cntr' :
            self.show_contour(g, **kwargs)


            
    def show_colorscale(self, g) :
        """Show a colorscale of the dataset."""
            
        pixrg_args = {}
        if self.pixrange.smooth:
            pixrg_args['smooth'] = self.pixrange.smooth

        if self.pixrange.kernel:
            pixrg_args['kernel'] = self.pixrange.kernel

            
        g.show_colorscale(vmin=self.pixrange.range[0],
                          vmax=self.pixrange.range[1],
                          stretch=self.pixrange.stretch,
                          cmap=self.pixrange.colormap,
                          aspect='equal', **pixrg_args)


        
    def show_contour(self, g, coords='abs', ref='') :
        """Show a contour map of the dataset."""
        
        if coords == 'off' :
            hdu, ra, dec = self.to_offsets(ref)
        else :
            hdu = self.filename

        if self.slices :
            hdu = self.get_slice(hdu=hdu)

        cntr_args = {}
        if self.contour.smooth:
            cntr_args['smooth'] = self.contour.smooth

        if self.contour.kernel:
            cntr_args['kernel'] = self.contour.kernel

        
        g.show_contour(data=hdu,levels=self.contour.levels,
                       colors=self.contour.colors,
                       linewidths=self.contour.linewidth,
                       linestyles=self.contour.linestyle,
                       **cntr_args)




    def get_slice(self, hdu):
        """gets a slice out of a data cube"""

        hh = fits.open(hdu)[0]
        w = wcs.WCS(hh.header)

        save_type = []
        naxis = w.naxis
        ssl = copy.deepcopy(self.slices)
                
        for dx in range(naxis):
            rev_idx = naxis - dx - 1
            if dx in self.dims:
                save_type.append(w.wcs.ctype[dx])

            else:
                hh.data = hh.data.take(indices=ssl.pop(0), axis=rev_idx)

        ct = 0
        for dx in range(naxis):
            if not w.wcs.ctype[ct] in save_type:
                w = w.dropaxis(ct)
            else:
                ct += 1

        hh.header = w.to_header()
        
        return hh


        
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


        if cnfg != None and 'smooth' in cnfg:
            self.smooth = cnfg['smooth']
        else:
            try:
                self.smooth = parent.pixrange.smooth

            except AttributeError:
                self.smooth = None


        if cnfg != None and 'kernel' in cnfg:
            self.kernel = cnfg['kernel']
        else:
            try:
                self.kernel = parent.pixrange.kernel
                
            except AttributeError:
                self.kernel = None

            

            
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


        if cnfg != None and 'smooth' in cnfg:
            self.smooth = cnfg['smooth']
        else:
            try:
                self.smooth = parent.contour.smooth

            except AttributeError:
                self.smooth = None


        if cnfg != None and 'kernel' in cnfg:
            self.kernel = cnfg['kernel']
        else:
            try:
                self.kernel = parent.contour.kernel
                
            except AttributeError:
                self.kernel = None

            
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




class Colorbar(object):
    """creat colorbar"""

    property_list = ['location', 'width', 'label_pad', 'label_rotation',
                     'label_text', 'box', 'frame_color', 'frame_linewidth',
                     'labels', 'pad', 'ticks', 'hide']

    # not implemented yet:
    # box
    # frame_color
    # frame_linewidth
    # label_properties ()
    # hide
    # show
    # log_format

    def __init__(self, cnfg, parent, fonts=None) :
        """Initialization of a colorbar object."""

        self.add = 'y'

        for pr in self.property_list:
            if cnfg != None and pr in cnfg:
                setattr(self, pr, cnfg[pr])
            else :
                try:
                    setattr(self, pr, getattr(parent.colorbar, pr))
                except AttributeError:
                    pass


        for pr in ['font', 'label_font']:
            if cnfg != None and pr in cnfg:

                if hasattr(parent, pr):
                    setattr(self, pr,
                            copy.deepcopy(getattr(parent.colorbar, pr)))
                elif fonts != None:
                    setattr(self, pr, self.copyfonts(fonts))
                else :
                    setattr(self, pr, {})

                ff = cnfg[pr]

                for prop in ['family', 'style', 'size', 'variant', 'stretch',
                             'weight']:
                    dct = getattr(self, pr)
                    if prop in ff :
                        dct[prop] = ff[prop]

                    setattr(self, pr, dct)
            else:
                try:
                    setattr(self, pr,
                            copy.deepcopy(getattr(parent.colorbar, pr)))
                except AttributeError:
                    if fonts != None:
                        setattr(self, pr, self.copyfonts(fonts))



    @staticmethod
    def copyfonts(fnt):
        dct = {}
        for p in fnt:
            dct[p] = fnt[p]
        return dct



    def set_colorbar(self, g):
        """Set the colobar in the panel."""

        g.add_colorbar()
        if hasattr(self, 'width'):
            g.colorbar.set_width(float(self.width))
        if hasattr(self, 'location'):
            g.colorbar.set_location(self.location)
        if hasattr(self, 'label_text'):
            g.colorbar.set_axis_label_text(self.label_text)
        if hasattr(self, 'font'):
            g.colorbar.set_font(**self.font)

        if hasattr(self, 'label_rotation'):
            g.colorbar.set_axis_label_rotation(self.label_rotation)
        if hasattr(self, 'label_pad'):
            g.colorbar.set_axis_label_pad(self.label_pad)
        if hasattr(self, 'label_font'):
            g.colorbar.set_axis_label_font(**self.label_font)

        if hasattr(self, 'labels'):
            g.colorbar.set_labels(self.labels)

        if hasattr(self, 'pad'):
            g.colorbar.set_pad(self.pad)

        if hasattr(self, 'ticks'):
            g.colorbar.set_ticks(self.ticks)

        if hasattr(self, 'hide'):
            if self.hide: 
                g.colorbar.hide()

