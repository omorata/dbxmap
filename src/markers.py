#!/usr/bin/env python3
##
## markers.py
##
## O. Morata 2020
##
## module containing classes related to makers and labels
##

import copy
import os
import sys

import astropy.coordinates as coord
from astropy.io import ascii
import astropy.units as u
import numpy as np

import framework as fw


class Marker(object) :
    """Class to define markers (including polygons)."""

    style_props = ['edgecolor', 'linewidth', 'linestyle', 'facecolor', 'alpha',
                   'zorder'] 
    lstyle_props = ['color', 'family', 'style', 'variant', 'stretch', 'weight',
                    'size']

    props = ['type', 'show_label', 'lpad', 'alpha', 'zorder', 'linecolor']
    
    ini_properties = ['show_label', 'lpad', 'zorder']
    defaults = [False, [0,0], 5]

    
    def __init__ (self, cfg, parent, fonts=None):

        property_list = self.props[:]
        property_list.extend(self.style_props)
        property_list.extend(self.lstyle_props)

        self.wkdir = parent.wkdir
        
        if hasattr(parent, 'markers'):
            self.marklist = copy.deepcopy(parent.markers.marklist)
        else:
            self.marklist = []

        if hasattr(parent, 'markers'):
            self.marker_props = copy.deepcopy(parent.markers.marker_props)
        elif hasattr(parent, 'marker_props'):
            self.marker_props = copy.deepcopy(parent.marker_props)
        else :
            self.marker_props = self.default_marker_props()

            if fonts != None:
                for p in fonts:
                    self.marker_props[p] = fonts[p]

        if cfg != None :
            for prop in property_list :
                if prop in cfg:
                    self.marker_props[prop] = cfg[prop]

            if 'file' in cfg:
                fname = os.path.join(self.wkdir, cfg['file'])

                self.marker = self.read_markers(fname)

            marker_str = [k for k in cfg if 'marker' in k]

            if marker_str :
                    
                for mrk in marker_str :
                    self.marklist.append(Marker(cfg[mrk], self))


                
    def default_marker_props(self):
        """Define default values for marker properties."""

        props = { key:value for key, value in
                  zip(self.ini_properties, self.defaults) }
        return props


    
    def add_markers(self, g):
        """Adds markers to the panel."""

        marker_fields = self.marker
        
        for mk in marker_fields:

            if mk['type'] == 'polygon' :
                g.show_polygons(mk['corners'], **mk['style'])

            elif mk['type'] == "ellipse" :
                g.show_ellipses(mk['x'], mk['y'],
                                mk['maj'], mk['min'], angle=mk['pa'],
                                **mk['style'])

            elif mk['type'] == "line" :
                g.show_lines(mk['line'], color=mk['linecolor'],
                             **mk['style'])

            else :
                g.show_markers(mk['x'], mk['y'], marker=mk['sym'], s=mk['size'],
                               c=mk['c'], **mk['style'])

            if self.marker_props['show_label'] :
                x = mk['x'] + mk['lpad'][0] / 3600.
                y = mk['y'] + mk['lpad'][1] / 3600.

                g.add_label(x, y, mk['id'], **mk['l_style'], clip_on=True)



    def read_markers(self, fname):
        """Reads a markers file."""

        tbl = ascii.read(fname, delimiter=" ", format="basic")

        mark_list = []
        properties = None
        for marker in tbl:

            if self.marker_props['type'] == 'polygon' :
                if marker['type'] == 'Polygon' :
                    properties = self.read_polygon(marker)
                else :
                    print("\n  +++ ERROR: marker type polygon reading a",
                          "file without Polygons\n")
                    sys.exit(1)

            elif marker['type'] == 'ellipse' :
                properties = self.read_ellipse(marker)
            elif marker['type'] == 'line' :
                properties = self.read_line(marker)
            else:
                properties = self.read_symbol(marker)

            if properties != None:
                properties = self.read_styles(properties)
                mark_list.append(properties)
            else :
                print("ERROR: unknown marker")
                sys.exit(1)

        return mark_list



    def read_styles(self, attrib):
        """Adds the line and text styles to the marker."""

        style = {}
        for ss in self.style_props:
            if ss in self.marker_props:
                style[ss] =  self.marker_props[ss]

        attrib['style'] = style

        labelstyle = {}
        for ls in self.lstyle_props:
            if ls in self.marker_props :
                labelstyle[ls] =  self.marker_props[ls]

        attrib['l_style'] = labelstyle

        attrib['lpad'] = self.marker_props['lpad']

        return attrib


    
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
        attrib['x'] = np.average(array_corners[:,0])
        attrib['y'] = np.average(array_corners[:,1])

        return attrib



    def read_symbol(self, it):
        """Read the definition of a symbol marker.

        List of valid symbols in matplotlib:
        https://matplotlib.org/api/markers_api.html#module-matplotlib.markers
        """

        attrib = {'id' : it['id'], 'type' : self.marker_props['type'] }

        if '(' and ')' in self.marker_props['type']:
            attrib['sym'] = self.read_tuple('type')
        else :
            attrib['sym'] = self.marker_props['type']

        center = it['center'].split(" ")
        if it['coords'] == 'world_deg' :
            attrib['x'] = float(center[0])
            attrib['y'] = float(center[1])
        elif it['coords'] == 'radec' :
            attrib['x'], attrib['y'] = self.readangles(center, 'radec')

        attrib['size'] = float(it['size'])

        if self.marker_props['facecolor'] == 'None' :
            attrib['c'] = None
        else :
            attrib['c'] = self.marker_props['facecolor']

        return attrib



    def read_ellipse(self, it):
        """Read the definition of a ellipse.

        It will transform, if necessary, units to degrees.
        """

        attrib = {'id' : it['id'], 'type' : it['type'] }

        center = it['center'].split(" ")
        if it['coords'] == 'world_deg' :
            attrib['x'] = float(center[0])
            attrib['y'] = float(center[1])
        elif it['coords'] == 'radec' :
            attrib['x'], attrib['y'] = self.readangles(center, 'radec')

        size = it['size'].split(" ")
        attrib['maj'] = fw.View.read_units(size[0])
        attrib['min'] = fw.View.read_units(size[1])
        attrib['pa'] = fw.View.read_units(size[2])

        return attrib



    def read_line(self, it):
        """Read the definition of a line.

        It will transform, if necessary, units to degrees.
        """

        attrib = {'id' : it['id'], 'type' : it['type'] }

        center = it['center'].split(" ")
        size = it['size'].split(" ")

        if it['coords'] == 'world_deg' :
            attrib['x'] = float(center[0])
            attrib['y'] = float(center[1])
        elif it['coords'] == 'radec' :
            attrib['x'], attrib['y'] = self.readangles(center, 'radec')

        array = np.array([ [attrib['x'], float(size[0])],
                          [attrib['y'], float(size[1])] ])

        attrib['line'] = [array]

        attrib['linecolor'] = self.marker_props['linecolor']

        return attrib



    def read_tuple(self, key):
        """Read a tuple from the file."""

        tuples = []

        tval = self.marker_props[key].strip('()').split(',')
        
        if len(tval) == 2 :
            tval.append(0)
        elif len(tval) == 1 :
            tval.extend([0,0])

        for i in range(3):
            if tval[i] == '':
                tval[i] = 0

            tuples.append(int(tval[i]))

        return tuples



    @staticmethod
    def readangles(c, angle_units):
        """Convert input angle units to degrees.

        Right now, only "radec" is implemented.
        To be implemented: "offset" "lb"
        """

        if angle_units == 'radec':
            n = np.shape(c)[0]
            if n == 2:
                if ':' in c[0]:
                    ra = (coord.Angle(c[0]+' hour')).to(u.deg).degree
                else :
                    ra = (coord.Angle(c[0])).to(u.deg).degree

                if ':' in c[1]:
                    dec = (coord.Angle(c[1]+' degrees')).degree
                else :
                    dec = (coord.Angle(c[1])).degree

            elif n == 6:
                jra = " ".join(c[:3])
                jdec = " ".join(c[3:])
                ra = (coord.Angle(jra+' hour')).to(u.deg).degree
                dec = (coord.Angle(jdec+' degrees')).degree
            else:
                print(" >>> ERROR:Wrong format for the input coordinates")
                sys.exit(1)

            return ra, dec




class Label(object):
    """Create a label."""

    property_list = ['text', 'relative', 'position', 'color', 'size',
                     'style', 'family', 'variant', 'stretch', 'weight']
    
    ini_properties = [ 'relative']
    defaults = [False]

    def __init__(self, cfg, parent, fonts=None):

        if hasattr(parent, 'labels') :
            self.label_list = copy.deepcopy(parent.labels.label_list)
        else:
            self.label_list = []

        if hasattr(parent, 'labels') :
            self.label_props = copy.deepcopy(parent.labels.label_props)
        elif hasattr(parent, 'label_props'):
            self.label_props = copy.deepcopy(parent.label_props)
        else:
            self.label_props = self.default_label_props()

            if fonts != None:
                for p in fonts:
                    self.label_props[p] = fonts[p]



        if cfg != None:
            for prop in self.property_list :
                if prop in cfg:
                    self.label_props[prop] = cfg[prop]

            file_str = [k for k in cfg if 'label' in k]

            if file_str :
                for lab in file_str:
                    self.label_list.append(Label(cfg[lab], self))
        #else :
        #    print("\n  +++ WARNING: naked 'label' field in the configuration",
        #          "file +++\n")



    def add_label(self, p) :
        """Adds label to the panel."""

        p.add_label(self.label_props['position'][0],
                    self.label_props['position'][1], **self.label_props)



    def default_label_props(self):
        """Define default values for label properties."""

        props = { key:value for key, value in
                  zip(self.ini_properties, self.defaults) }

        return props


