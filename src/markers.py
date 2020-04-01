#!/usr/bin/env python3
##
## markers.py
##
## O. Morata 2020
##
## module containing classes related to makers and labels
##

import os
import sys
import numpy as np

from astropy.io import ascii


class Marker(object) :
    """Class to define markers (including polygons)."""
    
    properties = ['type', 'edgecolor', 'linewidth', 'linestyle',
                  'facecolor', 'show_label', 'color', 'weight', 'size',
                  'lpad', 'alpha', 'zorder', 'linecolor']
    defaults = ['', 'black', 1.0, 'solid', 'none', False, 'black',
                'normal', 12, [0,0], 1.0, 50, 'black']

    # split them style, l_style, and the rest
    
    def __init__ (self, cfg, parent):

        self.wkdir = parent.wkdir
        
        if hasattr(parent, 'markers'):
            self.marklist = parent.markers.marklist.copy()
        else:
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

                self.marker = self.read_markers(fname)

            marker_str = [k for k in cfg if 'marker' in k]

            if marker_str :
                    
                for mrk in marker_str :
                    self.marklist.append(Marker(cfg[mrk], self))

                

                
    def default_marker_props(self):
        """Define default values for marker properties."""

        props = { key:value for key, value in
                  zip(self.properties, self.defaults) }
        return props


    
    def add_markers(self, gc, i):
        """Adds markers to the panel."""

        marker_fields = self.marker
        
        for mk in marker_fields:

            if mk['type'] == 'polygon' :
                gc[i].show_polygons(mk['corners'], **mk['style'])

            elif mk['type'] == "ellipse" :
                gc[i].show_ellipses(mk['x'], mk['y'],
                                   mk['maj'], mk['min'], angle=mk['pa'],
                                   **mk['style'])

            elif mk['type'] == "line" :
                gc[i].show_lines(mk['line'], color=mk['linecolor'],
                                 **mk['style'])

            else :
                gc[i].show_markers(mk['x'], mk['y'],
                                   marker=mk['sym'],
                                   s=mk['size'],
                                   c=mk['c'],
                                   **mk['style'])

            if self.marker_props['show_label'] :
                x = mk['x'] + mk['lpad'][0] / 3600.
                y = mk['y'] + mk['lpad'][1] / 3600.

                gc[i].add_label(x, y, mk['id'], **mk['l_style'])



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
        for ss in ['edgecolor', 'linewidth', 'linestyle', 'facecolor',
                   'alpha', 'zorder']:
            style[ss] =  self.marker_props[ss]

        attrib['style'] = style

        labelstyle = {}
        for ls in ['color', 'weight', 'size']:
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

        attrib['size'] = float(it['size'])

        if self.marker_props['facecolor'] == 'None' :
            attrib['c'] = None
        else :
            attrib['c'] = self.marker_props['facecolor']

        return attrib



    def read_ellipse(self, it):
        """Read the definition of a ellipse.

        it assumes the major and minor axes are in arcsec.
        """

        attrib = {'id' : it['id'], 'type' : it['type'] }

        center = it['center'].split(" ")
        if it['coords'] == 'world_deg' :
            attrib['x'] = float(center[0])
            attrib['y'] = float(center[1])

        size = it['size'].split(" ")
        if it['coords'] == 'world_deg' :
            attrib['maj'] = float(size[0]) / 3600.
            attrib['min'] = float(size[1]) / 3600.
            attrib['pa'] = float(size[2])

        return attrib



    def read_line(self, it):
        """Read the definition of a line.

        it assumes the major and minor axes are in arcsec.
        """

        attrib = {'id' : it['id'], 'type' : it['type'] }

        center = it['center'].split(" ")
        size = it['size'].split(" ")

        if it['coords'] == 'world_deg' :
            attrib['x'] = float(center[0])
            attrib['y'] = float(center[1])

            array = np.array([
                [attrib['x'], float(size[0])],
                [attrib['y'], float(size[1])]
            ])

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

