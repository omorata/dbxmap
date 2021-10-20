#!/usr/bin/env python3
##
## legend.py
## O. Morata 2021
##
## module containing the class Legend to show legends of the lines and
## markers in the plot
##

import copy

import matplotlib.lines as mlines


class Legend(object) :
    """ create legend """

    def __init__(self, cnfg, parent):

        self.legend_props = {}
        
        if hasattr(parent, 'legend'):
            if hasattr(parent.legend, 'legend_props') :
                self.legend_props = copy.deepcopy(parent.legend.legend_props)


        if cnfg != None:
            for icf in cnfg.keys():
                self.legend_props[icf] = cnfg[icf]

        if 'bbox_to_anchor' in self.legend_props:
            self.legend_props['bbox_to_anchor'] =  self.read_tuple('float',
                                                            'bbox_to_anchor')
            print(self.legend_props['bbox_to_anchor'])
        self.handles = []



    def read_tuple(self, valtype, key):
        """Process a tuple from the yaml file."""

        tuples = []

        tval = self.legend_props[key].strip('()').split(',')
        for i, val in enumerate(tval):
            if valtype == 'float':
                tuples.append(float(tval[i]))

        return tuples



    def add_legend(self, g):
        """ add legend to panel"""
        
        self.read_layers(g)
        self.show_legend(g)


        
    def read_layers(self, g):
        """ reads the defined layers in the panel

            Currently, it recognizes two types of layers:
              - contours or lines, with name starting with c_
              - markers, with name starting with m_
        """

        layers = g._layers
        elements = list(layers.keys())
        
        for l in elements:

            if "c_" in l:
                b = g.get_layer(l)
                self.handles.append(mlines.Line2D([],[], color=b.colors,
                                                  linewidth=b.linewidths,
                                                  linestyle=b.linestyles,
                                                  label=l[2:]))

            elif "m_" in l:
                b = g.get_layer(l)

                pppath = b._paths[0]
                self.handles.append(mlines.Line2D([],[], color=b._facecolors,
                                                  marker=pppath,
                                                  linestyle="",
                                                  linewidth=b._linewidths,
                                                  label=l[2:]))

            

    def show_legend(self, g):
        """create the legend"""
        
        if self.handles :
            g.ax.legend(handles=self.handles, **self.legend_props )
    
