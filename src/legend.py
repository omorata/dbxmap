#!/usr/bin/env python3
##
## legend.py
## O. Morata 2021
##
## module containing the class Legend to show legends of the lines and
## markers in the plot
##


import matplotlib.lines as mlines
from pprint import pprint


class Legend(object) :
    """ create legend """

    def __init__(self, cnfg, parent):
    

        if cnfg != None and 'loc' in cnfg:
            self.loc = cnfg['loc']
        else :
            try:
                self.loc = parent.legend.loc

            except AttributeError :
                self.loc = 'lower right'

                
        if cnfg != None and 'fontsize' in cnfg:
            self.fontsize = cnfg['fontsize']
        else :
            try:
                self.fontsize = parent.legend.fontsize

            except AttributeError :
                self.fontsize = 'x-small'

        if cnfg != None and 'marker_size' in cnfg:
            self.markersize = cnfg['marker_size']
        else :
            try:
                self.markersize = parent.legend.markersize

            except AttributeError :
                self.markersize = 50

                
        if cnfg != None and 'frameon' in cnfg:
            self.frameon = cnfg['frameon']
        else :
            try:
                self.frameon = parent.legend.frameon

            except AttributeError :
                self.loc = True

        self.handles = []



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
                                                  markersize=self.markersize,
                                                  linestyle="",
                                                  linewidth=b._linewidths,
                                                  label=l[2:]))

            

    def show_legend(self, g):
        """create the legend"""
        
        if self.handles :
            g.ax.legend(handles=self.handles, loc=self.loc,
                         fontsize=self.fontsize, frameon=self.frameon,
                         fancybox=True, shadow=False)
    
