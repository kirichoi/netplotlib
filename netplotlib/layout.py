# -*- coding: utf-8 -*-
"""
Layout extension support for netplorlib library

Kiri Choi (c) 2021 
"""

import numpy as np

def getPos(layoutPlugin):
    pos = dict()
    
    layout = layoutPlugin.getLayout(0) # TODO: currently support only the first layout
    layoutdim = layout.getDimensions()
    for i in range(layout.getNumSpeciesGlyphs()):
        sg = layout.getSpeciesGlyph(i)
        bbox = sg.getBoundingBox()
        sgpos = bbox.getPosition()
        pos[sg.getSpeciesId()] = np.array([sgpos.x_offset, sgpos.y_offset])
        
    return pos

def readLayout(layoutPlugin):
    pos = dict()
    size = dict()
    
    layout = layoutPlugin.getLayout(0) # TODO: currently support only the first layout
    layoutdim = layout.getDimensions()
    for i in range(layout.getNumSpeciesGlyphs()):
        sg = layout.getSpeciesGlyph(i)
        bbox = sg.getBoundingBox()
        sgpos = bbox.getPosition()
        pos[sg.getSpeciesId()] = np.array([sgpos.x_offset, sgpos.y_offset])
        
        dim = bbox.getDimensions()
        size[sg.getSpeciesId()] = np.array([dim.getWidth(), dim.getHeight()])
    return pos, size