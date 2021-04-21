# -*- coding: utf-8 -*-
"""
Layout extension support for netplorlib library

Kiri Choi (c) 2021 
"""

def readLayout(layoutPlugin):
    pos = dict()
    
    layout = layoutPlugin.getLayout(0) # TODO: currently support only the first layout
    layoutdim = layout.getDimensions()
    for i in range(len(layout.getNumSpeciesGlyphs())):
        sg = layout.getSpeciesGlyph(i)
        bbox = sg.getBoundingBox()
        sgpos = bbox.getPosition()
        pos[] = np.array([sgpos.x_offset(), sgpos.y_offset()])
        
        dim = bbox.getDimensions()
        dim.getWidth()