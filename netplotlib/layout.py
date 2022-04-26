# -*- coding: utf-8 -*-
"""
Layout extension support for netplorlib library

Kiri Choi (c) 2021 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyArrowPatch, FancyBboxPatch, ArrowStyle, 
                                Arc, RegularPolygon, Rectangle, PathPatch,
                                Circle)
from matplotlib import gridspec, cm, colors
from matplotlib.path import Path
from matplotlib.transforms import Bbox
from matplotlib.text import Text
from matplotlib.patheffects import AbstractPathEffect
import libsbml

def getPos(layoutPlugin):
    """
    Return the layout of the model
    
    :param layoutPlugin: SBML layout plugin
    :returns pos: Dictionary of all nodes and corresponding coordinates
    """
    
    pos = dict()
    
    layout = layoutPlugin.getLayout(0) # TODO: currently support only the first layout
    
    for i in range(layout.getNumCompartmentGlyphs()):
        cg = layout.getCompartmentGlyph(i)
        bbox = cg.getBoundingBox()
        cgpos = bbox.getPosition()
        pos[cg.getCompartmentId()] = np.array([cgpos.x_offset, cgpos.y_offset])
    
    for i in range(layout.getNumSpeciesGlyphs()):
        sg = layout.getSpeciesGlyph(i)
        bbox = sg.getBoundingBox()
        sgpos = bbox.getPosition()
        pos[sg.getSpeciesId()] = np.array([sgpos.x_offset, sgpos.y_offset])
    
    return pos

def stretch_font(text, width, height):
    
    class TextScaler(AbstractPathEffect):
        
        def __init__(self, text, width, height):
            self._text = text
            self._width = width
            self._height = height

        def draw_path(self, renderer, gc, tpath, affine, rgbFace=None):
            ax = self._text.axes
            renderer = ax.get_figure().canvas.get_renderer()
            bbox = text.get_window_extent(renderer=renderer)
            bbox = Bbox(ax.transData.inverted().transform(bbox))

            scale_x = self._width/bbox.width
            scale_y = self._height/bbox.height

            affine = affine.identity().scale(scale_x, scale_y) + affine
            renderer.draw_path(gc, tpath, affine, rgbFace)

    text.set_path_effects([TextScaler(text, width, height)])

def draw(NetworkClass, show=True, savePath=None, dpi=150):
    """
    Draw network diagram using SBML Layout package
    
    :param show: flag to show the diagram
    :param savePath: path to save the diagram
    :param dpi: dpi settings for the diagram
    """
    
    size = dict()
    
    fig = plt.figure()

    if NetworkClass.customAxis != None:
        ax = NetworkClass.customAxis
    else:
        ax = plt.gca()
    
    layout = NetworkClass._Var.layoutPlugin.getLayout(0) # TODO: currently support only the first layout
    
    layoutCompartmentGlyphIds = []
    layoutSpeciesGlyphIds = []
    layoutReactionGlyphIds = []
    layoutTextGlyphIds = []

    for n in range(layout.getNumCompartmentGlyphs()):
        cg = layout.getCompartmentGlyph(n)
        bbox = cg.getBoundingBox()
        cgpos = bbox.getPosition()
        layoutCompartmentGlyphIds.append(cg.getCompartmentId())
        
        dim = bbox.getDimensions()
        size[cg.getCompartmentId()] = np.array([dim.getWidth(), dim.getHeight()])
    
        node_color = NetworkClass.compartmentColor
    
        c = FancyBboxPatch((cgpos.x_offset,
                            cgpos.y_offset),
                            dim.getWidth(), 
                            dim.getHeight(),
                            boxstyle="round,pad=0.01, rounding_size=0.02",
                            linewidth=NetworkClass.compartmentEdgelw, 
                            edgecolor=NetworkClass.compartmentEdgeColor, 
                            facecolor=node_color,
                            alpha=0.5)
        
        ax.add_patch(c)
        
    for n in range(layout.getNumSpeciesGlyphs()):
        sg = layout.getSpeciesGlyph(n)
        bbox = sg.getBoundingBox()
        sgpos = bbox.getPosition()
        layoutSpeciesGlyphIds.append(sg.getSpeciesId())
        
        dim = bbox.getDimensions()
        size[sg.getSpeciesId()] = np.array([dim.getWidth(), dim.getHeight()])
        
        if (sg.getSpeciesId() in NetworkClass._Var.boundaryId):
            node_color = NetworkClass.boundaryColor
        else:
            node_color = NetworkClass.nodeColor
        
        c = FancyBboxPatch((sgpos.x_offset,
                            sgpos.y_offset),
                            dim.getWidth(), 
                            dim.getHeight(),
                            boxstyle="round,pad=0.01, rounding_size=0.02",
                            linewidth=NetworkClass.nodeEdgelw, 
                            edgecolor=NetworkClass.nodeEdgeColor, 
                            facecolor=node_color)
        
        NetworkClass._Var.G.nodes[sg.getSpeciesId()]['patch'] = c
    
        ax.add_patch(c)
        
    for n in range(layout.getNumReactionGlyphs()):
        rg = layout.getReactionGlyph(n)
        
        cg = rg.getCurve()
        
        for c in range(cg.getNumCurveSegments()):
            csg = cg.getCurveSegment(c)
            if type(csg) == libsbml.LineSegment:
                csgs = csg.getStart()
                csge = csg.getEnd()
                
                e = PathPatch(Path([(csgs.getXOffset(), csgs.getYOffset()), 
                                    (csge.getXOffset(), csge.getYOffset())], 
                                   [Path.MOVETO, Path.LINETO]), 
                              color=NetworkClass.reactionColor,
                              lw=(1+NetworkClass.edgelw))
            else:
                csgs = csg.getStart()
                csge = csg.getEnd()
                csgcp1 = csg.getBasePoint1()
                csgcp2 = csg.getBasePoint2()
                
                e = FancyArrowPatch(path=Path([(csgs.getXOffset(), csgs.getYOffset()),
                                        (csgcp1.getXOffset(), csgcp1.getYOffset()),
                                        (csgcp2.getXOffset(), csgcp2.getYOffset()),
                                        (csge.getXOffset(), csge.getYOffset())], 
                                       [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]), 
                                  color=NetworkClass.reactionColor,
                                  lw=(1+NetworkClass.edgelw))
            ax.add_patch(e)
        
        for b in range(rg.getNumSpeciesReferenceGlyphs()):
            srg = rg.getSpeciesReferenceGlyph(b)
            srgc = srg.getCurve()
            
            if srg.getRoleString() == 'substrate':
                arrowstyle = ArrowStyle.Curve()
            elif srg.getRoleString() == 'sidesubstrate':
                arrowstyle = ArrowStyle.Curve()
            elif srg.getRoleString() == 'product':
                arrowstyle = ArrowStyle.CurveFilledB(head_length=2.5*NetworkClass.edgelw, head_width=2*NetworkClass.edgelw)
            elif srg.getRoleString() == 'sideproduct':
                arrowstyle = ArrowStyle.CurveFilledB(head_length=2.5*NetworkClass.edgelw, head_width=2*NetworkClass.edgelw)
            elif srg.getRoleString() == 'activator':
                arrowstyle = ArrowStyle.Curve()
            elif srg.getRoleString() == 'inhibitor':
                arrowstyle = ArrowStyle.BarAB(widthA=0.0, angleA=None, widthB=3*NetworkClass.edgelw, angleB=None)
            else:
                arrowstyle = ArrowStyle.Curve()
                
            for bc in range(srgc.getNumCurveSegments()):
                srgcs = srgc.getCurveSegment(bc)
                if type(srgcs) == libsbml.LineSegment:
                    srgcss = srgcs.getStart()
                    srgcse = srgcs.getEnd()
                    
                    e = FancyArrowPatch(path=Path([(srgcss.getXOffset(), srgcss.getYOffset()), 
                                        (srgcse.getXOffset(), srgcse.getYOffset())], 
                                       [Path.MOVETO, Path.LINETO]),
                                  arrowstyle=arrowstyle,
                                  color=NetworkClass.reactionColor,
                                  lw=(1+NetworkClass.edgelw))
                
                else:
                    srgcss = srgcs.getStart()
                    srgcse = srgcs.getEnd()
                    srgcscp1 = srgcs.getBasePoint1()
                    srgcscp2 = srgcs.getBasePoint2()
                    
                    e = FancyArrowPatch(path=Path([(srgcss.getXOffset(), srgcss.getYOffset()),
                                            (srgcscp1.getXOffset(), srgcscp1.getYOffset()),
                                            (srgcscp2.getXOffset(), srgcscp2.getYOffset()),
                                            (srgcse.getXOffset(), srgcse.getYOffset())], 
                                           [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]), 
                                        arrowstyle=arrowstyle,
                                        color=NetworkClass.reactionColor,
                                        lw=(1+NetworkClass.edgelw))
                    
                ax.add_patch(e)
                
            if srg.getRoleString() == 'activator':
                e = Circle((srgcse.getXOffset(), srgcse.getYOffset()),
                           radius=3*NetworkClass.edgelw,
                           color=NetworkClass.reactionColor,
                           lw=(1+NetworkClass.edgelw))
                ax.add_patch(e)
    
    for n in range(layout.getNumTextGlyphs()):
        tg = layout.getTextGlyph(n)
        bbox = tg.getBoundingBox()
        tgpos = bbox.getPosition()
        speciesName = NetworkClass._Var.sbmlmodel.getSpecies(tg.getOriginOfTextId())
        if speciesName.getName() != '':
            speciesName = speciesName.getName()
        else:
            speciesName = speciesName.getId()
        layoutTextGlyphIds.append(speciesName)
        
        dim = bbox.getDimensions()
        
        mattext = ax.text(tgpos.x_offset, tgpos.y_offset, speciesName,
                          horizontalalignment='left', verticalalignment='bottom',
                          color=NetworkClass.labelColor)
        
        stretch_font(mattext, dim.getWidth(), dim.getHeight())
    
    layoutDim = layout.getDimensions()
    plt.xlim(0, layoutDim.getWidth())
    plt.ylim(0, layoutDim.getHeight())
    plt.axis('off')
    plt.show()
    
    if savePath != None:
        try:
            fig.savefig(savePath, bbox_inches='tight', dpi=dpi)
        except IOError as e:
            raise Exception("Error saving diagram: " + str(e))
            
    if show and NetworkClass.customAxis == None:
        plt.show()
    plt.close()
