# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like

def getListOfAlgorithms():
    """
    Print list of supported layout algorithms
    """
    
    algList = ['kamada-kawai', 'spring', 'twopi', 'neato', 'dot']
    algList.sort()
    
    return algList

def checkValidity(self):
    
    if not isinstance(self.scale, (int, float)):
        raise Exception('scale paramter only accepts number')
    elif not isinstance(self.fontsize, (int, float)):
        raise Exception('fontsize paramter only accepts number')
    elif not isinstance(self.edgelw, (int, float)):
        raise Exception('edgelw paramter only accepts number')
    elif not is_color_like(self.nodeColor):
        raise Exception('nodeColor paramter does not look like a color')
    elif not is_color_like(self.reactionNodeColor):
        raise Exception('reactionNodeColor paramter does not look like a color')
    elif not is_color_like(self.labelColor):
        raise Exception('labelColor paramter does not look like a color')
    elif type(self.labelReactionIds) is not bool:
        raise Exception('labelReactionIds paramter only accepts boolean')
    elif not is_color_like(self.reactionColor):
        raise Exception('reactionColor paramter does not look like a color')
    elif not is_color_like(self.modifierColor):
        raise Exception('modifierColor paramter does not look like a color')
    elif not is_color_like(self.boundaryColor):
        raise Exception('boundaryColor paramter does not look like a color')
    elif not is_color_like(self.nodeEdgeColor):
        raise Exception('nodeEdgeColor paramter does not look like a color')
    elif not isinstance(self.nodeEdgelw, (int, float)):
        raise Exception('nodeEdgelw paramter only accepts number')
    elif self.edgeType != 'default' or self.edgeType != 'bezier':
        raise Exception('unknown edgetype')
    elif not is_color_like(self.compartmentColor):
        raise Exception('compartmentColor paramter does not look like a color')
    elif not is_color_like(self.compartmentEdgeColor):
        raise Exception('compartmentEdgeColor paramter does not look like a color')
    elif not isinstance(self.compartmentEdgelw, (int, float)):
        raise Exception('compartmentEdgelw paramter only accepts number')
    elif type(self.highlight) is not list:
        raise Exception('highlight paramter only accepts list')
    elif not is_color_like(self.hlNodeColor):
        raise Exception('hlNodeColor paramter does not look like a color')
    elif not is_color_like(self.hlNodeEdgeColor):
        raise Exception('hlNodeEdgeColor paramter does not look like a color')
    elif type(self.drawReactionNode) is not bool:
        raise Exception('drawReactionNode paramter only accepts boolean')
    elif type(self.breakBoundary) is not bool:
        raise Exception('breakBoundary paramter only accepts boolean')
    elif type(self.tightLayout) is not bool:
        raise Exception('tightLayout paramter only accepts boolean')
    elif type(self.analyzeFlux) is not bool:
        raise Exception('analyzeFlux paramter only accepts boolean')
    elif type(self.analyzeRates) is not bool:
        raise Exception('analyzeRates paramter only accepts boolean')
    elif not is_color_like(self.analyzeColorHigh):
        raise Exception('analyzeColorHigh paramter does not look like a color')
    elif not is_color_like(self.analyzeColorLow):
        raise Exception('analyzeColorLow paramter does not look like a color')
    elif self.analyzeColorMap not in plt.colormaps():
        raise Exception('analyzeColorMap paramter does not look like a colormap')
    elif type(self.analyzeColorScale) is not bool:
        raise Exception('analyzeColorScale paramter only accepts boolean')
    elif type(self.drawInlineTimeCourse) is not bool:
        raise Exception('drawInlineTimeCourse paramter only accepts boolean')
    elif not isinstance(self.nodeEdgelw, (int, float)):
        raise Exception('nodeEdgelw paramter only accepts number')
    elif not isinstance(self.nodeEdgelw, (int, float)):
        raise Exception('nodeEdgelw paramter only accepts number')
    elif not isinstance(self.simulationStartTime, (int, float)):
        raise Exception('simulationStartTime paramter only accepts number')
    elif not isinstance(self.simulationEndTime, (int, float)):
        raise Exception('simulationEndTime paramter only accepts number')
    elif not isinstance(self.numPoints, (int)):
        raise Exception('numPoints paramter only accepts integer')
    elif type(self.plotStatistics) is not bool:
        raise Exception('plotStatistics paramter only accepts boolean')
    elif type(self.forceAnalysisAtEndTime) is not bool:
        raise Exception('forceAnalysisAtEndTime paramter only accepts boolean')
    elif type(self.plotColorbar) is not bool:
        raise Exception('plotColorbar paramter only accepts boolean')
    elif type(self.inlineTimeCourseSelections) is not list:
        raise Exception('inlineTimeCourseSelections paramter only accepts list')
    elif self.layoutAlgorithm != 'default' or self.edgeType != 'bezier':
        raise Exception('unknown edgetype')
    elif type(self.ignoreLayout) is not bool:
        raise Exception('ignoreLayout paramter only accepts boolean')
    
