=============
netplotlib API
=============

--------------
Single Network
--------------

.. method:: netplotlib.Network(model)

    :param model: SBML or Antimony string of a model
        
    An instance of Network object accepts following properties:
    
        - scale: scaling factor for layout algorithm
        - fontsize: fontsize for labels
        - edgelw: linewidth of edges
        - nodeColor: node color
        - reactionNodeColor: reaction node color
        - labelColor: label color
        - labelReactionIds: boolean flag for labeling reaction ids
        - reactionColor: edge color
        - modifierColor: modifier edge color
        - boundaryColor: boundary node color
        - nodeEdgeColor: node edge color
        - nodeEdgelw: linewidth of node edges
        - highlight: list of species ids or reaction ids to highlight
        - hlNodeColor: node color of highlighted nodes
        - hlNodeEdgeColor: node edge color of highlighted nodes
        - drawReactionNode: boolean flag for drawing reaction nodes
        - breakBoundary: boolean flag for breaking all boundary species into separate nodes
        - analyzeFlux: boolean flag for visualizing flux
        - analyzeRates: boolean flag for visualizing species rate of changes
        - analyzeColorHigh: color to use for higher values during analysis
        - analyzeColorLow: color to use for lower values during analysis
        - analyzeColorMap: colormap to use for analysis. 
        - analyzeColorScale: boolean flag for using colormaps. Setting this true ignore analyzeColorHigh and analyzeColorLow
        - drawInlineTimeCourse: boolean flag for plotting inline time-cource simulation
        - simTime: value for simulation duration
        - forceAnalysisAtSimTime: boolean flag for running analysis at the end of simTime
        - plotColorbar: boolean flag for visualizing color bar
        - inlineTimeCourseSelections: list of species to plot for time-course simulation
    
    .. method:: Network.draw(show=True, savePath=None)
    
        :param show: boolean flag to show the diagram
        :param savePath: path to save the diagram
    
        Draws network diagram
        
    .. method:: Network.getLayout(returnState=False)
    
        :param returnState: boolean flag for returning the networkx.Graph object
        
        Returns the layout
        
    .. method:: Network.reset()
    
        Resets all propertiess
    
    .. method:: Network.savefig(path)
    
        :param path: path to save the diagram
        
        Save a network diagram to specified location
    

----------------
Network Ensemble
----------------

.. method:: netplotlib.NetworkEnsemble(models)

    :param models: list of SBML or Antimony strings of models
        
    An instance of NetworkEnsemble object accepts following properties:
    
        - scale: scaling factor for layout algorithm
        - fontsize: fontsize for labels
        - edgelw: linewidth of edges
        - nodeColor: node color
        - reactionNodeColor: reaction node color
        - labelColor: label color
        - labelReactionIds: boolean flag for labeling reaction ids
        - reactionColor: edge color
        - modifierColor: modifier edge color
        - boundaryColor: boundary node color
        - nodeEdgeColor: node edge color
        - nodeEdgelw: linewidth of node edges
        - highlight: list of species ids or reaction ids to highlight
        - hlNodeColor: node color of highlighted nodes
        - hlNodeEdgeColor: node edge color of highlighted nodes
        - edgeLabel: boolean flag for displaying edge weights
        - edgeLabelFontSize: fontsize of edge weight labels
        - drawReactionNode: flag for drawing reaction nodes
        - breakBoundary: boolean flag for breaking all boundary species into separate nodes
        - weights: list of custom weights to override
        - edgeTransparency: boolean flag for changing the transparency of the edges accordin to edge weights
        - plottingThreshold: value of threshold to prevent from displaying weighted edges
        - removeBelowThreshold: boolean flag for preventing weighted edges below plottingThreshold from displaying
        - analyzeFlux: boolean flag for visualizing flux
    
    .. method:: Network.drawWeightedDiagram(show=True, savePath=None)
    
        :param show: boolean flag to show the diagram
        :param savePath: path to save the diagram
        
        Draw weighted reaction network based on frequency of reactions
       
    .. method:: Network.drawNetworkGrid(nrows, ncols, auto=False, show=True, savePath=None)
    
        :param nrows: number of rows
        :param ncols: number of columns
        :param auto: boolean flag to automatically figure out optimal number of rows and columns. Setting this to True will ignore nrows and ncols.
        :param show: boolean flag to show the diagram
        :param savePath: path to save the diagram
        
        Plot a grid of network diagrams
        
    .. method:: Network.getLayout()
    
        Return the layout
        
    .. method:: Network.reset()
    
        Resets all propertiess
        
    .. method:: Network.savefig(path)
    
        :param path: path to save the diagram
        
        Save a weighted network diagram to specified location


-------------
Miscellaneous
-------------

.. autofunction:: netplotlib.netplotlib.getVersion

