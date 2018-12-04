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
        - reactionColor: edge color
        - modifierColor: modifier edge color
        - boundaryColor: boundary node color
        - nodeEdgeColor: node edge color
        - nodeEdgelw: linewidth of node edges
        - highlight: list of species ids or reaction ids to highlight
        - hlNodeColor: node color of highlighted nodes
        - hlNodeEdgeColor: node edge color of highlighted nodes
        - drawReactionNode: flag for drawing reaction nodes
        - breakBoundary: flag for breaking all boundary species into separate nodes
    
    .. method:: Network.draw()
    
        Draw network diagram
        
    .. method:: Network.getLayout()
    
        Return the layout
        
    .. method:: Network.reset()
    
        Resets all propertiess
    

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
        - reactionColor: edge color
        - modifierColor: modifier edge color
        - boundaryColor: boundary node color
        - nodeEdgeColor: node edge color
        - nodeEdgelw: linewidth of node edges
        - highlight: list of species ids or reaction ids to highlight
        - hlNodeColor: node color of highlighted nodes
        - hlNodeEdgeColor: node edge color of highlighted nodes
        - drawReactionNode: flag for drawing reaction nodes
        - breakBoundary: flag for breaking all boundary species into separate nodes
    
    .. method:: Network.draw()
    
        Draw weighted reaction network based on frequency of reactions
        
    .. method:: Network.getLayout()
    
        Return the layout
        
    .. method:: Network.reset()
    
        Resets all propertiess


-------------
Miscellaneous
-------------

.. autofunction:: netplotlib.netplotlib.getVersion

