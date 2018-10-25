# -*- coding: utf-8 -*-
"""
Simple module to visualize reaction networks and weighted reaction networks
from ensemble of models.

Kiri Choi (c) 2018
"""

import tellurium as te
import antimony
import networkx as nx
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle, FancyBboxPatch
import matplotlib.pyplot as plt
import numpy as np
import sympy
import tesbml
import itertools
from testmodels import testmodels as ts

def plotNetworkFromAntimony(model, scale=1.5, fontsize=20, lw=3, node='tab:blue',
                        reaction='tab:blue', label='w', edge='k', 
                        modifier='tab:red', boundary='tab:gray',
                        break_boundary=False):
    """
    plot reaction network from an Antimony string
    
    :param model: Antimony string of a model to plot
    :param scale: scaling factor for layout algorithm
    :param fontsize: fontsize for labels
    :param lw: linewidth of edges
    :param node: node color
    :param reaction: reaction node color
    :param label: label color
    :param edge: edge color
    :param modifier: modifier edge color
    :param boundary: boundary node color
    :param break_boundary: flag for breaking all boundary species into separate nodes
    """
    
    r = te.loada(model)
    
    plotNetworkFromSBML(r.getSBML(), scale=scale, fontsize=fontsize, lw=lw, 
                node=node, reaction=reaction, label=label, edge=edge,
                modifier=modifier, boundary=boundary, break_boundary=break_boundary)
    

def plotNetworkFromSBML(model, scale=1.5, fontsize=20, lw=3, node='tab:blue',
                reaction='tab:blue', label='w', edge='k', modifier='tab:red', 
                boundary='tab:gray', break_boundary=False):
    """     
    plot reaction network from an SBML string
    
    :param model: SBML string of a model to plot
    :param scale: scaling factor for layout algorithm
    :param fontsize: fontsize for labels
    :param lw: linewidth of edges
    :param node: node color
    :param reaction: reaction node color
    :param label: label color
    :param edge: edge color
    :param modifier: modifier edge color
    :param boundary: boundary node color
    :param break_boundary: flag for breaking all boundary species into separate nodes
    """
    
    r = te.loadSBMLModel(model)
    numBnd = r.getNumBoundarySpecies()
    numFlt = r.getNumFloatingSpecies()
    boundaryId = r.getBoundarySpeciesIds()
    floatingId = r.getFloatingSpeciesIds()
    speciesId = boundaryId + floatingId
    rid = r.getReactionIds()
    
    # prepare symbols for sympy
    paramIdsStr = ' '.join(r.getGlobalParameterIds())
    floatingIdsStr = ' '.join(r.getFloatingSpeciesIds())
    boundaryIdsStr = ' '.join(r.getBoundarySpeciesIds())
    comparmentIdsStr = ' '.join(r.getCompartmentIds())
    
    allIds = paramIdsStr + ' ' + floatingIdsStr + ' ' + boundaryIdsStr + ' ' + comparmentIdsStr
    
    avsym = sympy.symbols(allIds)
    
    # extract reactant, product, modifiers, and kinetic laws
    rct = []
    prd = []
    mod = []
    mod_target = []
    kineticLaw = []
    mod_type = []
    
    doc = tesbml.readSBMLFromString(r.getSBML())
    sbmlmodel = doc.getModel()

    for slr in sbmlmodel.getListOfReactions():
        temprct = []
        tempprd = []
        tempmod = []
        
        sbmlreaction = sbmlmodel.getReaction(slr.getId())
        for sr in range(sbmlreaction.getNumReactants()):
            sbmlrct = sbmlreaction.getReactant(sr)
            temprct.append(sbmlrct.getSpecies())
        for sp in range(sbmlreaction.getNumProducts()):
            sbmlprd = sbmlreaction.getProduct(sp)
            tempprd.append(sbmlprd.getSpecies())
        for sm in range(sbmlreaction.getNumModifiers()):
            sbmlmod = sbmlreaction.getModifier(sm)
            tempmod.append(sbmlmod.getSpecies())
        kl = sbmlreaction.getKineticLaw()
        
        rct.append(temprct)
        prd.append(tempprd)
        mod.append(tempmod)
        kineticLaw.append(kl.getFormula())
    
    # use sympy for analyzing modifiers weSmart
    for ml in range(len(mod)):
        mod_type_temp = []
        expression = kineticLaw[ml]
        n,d = sympy.fraction(expression)
        for ml_i in range(len(mod[ml])):
            if n.has(mod[ml][ml_i]):
                mod_type_temp.append('activator')
            elif d.has(mod[ml][ml_i]):
                mod_type_temp.append('inhibitor')
            else:
                continue
        mod_type.append(mod_type_temp)
    
    for i in range(len(mod)):
        mod_target_temp = []
        if len(mod[i]) > 0:
            mod_target_temp.append(rid[i])
        mod_target.append(mod_target_temp)
    
    mod_flat = [item for sublist in mod for item in sublist]
    modtype_flat = [item for sublist in mod_type for item in sublist]
    modtarget_flat = [item for sublist in mod_target for item in sublist]
    
    # initialize directional graph
    G = nx.DiGraph()

    # add edges
    # TODO: Separate boundary species?
    for i in range(sbmlmodel.getNumReactions()):
        for k in range(len(rct[i])):
            G.add_edges_from([(rct[i][k], rid[i])], weight=(1+lw))
        for j in range(len(prd[i])):
            G.add_edges_from([(rid[i], prd[i][j])], weight=(1+lw))
                    
        if len(mod[i]) > 0:
            if mod_type[i][0] == 'inhibitor':
                G.add_edges_from([(mod[i][0], rid[i])], weight=(1+lw))
            elif mod_type[i][0] == 'activator':
                G.add_edges_from([(mod[i][0], rid[i])], weight=(1+lw))
        
    # calcutate positions
    thres = 0.1
    shortest_dist = dict(nx.shortest_path_length(G, weight='weight'))
    pos = nx.kamada_kawai_layout(G, dist=shortest_dist, scale=scale)
    
    dist_flag = True
    
    while dist_flag:
        dist_flag = False
        for i in itertools.combinations(pos.keys(), 2):
            pos_dist = np.linalg.norm(pos[i[0]] - pos[i[1]])
            if pos_dist < thres:
                dist_flag = True
                shortest_dist[i[0]][i[1]] = 4
        pos = nx.kamada_kawai_layout(G, dist=shortest_dist, scale=scale)
    
    # check the range of x and y positions
    max_width = []
    max_height = []
    for key, value in pos.items():
        max_width.append(value[0])
        max_height.append(value[1])
    
    max_width = [min(max_width), max(max_width)]
    max_height = [min(max_height), max(max_height)]
    
    # initialize figure
    fig = plt.figure()
    ax = plt.gca()
    
    # add nodes to the figure
    for n in G:
        if n in rid:
            rec_rad = max(0.01*(len(n)+2), 0.055)
            c = Circle(pos[n], radius=rec_rad, color=reaction)
            plt.text(pos[n][0], pos[n][1], n, fontsize=fontsize, 
                 horizontalalignment='center', verticalalignment='center', color=label)
        else:
            # TODO: if the label is too long, increase the height and change line/abbreviate?
            rec_width = max(0.04*(len(n)+2), 0.17)
            rec_height = 0.12
            if n in boundaryId:
                node_color = boundary
            else:
                node_color = node
            c = FancyBboxPatch((pos[n][0]-rec_width/2, pos[n][1]-rec_height/2),
                               rec_width, rec_height,
                               boxstyle="round,pad=0.01, rounding_size=0.02",
                               linewidth=0, color=node_color)
            plt.text(pos[n][0], pos[n][1], n, 
                     fontsize=fontsize, horizontalalignment='center', 
                     verticalalignment='center', color=label)
        ax.add_patch(c)
        G.node[n]['patch'] = c
    
    # add edges to the figure
    seen={}
    for (u,v,d) in G.edges(data=True):
        n1 = G.node[u]['patch']
        n2 = G.node[v]['patch']
        rad = 0.1
        if (u,v) in seen:
            rad = seen.get((u,v)) # TODO: No curvature when there is just a single line between two nodes
            rad = (rad+np.sign(rad)*0.1)*-1 # TODO: Change curvature
        
        if u in rid or v in rid:
            if u in rid:
                X1 = (n1.center[0],n1.center[1])
            else:
                X1 = (n1.get_x()+n1.get_width()/2,n1.get_y()+n1.get_height()/2)
            if v in rid:
                X2 = (n2.center[0],n2.center[1])
            else:
                X2 = (n2.get_x()+n2.get_width()/2,n2.get_y()+n2.get_height()/2)
        else:
            X1 = (n1.get_x()+n1.get_width()/2,n1.get_y()+n1.get_height()/2)
            X2 = (n2.get_x()+n2.get_width()/2,n2.get_y()+n2.get_height()/2)
        
        if u not in rid and v in rid and u not in mod_flat and v not in modtarget_flat: # species node to reaction node
            color=edge
            arrowstyle='-'
        elif u not in rid and v in rid and u in mod_flat and v in modtarget_flat: # modifiers
            uind = [i for i, e in enumerate(mod_flat) if e == u]
            vind = [i for i, e in enumerate(modtarget_flat) if e == v]
            if modtype_flat[list(set(uind).intersection(vind))[0]] == 'inhibitor': # inhibition
                color=modifier
                arrowstyle='-['
            else: # activation
                color=modifier
                arrowstyle='-|>'
        else: # reaction node to species node
            color=edge
            arrowstyle='-|>'
        e = FancyArrowPatch(X1,
                            X2,
                            patchA=n1,
                            patchB=n2,
                            arrowstyle=arrowstyle,
                            connectionstyle='arc3,rad=%s'%rad,
                            mutation_scale=10.0,
                            lw=G[u][v]['weight'],
                            color=color)
        seen[(u,v)]=rad
        ax.add_patch(e)
    
    # reset width and height
    ax.autoscale()
    fig.set_figwidth((abs(max_width[0] - max_width[1])+0.5)*5)
    fig.set_figheight((abs(max_height[0] - max_height[1])+0.5)*5)
    plt.axis('off')
    plt.axis('equal')
    
    plt.show()

def plotWeightedNetworkFromAntimony(models, scale=1.5, fontsize=20, lw=10, 
                node='tab:blue', reaction='tab:blue', label='w', edge='k', 
                modifier='tab:red', boundary='tab:gray', break_boundary=False):
    """     
    plot weighted reaction network from a list of Antimony strings
    
    :param models: a list of antimony string of models to plot
    :param scale: scaling factor for layout algorithm
    :param fontsize: fontsize for labels
    :param lw: weighting constant for the thickness of the edges
    :param node: node color
    :param reaction: reaction node color
    :param label: label color
    :param edge: edge color
    :param modifier: modifier edge color
    :param boundary: boundary node color
    :param break_boundary: flag for breaking all boundary species into separate nodes
    :returns allRxn: list of pairs of all reactions present throughout the ensemble
    :returns count: weighted frequency of each reaction
    """
    
    for i in models:
        r = te.loada(i)
        models[i] = r.getSBML()
        
    plotWeightedNetworkFromSBML(models, scale=scale, fontsize=fontsize, lw=lw, 
                node=node, reaction=reaction, label=label, edge=edge,
                modifier=modifier, boundary=boundary, break_boundary=break_boundary)
    
    
def plotWeightedNetworkFromSBML(models, scale=1.5, fontsize=20, lw=10, node='tab:blue',
                reaction='tab:blue', label='w', edge='k', modifier='tab:red', 
                boundary='tab:gray', break_boundary=False):
    """     
    plot weighted reaction network from a list of SBML strings
    
    :param models: a list of antimony string of models to plot
    :param scale: scaling factor for layout algorithm
    :param fontsize: fontsize for labels
    :param lw: weighting constant for the thickness of the edges
    :param node: node color
    :param reaction: reaction node color
    :param label: label color
    :param edge: edge color
    :param modifier: modifier edge color
    :param boundary: boundary node color
    :param break_boundary: flag for breaking all boundary species into separate nodes
    :returns allRxn: list of pairs of all reactions present throughout the ensemble
    :returns count: weighted frequency of each reaction
    """
    
    # extract reactant, product, modifiers, and kinetic laws
    allRxn = []
    count = []
    rid = []
    mod = []
    mod_target = []
    mod_type = []
    rid_ind = 0

    for m in models:
        rct = []
        prd = []
        mod_m = []
        mod_target_m = []
        kineticLaw = []
        mod_type_m = []
        
        r = te.loadSBMLModel(m)
        numBnd = r.getNumBoundarySpecies()
        numFlt = r.getNumFloatingSpecies()
        boundaryId = r.getBoundarySpeciesIds()
        floatingId = r.getFloatingSpeciesIds()
        speciesId = boundaryId + floatingId
        
        # prepare symbols for sympy
        paramIdsStr = ' '.join(r.getGlobalParameterIds())
        floatingIdsStr = ' '.join(r.getFloatingSpeciesIds())
        boundaryIdsStr = ' '.join(r.getBoundarySpeciesIds())
        comparmentIdsStr = ' '.join(r.getCompartmentIds())
        
        allIds = paramIdsStr + ' ' + floatingIdsStr + ' ' + boundaryIdsStr + ' ' + comparmentIdsStr
        
        avsym = sympy.symbols(allIds)
        
        doc = tesbml.readSBMLFromString(r.getSBML())
        sbmlmodel = doc.getModel()
    
        for slr in sbmlmodel.getListOfReactions():
            temprct = []
            tempprd = []
            tempmod = []
            
            sbmlreaction = sbmlmodel.getReaction(slr.getId())
            for sr in range(sbmlreaction.getNumReactants()):
                sbmlrct = sbmlreaction.getReactant(sr)
                temprct.append(sbmlrct.getSpecies())
            for sp in range(sbmlreaction.getNumProducts()):
                sbmlprd = sbmlreaction.getProduct(sp)
                tempprd.append(sbmlprd.getSpecies())
            for sm in range(sbmlreaction.getNumModifiers()):
                sbmlmod = sbmlreaction.getModifier(sm)
                tempmod.append(sbmlmod.getSpecies())
            kl = sbmlreaction.getKineticLaw()
            
            rct.append(temprct)
            prd.append(tempprd)
            mod_m.append(tempmod)
            kineticLaw.append(kl.getFormula())
        
        # use sympy for analyzing modifiers weSmart
        for ml in range(len(mod_m)):
            mod_type_temp = []
            expression = kineticLaw[ml]
            n,d = sympy.fraction(expression)
            for ml_i in range(len(mod_m[ml])):
                if n.has(mod_m[ml][ml_i]):
                    mod_type_temp.append('activator')
                elif d.has(mod_m[ml][ml_i]):
                    mod_type_temp.append('inhibitor')
                else:
                    continue
            mod_type_m.append(mod_type_temp)
        
        for i in range(len(mod_m)):
            mod_target_temp = []
            if len(mod_m[i]) > 0:
                mod_target_temp.append(rid[i])
            mod_target_m.append(mod_target_temp)
            
        mod_flat = [item for sublist in mod_m for item in sublist]
        modtype_flat = [item for sublist in mod_type_m for item in sublist]
        modtarget_flat = [item for sublist in mod_target_m for item in sublist]
        
        for t in range(sbmlmodel.getNumReactions()):
            if [rct[t], prd[t]] not in allRxn:
                allRxn.append([rct[t], prd[t]])
                count.append(1)
                rid.append("J" + str(rid_ind))
                mod.append(mod_flat)
                mod_type.append(modtype_flat)
                mod_target.append(modtarget_flat)
                rid_ind += 1
            else:
                count[allRxn.index([rct[t], prd[t]])] += 1
    
    count = count/np.sum(len(models))

    # initialize directional graph
    G = nx.DiGraph()

    # add edges
    # TODO: Separate boundary species?
    for i in range(len(allRxn)):
        for k in range(len(allRxn[i][0])):
            G.add_edges_from([(allRxn[i][0][k], rid[i])], weight=(count[i]*lw))
        for j in range(len(allRxn[i][1])):
            G.add_edges_from([(rid[i], allRxn[i][1][j])], weight=(count[i]*lw))
                    
        if len(mod[i]) > 0:
            if mod_type[i][0] == 'inhibitor':
                G.add_edges_from([(mod[i][0], rid[i])], weight=(count[i]*lw))
            elif mod_type[i][0] == 'activator':
                G.add_edges_from([(mod[i][0], rid[i])], weight=(count[i]*lw))

    # calcutate positions
    thres = 0.1
    shortest_dist = dict(nx.shortest_path_length(G, weight='weight'))
    pos = nx.kamada_kawai_layout(G, dist=shortest_dist, scale=scale)
    
    dist_flag = True
    
    while dist_flag:
        dist_flag = False
        for i in itertools.combinations(pos.keys(), 2):
            pos_dist = np.linalg.norm(pos[i[0]] - pos[i[1]])
            if pos_dist < thres:
                dist_flag = True
                shortest_dist[i[0]][i[1]] = 4
        pos = nx.kamada_kawai_layout(G, dist=shortest_dist, scale=scale)
        
    # check the range of x and y positions
    max_width = []
    max_height = []
    for key, value in pos.items():
        max_width.append(value[0])
        max_height.append(value[1])
    
    max_width = [min(max_width), max(max_width)]
    max_height = [min(max_height), max(max_height)]
    
    # initialize figure
    fig = plt.figure()
    ax = plt.gca()
    
    # add nodes to the figure
    for n in G:
        if n in rid:
            rec_rad = max(0.01*(len(n)+2), 0.055)
            c = Circle(pos[n], radius=rec_rad, color=reaction)
            plt.text(pos[n][0], pos[n][1], n, fontsize=fontsize, 
                 horizontalalignment='center', verticalalignment='center', color=label)
        else:
            # TODO: if the label is too long, increase the height and change line/abbreviate?
            rec_width = max(0.04*(len(n)+2), 0.17)
            rec_height = 0.12
            if n in boundaryId:
                node_color = boundary
            else:
                node_color = node
            c = FancyBboxPatch((pos[n][0]-rec_width/2, pos[n][1]-rec_height/2),
                               rec_width, rec_height,
                               boxstyle="round,pad=0.01, rounding_size=0.02",
                               linewidth=0, color=node_color)
            plt.text(pos[n][0], pos[n][1], n, 
                     fontsize=fontsize, horizontalalignment='center', 
                     verticalalignment='center', color=label)
        ax.add_patch(c)
        G.node[n]['patch'] = c
    
    # add edges to the figure
    seen={}
    for (u,v,d) in G.edges(data=True):
        n1 = G.node[u]['patch']
        n2 = G.node[v]['patch']
        rad = 0.1
        if (u,v) in seen:
            rad = seen.get((u,v)) # TODO: No curvature when there is just a single line between two nodes
            rad = (rad+np.sign(rad)*0.1)*-1 # TODO: Change curvature
        
        if u in rid or v in rid:
            if u in rid:
                X1 = (n1.center[0],n1.center[1])
            else:
                X1 = (n1.get_x()+n1.get_width()/2,n1.get_y()+n1.get_height()/2)
            if v in rid:
                X2 = (n2.center[0],n2.center[1])
            else:
                X2 = (n2.get_x()+n2.get_width()/2,n2.get_y()+n2.get_height()/2)
        else:
            X1 = (n1.get_x()+n1.get_width()/2,n1.get_y()+n1.get_height()/2)
            X2 = (n2.get_x()+n2.get_width()/2,n2.get_y()+n2.get_height()/2)
        
        if u not in rid and v in rid and u not in mod_flat and v not in modtarget_flat: # species node to reaction node
            color=edge
            arrowstyle='-'
        elif u not in rid and v in rid and u in mod_flat and v in modtarget_flat: # modifiers
            uind = [i for i, e in enumerate(mod_flat) if e == u]
            vind = [i for i, e in enumerate(modtarget_flat) if e == v]
            if modtype_flat[list(set(uind).intersection(vind))[0]] == 'inhibitor': # inhibition
                color=modifier
                arrowstyle='-['
            else: # activation
                color=modifier
                arrowstyle='-|>'
        else: # reaction node to species node
            color=edge
            arrowstyle='-|>'
        e = FancyArrowPatch(X1,
                            X2,
                            patchA=n1,
                            patchB=n2,
                            arrowstyle=arrowstyle,
                            connectionstyle='arc3,rad=%s'%rad,
                            mutation_scale=10.0,
                            lw=G[u][v]['weight'],
                            alpha=G[u][v]['weight']/lw,
                            color=color)
        seen[(u,v)]=rad
        ax.add_patch(e)

    # Edge labels
    edgeLabels = {}
    for i in range(len(allRxn)):
        for j in range(len(allRxn[i][1])):
            edgeLabels[(rid[i], allRxn[i][1][j])] = round(count[i], 3)
            
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edgeLabels, font_size=12)

    # reset width and height
    ax.autoscale()
    fig.set_figwidth((abs(max_width[0] - max_width[1])+0.5)*5)
    fig.set_figheight((abs(max_height[0] - max_height[1])+0.5)*5)
    plt.axis('off')
    plt.axis('equal')
    
    plt.show()
    
    return allRxn, count

