# -*- coding: utf-8 -*-
"""
Simple module to visualize reaction networks and weighted reaction networks
from ensemble of models.

Kiri Choi (c) 2018
"""

import tellurium as te
import antimony
import networkx as nx
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle
import matplotlib.pyplot as plt
import numpy as np
import sympy
import tesbml
from test_models import testmodels as ts

def plotNetworkFromSBML(model):
    """
    plot reaction network from a single SBML model
    
    :param model: SBML file path or string    
    """
    
    r = te.loadSBMLModel(model)
    
    plotNetwork(r.getAntimony())
    

def plotNetwork(model, scale=1.5, fontsize=20, lw=3):
    """     
    plot reaction network from a single model
    
    :param model: antimony string of a model to plot
    :param scale: scaling factor for layout algorithm
    :param fontsize: fontsize for labels
    :param lw: linewidth of edges
    :param node: node color
    :param reaction: reaction node color
    :param label: label color
    :param edge: edge color
    :param modifier: modifier edge color
    """
    
    r = te.loada(model)
    numBnd = r.getNumBoundarySpecies()
    numFlt = r.getNumFloatingSpecies()
    rid = r.getReactionIds()
    
    paramIdsStr = ' '.join(r.getGlobalParameterIds())
    floatingIdsStr = ' '.join(r.getFloatingSpeciesIds())
    boundaryIdsStr = ' '.join(r.getBoundarySpeciesIds())
    comparmentIdsStr = ' '.join(r.getCompartmentIds())
    
    allIds = paramIdsStr + ' ' + floatingIdsStr + ' ' + boundaryIdsStr + ' ' + comparmentIdsStr
    
    avsym = sympy.symbols(allIds)
    
    rct = []
    prd = []
    mod = []
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
    
    # Sympy weSmart
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
    
    G = nx.DiGraph()

    # TODO: Support modifiers
    # TODO: Separate boundary species?
    for i in range(sbmlmodel.getNumReactions()):
        if len(rct[i]) == 1:
            if len(prd[i]) == 1:
                G.add_edges_from([(rct[i][0], rid[i])], weight=(1+lw))
                G.add_edges_from([(rid[i], prd[i][0])], weight=(1+lw))
            else:
                G.add_edges_from([(rct[i][0], rid[i])], weight=(1+lw))
                for j in range(len(prd[i])):
                    G.add_edges_from([(rid[i], prd[i][j])], weight=(1+lw))
        else:
            if len(prd[i]) == 1:
                for k in range(len(rct[i])):
                    G.add_edges_from([(rct[i][k], rid[i])], weight=(1+lw))
                G.add_edges_from([(rid[i], prd[i][0])], weight=(1+lw))
            else:
                for k in range(len(rct[i])):
                    G.add_edges_from([(rct[i][k], rid[i])], weight=(1+lw))
                for j in range(len(prd[i])):
                    G.add_edges_from([(rid[i], prd[i][j])], weight=(1+lw))
                    
        if len(mod[i]) > 0:
            if mod_type[i][0] == 'inhibitor':
                G.add_edges_from([(mod[i][0], rid[i])], weight=(1+lw))
            elif mod_type[i][0] == 'activator':
                G.add_edges_from([(mod[i][0], rid[i])], weight=(1+lw))
    
    pos = nx.kamada_kawai_layout(G, scale=scale)
    
    max_width = []
    max_height = []
    for key, value in pos.items():
        max_width.append(value[0])
        max_height.append(value[1])
    
    max_width = [min(max_width), max(max_width)]
    max_height = [min(max_height), max(max_height)]
    
    flat_rct = [item for sublist in rct for item in sublist]
    flat_prd = [item for sublist in prd for item in sublist]

    strMaxLen = len(max(flat_rct+flat_prd, key=len))
    
    fig = plt.figure()
    ax = plt.gca()
    
    for n in G:
        if n in rid:
            c = Circle(pos[n], radius=0.05)
            plt.text(pos[n][0], pos[n][1], n, fontsize=fontsize, 
                 horizontalalignment='center', verticalalignment='center')
        else:
            rec_width = max(0.05*strMaxLen, 0.2)
            rec_height = 0.15
            c = Rectangle(np.array([pos[n][0]-rec_width/2,pos[n][1]-rec_height/2]), width=rec_width, height=rec_height)
            plt.text(pos[n][0], pos[n][1], n, 
                     fontsize=fontsize, horizontalalignment='center', 
                     verticalalignment='center')
        ax.add_patch(c)
        G.node[n]['patch'] = c
       
    seen={}
    for (u,v,d) in G.edges(data=True):
        n1 = G.node[u]['patch']
        n2 = G.node[v]['patch']
        rad = 0.1
        if (u,v) in seen:
            rad = seen.get((u,v)) # TODO: No curvature when there is just a single line between two nodes
            rad = (rad+np.sign(rad)*0.1)*-1 # TODO: Change curvature
        color = 'k'
        
        if u in rid or v in rid:
            if u in rid:
                X1 = (n1.center[0],n1.center[1])
            else:
                X1 = (n1.xy[0]+n1.get_width()/2,n1.xy[1]+n1.get_height()/2)
            if v in rid:
                X2 = (n2.center[0],n2.center[1])
            else:
                X2 = (n2.xy[0]+n2.get_width()/2,n2.xy[1]+n2.get_height()/2)
        else:
            X1 = (n1.xy[0]+n1.get_width()/2,n1.xy[1]+n1.get_height()/2)
            X2 = (n2.xy[0]+n2.get_width()/2,n2.xy[1]+n2.get_height()/2)
        e = FancyArrowPatch(X1,
                            X2,
                            patchA=n1,
                            patchB=n2,
                            arrowstyle='-|>',
                            connectionstyle='arc3,rad=%s'%rad,
                            mutation_scale=10.0,
                            lw=G[u][v]['weight'],
                            color=color)
        seen[(u,v)]=rad
        ax.add_patch(e)
    
    ax.autoscale()
    fig.set_figwidth((abs(max_width[0] - max_width[1])+0.5)*5)
    fig.set_figheight((abs(max_height[0] - max_height[1])+0.5)*5)
    plt.axis('off')
    plt.axis('equal')
    
    plt.show()
    
    
def plotWeightedNetwork(models, lw=10):
    """     
    plot weighted reaction network from an ensemble of model
    
    :param model: a list of antimony string of models to plot
    :param lw: weighting constant for the thickness of the edges
    :returns allRxn: list of pairs of all reactions present throughout the ensemble
    :returns count: weighted frequency of each reaction
    """
    
    allRxn = []
    count = []

    for m in models:
        antimony.loadAntimonyString(m)
        module = antimony.getModuleNames()[-1]
        rct = np.array(antimony.getReactantNames(module)).tolist()
        prd = np.array(antimony.getProductNames(module)).tolist()
        
        for t in range(antimony.getNumReactions(module)):
            if [rct[t], prd[t]] not in allRxn:
                allRxn.append([rct[t], prd[t]])
                count.append(1)
            else:
                count[allRxn.index([rct[t], prd[t]])] += 1
        
        antimony.clearPreviousLoads()

    count = count/np.sum(len(models))

    G = nx.DiGraph()

    for i in range(len(allRxn)):
        if len(allRxn[i][0]) == 1:
            if len(allRxn[i][1]) == 1:
                G.add_edges_from([(allRxn[i][0][0], allRxn[i][1][0])], weight=count[i]*lw)
            else:
                G.add_edges_from([(allRxn[i][0], 'J')], weight=count[i]*lw)
                for j in range(len(allRxn[i][1])):
                    G.add_edges_from([('J', allRxn[i][1][j])], weight=count[i]*lw)
        else:
            if len(allRxn[i][1]) == 1:
                for k in range(len(allRxn[i][0])):
                    G.add_edges_from([(allRxn[i][0][k], 'J')], weight=count[i]*lw)
                G.add_edges_from([('J', allRxn[i][1][0])], weight=count[i]*lw)
            else:
                for k in range(len(allRxn[i][0])):
                    G.add_edges_from([(allRxn[i][0][k], 'J')], weight=count[i]*lw)
                for j in range(len(allRxn[i][1])):
                    G.add_edges_from([('J', allRxn[i][1][j])], weight=count[i]*lw)

    pos = nx.spring_layout(G)
    
    max_width = []
    max_height = []
    for key, value in pos.items():
        max_width.append(value[0])
        max_height.append(value[1])
    
    max_width = [min(max_width), max(max_width)]
    max_height = [min(max_height), max(max_height)]
    
#    edges = G.edges()
#    weights = [G[u][v]['weight'] for u,v in edges]
    
    edgeLabels = {}
    for i in range(len(allRxn)):
        if len(allRxn[i][0]) == 1:
            if len(allRxn[i][1]) == 1:
                edgeLabels[(allRxn[i][0][0], allRxn[i][1][0])] = round(count[i], 3)
            else:
                edgeLabels[(allRxn[i][0], 'J')] = round(count[i], 3)
                for j in range(len(allRxn[i][1])):
                    edgeLabels[('J', allRxn[i][1][j])] = round(count[i], 3)
        else:
            if len(allRxn[i][1]) == 1:
                for k in range(len(allRxn[i][0])):
                    edgeLabels[(allRxn[i][0][k], 'J')] = round(count[i], 3)
                edgeLabels[('J', allRxn[i][1][0])] = round(count[i], 3)
            else:
                for k in range(len(allRxn[i][0])):
                    edgeLabels[(allRxn[i][0][k], 'J')] = round(count[i], 3)
                for j in range(len(allRxn[i][1])):
                    edgeLabels[('J', allRxn[i][1][j])] = round(count[i], 3)
    
    flat_rct = [item for sublist in rct for item in sublist]
    flat_prd = [item for sublist in prd for item in sublist]

    strMaxLen = len(max(flat_rct+flat_prd, key=len))
#    nx.draw(G, pos, node_size=250*strMaxLen, with_labels=True, width=weights, 
#            node_shape='s')
    
    
    fontsize = 20
    
    fig = plt.figure()
    ax = plt.gca()
    
    for n in G:
        if n == '':
            c = Rectangle(pos[n], width=0, height=0)
        else:
            c = Rectangle(pos[n], width=max(0.05*strMaxLen, 0.2), height=0.15)
        ax.add_patch(c)
        G.node[n]['patch'] = c
        plt.text(pos[n][0], pos[n][1], n, fontsize=fontsize, 
                 horizontalalignment='center', verticalalignment='center')
    seen={}
    for (u,v,d) in G.edges(data=True):
        n1 = G.node[u]['patch']
        n2 = G.node[v]['patch']
        rad = 0.1
        if (u,v) in seen:
            rad = seen.get((u,v))
            rad = (rad+np.sign(rad)*0.1)*-1 # TODO: Change curvature
        color = 'k'

        e = FancyArrowPatch((n1.xy[0]+n1.get_width()/2,n1.xy[1]+n1.get_height()/2),
                            (n2.xy[0]+n2.get_width()/2,n2.xy[1]+n2.get_height()/2),
                            patchA=n1,
                            patchB=n2,
                            arrowstyle='-|>',
                            connectionstyle='arc3,rad=%s'%rad,
                            mutation_scale=10.0,
                            lw=G[u][v]['weight'],
                            color=color)
        seen[(u,v)]=rad
        ax.add_patch(e)
    
    ax.autoscale()
    fig.set_figwidth((abs(max_width[0] - max_width[1])+0.5)*5)
    fig.set_figheight((abs(max_height[0] - max_height[1])+0.5)*5)
    plt.axis('off')
    plt.axis('equal')
    
    
#    nx.draw(G, pos, node_size=500, with_labels=True, width=weights)
#    nx.draw_networkx_edge_labels(G, pos, edge_labels=edgeLabels, label_pos=0.3)
    
    plt.show()
    
    return allRxn, count

