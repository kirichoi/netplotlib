# -*- coding: utf-8 -*-
"""
Simple module to visualize reaction networks and weighted reaction networks
from ensemble of models.

Kiri Choi (c) 2018
"""

import tellurium as te
import antimony
import networkx as nx
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.pyplot as plt
import numpy as np

def plotNetworkFromSBML(model):
    """
    plot reaction network from a single SBML model
    
    :param model: SBML file path or string    
    """
    
    r = te.loadSBMLModel(model)
    
    plotNetwork(r.getAntimony())
    

def plotNetwork(model, lw=3):
    """     
    plot reaction network from a single model
    
    :param model: antimony string of a model to plot
    """
    
    r = te.loada(model)
    numBnd = r.getNumBoundarySpecies()
    numFlt = r.getNumFloatingSpecies()
    
    antimony.loadAntimonyString(model)
    module = antimony.getModuleNames()[-1]
    rct = np.array(antimony.getReactantNames(module)).tolist()
    prd = np.array(antimony.getProductNames(module)).tolist()
    
    G = nx.DiGraph()

    for i in range(antimony.getNumReactions(module)):
        if len(rct[i]) == 1:
            if len(prd[i]) == 1:
                G.add_edges_from([(rct[i][0], prd[i][0])], weight=(1+lw))
            else:
                G.add_edges_from([(rct[i][0], 'J')], weight=(1+lw))
                for j in range(len(prd[i])):
                    G.add_edges_from([('J', prd[i][j])], weight=(1+lw))
        else:
            if len(prd[i]) == 1:
                for k in range(len(rct[i])):
                    G.add_edges_from([(rct[i][k], 'J')], weight=(1+lw))
                G.add_edges_from([('J', prd[i][0])], weight=(1+lw))
            else:
                for k in range(len(rct[i])):
                    G.add_edges_from([(rct[i][k], 'J1')], weight=(1+lw))
                G.add_edges_from([("J1", 'J2')], weight=(1+lw))
                for j in range(len(prd[i])):
                    G.add_edges_from([('J2', prd[i][j])], weight=(1+lw))
    
    pos = nx.spring_layout(G)
    
    # TODO: Figure out the maximum position
    
#    edges = G.edges()
#    weights = [G[u][v]['weight'] for u,v in edges]
    
    flat_rct = [item for sublist in rct for item in sublist]
    flat_prd = [item for sublist in prd for item in sublist]

    strMaxLen = len(max(flat_rct+flat_prd, key=len))
#    nx.draw(G, pos, node_size=250*strMaxLen, with_labels=True, width=weights, 
#            node_shape='s')
    
    fontsize = 20
    
    fig = plt.figure()
    ax = plt.gca()
    
    for n in G:
        c = Rectangle(pos[n], width=0.1*(strMaxLen+strMaxLen*0.1), height=0.15)
        ax.add_patch(c)
        G.node[n]['patch'] = c
        plt.text(pos[n][0]+c.get_width()/2, pos[n][1]+c.get_height()/2, n, fontsize=fontsize, 
                 horizontalalignment='center', verticalalignment='center')
    seen={}
    for (u,v,d) in G.edges(data=True):
        n1 = G.node[u]['patch']
        n2 = G.node[v]['patch']
        rad = 0.1
        if (u,v) in seen:
            rad = seen.get((u,v))
            rad = (rad+np.sign(rad)*0.1)*-1
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
    fig.set_figheight(10)
    fig.set_figwidth(20)
    plt.axis('off')
    plt.axis('equal')
    
    
    plt.show()
    
    antimony.clearPreviousLoads()
    
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

    pos = nx.kamada_kawai_layout(G)
    
    edges = G.edges()
    weights = [G[u][v]['weight'] for u,v in edges]
    
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
    
    nx.draw(G, pos, node_size=500, with_labels=True, width=weights)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edgeLabels, label_pos=0.3)
    
    plt.show()
    
    return allRxn, count