# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os, re
import tellurium as te
import networkx as nx
from matplotlib.patches import (FancyArrowPatch, FancyBboxPatch, ArrowStyle, 
                                Arc, RegularPolygon, Rectangle)
from matplotlib.path import Path
import matplotlib.pyplot as plt
from matplotlib import gridspec, cm, colors
import numpy as np
from scipy import interpolate
import sympy
import itertools
import layout
import toolbox
import libsbml

def getVersion():
    """
    Print version string
    """
    
    try:
    	with open(os.path.join(os.path.dirname(__file__), '..', 'VERSION.txt'), 'r') as f:
    		version = f.read().rstrip()
    except:
    	with open(os.path.join(os.path.dirname(__file__), 'VERSION.txt'), 'r') as f:
    		version = f.read().rstrip()
    
    return version


class _Variable():
    
    def __init__(self):
        pass


class Network():
    
    def __init__(self, model):
        """
        Creates a new Network object. 
        
        :param model: SBML or Antimony string of a model
        :type name: str
        """
        
        try:
            self.rrInstance = te.loadSBMLModel(model)
        except:
            try:
                self.rrInstance = te.loadAntimonyModel(model)
            except:
                raise Exception("Input does not seem to be a valid SBML or Antimony string")
            
        self._Var = _Variable()
            
        doc = libsbml.readSBMLFromString(self.rrInstance.getSBML())
        self._Var.sbmlmodel = doc.getModel()
        self._Var.layoutPlugin = self._Var.sbmlmodel.getPlugin('layout')
        if ((self._Var.layoutPlugin != None) and (self._Var.layoutPlugin.getNumLayouts() > 0)):
            self._Var.layoutAvail = True
        else:
            self._Var.layoutAvail = False
        
        try:
            self._Var.boundaryId = self.rrInstance.getBoundarySpeciesIds()
            self._Var.floatingId = self.rrInstance.getFloatingSpeciesIds()
            self._Var.compartmentId = self.rrInstance.getCompartmentIds()
            self._Var.rid = self.rrInstance.getReactionIds()
            self._Var.stoch = self.rrInstance.getFullStoichiometryMatrix()
            self._Var.stoch_row = self._Var.stoch.rownames
        except:
            raise Exception("Failed to analyze the file: Check the file is valid")
        self.reset()
        

    def reset(self):
        """
        Resets all properties
        """
    
        self.scale = 1.
        self.fontsize = 10
        self.edgelw = 1.
        self.nodeColor = 'tab:blue'
        self.reactionNodeColor = 'tab:gray'
        self.labelColor = 'w'
        self.labelReactionIds = False
        self.reactionColor = 'k'
        self.modifierColor = 'tab:red'
        self.boundaryColor = 'tab:green'
        self.nodeEdgeColor = 'k'
        self.nodeEdgelw = 0
        self.edgeType = 'default'
        self.compartmentColor = 'tab:gray'
        self.compartmentEdgeColor = 'k'
        self.compartmentEdgelw = 2
        self.highlight = []
        self.hlNodeColor = 'tab:purple'
        self.hlNodeEdgeColor = 'tab:pink'
        self.drawReactionNode = True
        self.breakBoundary = False
        self.tightLayout = False
        self.analyzeFlux = False
        self.analyzeRates = False
        self.analyzeColorHigh = 'k'
        self.analyzeColorLow = 'k'
        self.analyzeColorMap = 'Reds'
        self.analyzeColorScale = False
        self.drawInlineTimeCourse = False
        self.simulationStartTime = 0
        self.simulationEndTime = 100
        self.numPoints = 100
        self.plotStatistics = False
        self.forceAnalysisAtEndTime = False
        self.plotColorbar = False
        self.inlineTimeCourseSelections = []
        self.customAxis = None
        self.layoutAlgorithm = 'kamada-kawai'
        self.ignoreLayout = False
        self._Var.pos = None
        self._Var.ignoreAnalysisFlags = False


    def _analyze(self):
        # Analyze the reaction rates
        if self.analyzeFlux:
            if self.forceAnalysisAtEndTime:
                self.rrInstance.simulate(self.simulationStartTime, self.simulationEndTime, self.numPoints)
                self._Var.flux = self.rrInstance.getReactionRates()
            else:
                try:
                    self.rrInstance.steadyState()
                    self._Var.flux = self.rrInstance.getReactionRates()
                except:
                    print("No steadyState is found - netplotlib will use the state at t=simTime")
                    self.rrInstance.reset()
                    self.rrInstance.simulate(self.simulationStartTime, self.simulationEndTime, self.numPoints)
                    self._Var.flux = self.rrInstance.getReactionRates()
                
        if self.analyzeRates:
            self.rrInstance.reset()
            self.rrInstance.simulate(self.simulationStartTime, self.simulationEndTime, self.numPoints)
            self._Var.reaction_rate = self.rrInstance.getRatesOfChange()


    def setLayout(self, pos):
        """
        Set custom layout and bypass whe is generated by the layout algorothm
        
        :param pos: Dictionary of all nodes and corresponding coordinates
        :type name: Dict
        """
        
        self._Var.pos = pos
    

    def getLayout(self, returnState=False):
        """
        Return the layout of the model
        
        :param returnState: boolean flag for returning the networkx.Graph object
        :returns pos: Dictionary of all nodes and corresponding coordinates
        """
        
        if self.layoutAlgorithm not in toolbox.getListOfAlgorithms():
            raise Exception("Unsupported layout algorithm: '" + str(self.layoutAlgorithm) + "'")
        
        avoid = ['C', 'CC', 'Ci', 'E1', 'EX', 'Ei', 'FF', 'GF', 'Ge', 'Gt', 'I', 'LC',
                 'LM', 'LT', 'Le', 'Li', 'Lt', 'N', 'Ne', 'O', 'Q', 'QQ', 'RR', 'S',
                 'Si', 'ZZ', 'ff', 'fu', 'im', 'jn', 'li', 'ln', 'oo', 'pi', 're',
                 'rf', 'yn']
        
        numBnd = self.rrInstance.getNumBoundarySpecies()
        numFlt = self.rrInstance.getNumFloatingSpecies()
        
        # prepare symbols for sympy
        boundaryId_sympy = [] 
        floatingId_sympy = []
        
        # Fix issues with reserved characters
        for i in range(numBnd):
            if self._Var.boundaryId[i] in avoid:
                boundaryId_sympy.append('_' + self._Var.boundaryId[i])
            else:
                boundaryId_sympy.append(self._Var.boundaryId[i])
        
        for i in range(numFlt):
            if self._Var.floatingId[i] in avoid:
                floatingId_sympy.append('_' + self._Var.floatingId[i])
            else:
                floatingId_sympy.append(self._Var.floatingId[i])
        
        paramIdsStr = ' '.join(self.rrInstance.getGlobalParameterIds())
        floatingIdsStr = ' '.join(floatingId_sympy)
        boundaryIdsStr = ' '.join(boundaryId_sympy)
        comparmentIdsStr = ' '.join(self.rrInstance.getCompartmentIds())
        
        allIds = paramIdsStr + ' ' + floatingIdsStr + ' ' + boundaryIdsStr + ' ' + comparmentIdsStr
        
        avsym = sympy.symbols(allIds)
        
        # extract reactant, product, modifiers, and kinetic laws
        self._Var.rct = []
        self._Var.prd = []
        self._Var.mod = []
        self._Var.r_type = []
        rr_type = []
        mod_target = []
        kineticLaw = []
        self._Var.mod_type = []
        
        for slr in self._Var.sbmlmodel.getListOfReactions():
            temprct = []
            tempprd = []
            tempmod = []
            
            sbmlreaction = self._Var.sbmlmodel.getReaction(slr.getId())
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
            
            if len(temprct) == 0:
                self._Var.rct.append(['Input'])
            else:
                self._Var.rct.append(sorted(temprct, key=lambda v: (v.upper(), v[0].islower())))
            if len(tempprd) == 0:
                self._Var.prd.append(['Output'])
            else:
                self._Var.prd.append(sorted(tempprd, key=lambda v: (v.upper(), v[0].islower())))
            
            if self._Var.rct[-1] == self._Var.prd[-1]:
                self._Var.mod.append(self._Var.rct[-1])
            else:
                self._Var.mod.append(sorted(tempmod, key=lambda v: (v.upper(), v[0].islower())))
            
            if sbmlreaction.getReversible():
                rr_type.append('reversible')
            else:
                rr_type.append('irreversible')
                
            if kl == None:
                kineticLaw.append(None)
            else:
                # Update kinetic law according to change in species name
                kl_split = re.split('( |\(|\))', kl.getFormula())
                for i in range(len(kl_split)):
                    if kl_split[i] in avoid:
                        kl_split[i] = '_' + kl_split[i]
                
                kineticLaw.append(''.join(kl_split))
            
        nkl = 0
        # use sympy for analyzing modifiers weSmart
        for ml in range(len(self._Var.mod)):
            mod_type_temp = []
            expression = kineticLaw[ml]
            if expression == None:
                for ml_i in range(len(self._Var.mod[ml])):
                    mod_type_temp.append('modifier')
                self._Var.r_type.append(rr_type[nkl])
                nkl += 1
            else:
                n,d = sympy.fraction(expression)
                for ml_i in range(len(self._Var.mod[ml])):
                    if n.has(sympy.symbols(self._Var.mod[ml][ml_i])):
                        mod_type_temp.append('activator')
                    elif d.has(sympy.symbols(self._Var.mod[ml][ml_i])):
                        mod_type_temp.append('inhibitor')
                    else:
                        mod_type_temp.append('modifier')
                n = '(' + str(n) + ')'
            
                # In case all products are in rate law, assume it is a reversible reaction
                if (all(ext in str(n) for ext in [s + '/' for s in self._Var.prd[ml]]) or
                    all(ext in str(n) for ext in [s + ')' for s in self._Var.prd[ml]]) or
                    all(ext in str(n) for ext in [s + '*' for s in self._Var.prd[ml]]) or
                    all(ext in str(n) for ext in [s + ';' for s in self._Var.prd[ml]]) or
                    all(ext in str(n) for ext in [s + '+' for s in self._Var.prd[ml]]) or
                    all(ext in str(n) for ext in [s + '-' for s in self._Var.prd[ml]]) or
                    all(ext in str(n) for ext in [s + ' ' for s in self._Var.prd[ml]])):
                    self._Var.r_type.append('reversible')
                else:
                    self._Var.r_type.append('irreversible')
            self._Var.mod_type.append(mod_type_temp)
        
        for i in range(len(self._Var.mod)):
            if len(self._Var.mod[i]) > 0:
                mod_target.append(np.repeat(self._Var.rid[i], len(self._Var.mod[i])).tolist())
        
        self._Var.mod_flat = [item for sublist in self._Var.mod for item in sublist]
        self._Var.modtype_flat = [item for sublist in self._Var.mod_type for item in sublist]
        self._Var.modtarget_flat = [item for sublist in mod_target for item in sublist]
        
        self._Var.speciesId = list(self._Var.rct + self._Var.prd)
        self._Var.speciesId = [item for sublist in self._Var.speciesId for item in sublist]
        self._Var.speciesId = list(set(self._Var.speciesId))
        
        if self.breakBoundary:
            speciesId_temp = []
            for i in range(len(self._Var.speciesId)):
                if self._Var.speciesId[i] not in self._Var.boundaryId + ['Input', 'Output']:
                    speciesId_temp.append(self._Var.speciesId[i])
                
            self._Var.speciesId = speciesId_temp
            
            boundaryId_temp = []
            bc = 0
            for i in range(len(self._Var.rid)):
                for j in range(len(self._Var.rct[i])):
                    if self._Var.rct[i][j] in self._Var.boundaryId + ['Input', 'Output']:
                        self._Var.rct[i][j] = self._Var.rct[i][j] + '_' + str(bc)
                        self._Var.speciesId.append(self._Var.rct[i][j])
                        boundaryId_temp.append(self._Var.rct[i][j])
                        bc += 1
                for k in range(len(self._Var.prd[i])):
                    if self._Var.prd[i][k] in self._Var.boundaryId + ['Input', 'Output']:
                        self._Var.prd[i][k] = self._Var.prd[i][k] + '_' + str(bc)
                        self._Var.speciesId.append(self._Var.prd[i][k])
                        boundaryId_temp.append(self._Var.prd[i][k])
                        bc += 1
            self._Var.boundaryId = boundaryId_temp

        # initialize directional graph
        self._Var.G = nx.DiGraph()
    
        # add edges
        for i in range(self._Var.sbmlmodel.getNumReactions()):
            for k in range(len(self._Var.rct[i])):
                self._Var.G.add_edges_from([(self._Var.rct[i][k], self._Var.rid[i])])
            
            for j in range(len(self._Var.prd[i])):
                self._Var.G.add_edges_from([(self._Var.rid[i], self._Var.prd[i][j])])
                        
            if len(self._Var.mod[i]) > 0:
                for l in range(len(self._Var.mod[i])):
                    self._Var.G.add_edges_from([(self._Var.mod[i][l], self._Var.rid[i])])
                
        if ((self._Var.layoutAvail) and not (self.ignoreLayout)):
            pos = layout.getPos(self._Var.layoutPlugin)
        else:
            # calcutate positions
            thres = 0.3
            if self.layoutAlgorithm == 'spring':
                pos = nx.spring_layout(self._Var.G, scale=self.scale, seed=1)
            elif self.layoutAlgorithm == 'kamada-kawai':
                shortest_dist = dict(nx.shortest_path_length(self._Var.G, weight='weight'))
                pos = nx.kamada_kawai_layout(self._Var.G, dist=shortest_dist, scale=self.scale)
                
                maxIter = 5
                maxIter_n = 0
                dist_flag = True
                
                if self.tightLayout:
                    comId = self._Var.speciesId
                else:
                    comId = self._Var.speciesId + self._Var.rid
                
                while dist_flag and (maxIter_n < maxIter):
                    dist_flag = False
                    for i in itertools.combinations(comId, 2):
                        pos_dist = np.sqrt((pos[i[0]][0] - pos[i[1]][0])**2 + (pos[i[0]][1] - pos[i[1]][1])**2)
                        if pos_dist < thres:
                            dist_flag = True
                            shortest_dist[i[0]][i[1]] = 2
                            shortest_dist[i[1]][i[0]] = 2
                    pos = nx.kamada_kawai_layout(self._Var.G, dist=shortest_dist, scale=self.scale)
                    maxIter_n += 1
            elif self.layoutAlgorithm == 'twopi' or self.layoutAlgorithm == 'neato' or self.layoutAlgorithm == 'dot':
                from networkx.drawing.nx_pydot import graphviz_layout
                try:
                    pos = graphviz_layout(self._Var.G, prog=self.layoutAlgorithm)
                except:
                    raise Exception("Error running graphviz: Please check graphviz is properly configured.")
                keylist = np.array(list(pos.keys()))
                poslist = np.array(list(pos.values()))
                poslist /= np.max(np.abs(poslist),axis=0)
                pos = dict(zip(keylist, poslist))
            else:
                raise Exception("Unsupported layout algorithm.")
        
        if returnState:
            return pos, self._Var
        else:
            return pos
    
    
    def generateTimelapse(self, start, end, points, visualize='flux', backend='PIL',
                          savePath=None, plotTime=False, dpi=150, duration=5):
        """
        Generate a timelapse of the specified variable
        
        :param start: start time
        :param end: end time
        :param points: number of timesteps
        :param visualize: variable to visualize. Either 'flux' or 'rate'
        :param backend: backend for file generation. Supports 'PIL' (.gif) or 'cv2' (.avi)
        :param savePath: path to save the diagram
        :param dpi: dpi settings for the diagram
        :param duration: timelapse duration
        """
        
        import tempfile
        if backend=='cv2':
            try:
                import cv2
            except:
                raise Exception("Cannot import opencv")
        elif backend=='PIL':
            try:
                from PIL import Image, ImageFont, ImageDraw
            except:
                raise Exception("Cannot import Pillow")
        else:
            raise Exception("Unsupported backend")
        
        if savePath == None:
            raise Exception("Specify the save path!")
        
        if visualize == 'flux':
            optbackup = self.analyzeFlux
            self.analyzeFlux = True
        elif visualize == 'rate':
            optbackup = self.analyzeRates
            self.analyzeRates = True
        
        self._Var.ignoreAnalysisFlags = True
        
        stbackup = self.simulationStartTime
        etbackup = self.simulationEndTime
        npbackup = self.numPoints
        aetbackup = self.forceAnalysisAtEndTime
        acsbackup = self.analyzeColorScale
        
        self.forceAnalysisAtEndTime = True
        self.analyzeColorScale = True
        pos = self.getLayout()
        self._Var.pos = pos
        
        tempdir = tempfile.TemporaryDirectory()
        
        self.rrInstance.reset()
        self._Var.flux = self.rrInstance.getReactionRates()
        self.draw(show=False, savePath=os.path.join(tempdir.name, "0.png"), dpi=dpi)
        
        for i in range(points-1):
            self.rrInstance.reset()
            self.rrInstance.simulate(start, start+(i+1)*(end-start)/(points-1), i+2)
            self._Var.flux = self.rrInstance.getReactionRates()
            self.draw(show=False, savePath=os.path.join(tempdir.name, "%s.png" % str(i+1)), dpi=dpi)
        
        fl = sorted(os.listdir(tempdir.name), key=lambda x: int(x.split(".")[0]))
        flf = [os.path.join(tempdir.name, f) for f in fl]
        
        if plotTime:
            td = (end-start)/points
            font = ImageFont.truetype("arial.ttf", 35) #TODO: Check if cross-platform compatible
            for i,f in enumerate(flf):
                img = Image.open(f)
                draw = ImageDraw.Draw(img)
                draw.text((10,10), 't = ' + str(i*td), (0,0,0), font=font)
                img.save(f)
        
        if backend == 'PIL':
            img, *imgs = [Image.open(f) for f in flf]
            img.save(fp=savePath, format='GIF', append_images=imgs, save_all=True, 
                     duration=duration, loop=0)
        elif backend == 'cv2':
            frame = cv2.imread(os.path.join(tempdir.name, flf[0]))
            height, width, layers = frame.shape
            
            video = cv2.VideoWriter(savePath, 0, int(len(fl)/duration), (width,height))
            
            for f in flf:
                video.write(cv2.imread(os.path.join(tempdir.name, f)))
            
            cv2.destroyAllWindows()
            video.release()
        
        tempdir.cleanup()
        
        if visualize == 'flux':
            self.analyzeFlux = optbackup
        elif visualize == 'rate':
            self.analyzeRates = optbackup
    
        self.simulationStartTime = stbackup
        self.simulationEndTime = etbackup
        self.numPoints = npbackup
        self.forceAnalysisAtEndTime = aetbackup
        self.analyzeColorScale = acsbackup
        self._Var.ignoreAnalysisFlags = False
    
    
    def draw(self, show=True, savePath=None, dpi=150):
        """
        Draw network diagram
        
        :param show: flag to show the diagram
        :param savePath: path to save the diagram
        :param dpi: dpi settings for the diagram
        """
        
        toolbox.checkValidity(self)
        
        if self._Var.pos == None:
            pos = self.getLayout()
        else:
            pos = self._Var.pos
            assert(len(self._Var.boundaryId)+len(self._Var.floatingId)+len(self._Var.rid) == len(pos))
        
        if (self.analyzeFlux or self.analyzeRates) and not self._Var.ignoreAnalysisFlags:
            self._analyze()
        
        if ((self._Var.layoutAvail) and not (self.ignoreLayout)):
            layout.draw(self, show=show, savePath=savePath, dpi=dpi)
        else:
            # initialize figure
            fig = plt.figure()
    
            if self.drawInlineTimeCourse:
                self.rrInstance.reset()
                if len(self.inlineTimeCourseSelections) == 0:
                    result = self.rrInstance.simulate(self.simulationStartTime, 
                                                      self.simulationEndTime, 
                                                      self.numPoints)
                else:
                    sel = self.inlineTimeCourseSelections
                    if 'time' not in sel:
                        sel = ['time'] + sel
                    result = self.rrInstance.simulate(self.simulationStartTime, 
                                                      self.simulationEndTime,
                                                      self.numPoints, 
                                                      selections=sel)
                
                plt.tight_layout()
                
                colorDict = dict(zip(self.rrInstance.getFloatingSpeciesIds(), 
                                     plt.rcParams['axes.prop_cycle'].by_key()['color'][:self.rrInstance.getNumFloatingSpecies()]))
                
                gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
                gs.update(wspace=0.025, hspace=0.05)
                ax = plt.subplot(gs[1])
                
                if len(self.inlineTimeCourseSelections) == 0:
                    ax.plot(result[:,0], result[:,1], lw=3)
                else:
                    for i in range(len(sel) - 1):
                        ax.plot(result[:,0], result[:,i+1], lw=3, c=colorDict.get(sel[i+1]))
                    
                ax = plt.subplot(gs[0])
            else:
                if self.customAxis != None:
                    ax = self.customAxis
                else:
                    ax = plt.gca()
        
            # check the range of x and y positions
            max_width = []
            max_height = []
            for key, value in pos.items():
                max_width.append(value[0])
                max_height.append(value[1])
            
            max_width = [min(max_width), max(max_width)]
            max_height = [min(max_height), max(max_height)]
            
            # add compartments first
            for com in self._Var.compartmentId:
                comBox = FancyBboxPatch((max_width[0]-(max_width[1]-max_width[0])/10,
                                         max_height[0]-(max_height[1]-max_height[0])/10), 
                                    max_width[1]-max_width[0]+2*(max_width[1]-max_width[0])/10, 
                                    max_height[1]-max_height[0]+2*(max_height[1]-max_height[0])/10,
                                    boxstyle="round,pad=0.01, rounding_size=0.01",
                                    linewidth=self.compartmentEdgelw, 
                                    edgecolor=self.compartmentEdgeColor, 
                                    facecolor=self.compartmentColor,
                                    alpha=0.5)
                ax.add_patch(comBox)
            
            hlInd = 0
            # add nodes to the figure
            for n in self._Var.G:
                if n in self._Var.rid:
                    rec_width = 0.05*(self.fontsize/20)
                    rec_height = 0.05*(self.fontsize/20)
                    if n in self.highlight:
                        if type(self.hlNodeEdgeColor) == list:
                            currHlNodeEdgeColor = self.hlNodeEdgeColor[hlInd]
                            hlInd += 1
                        else:
                            currHlNodeEdgeColor = self.hlNodeEdgeColor
                        if type(self.hlNodeColor) == list:
                            currHlNodeColor = self.hlNodeColor[hlInd]
                            hlInd += 1
                        else:
                            currHlNodeColor = self.hlNodeColor
                            
                        c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                            pos[n][1]-rec_height/2),
                                            rec_width, 
                                            rec_height,
                                            boxstyle="round,pad=0.01, rounding_size=0.01",
                                            linewidth=self.nodeEdgelw, 
                                            edgecolor=currHlNodeEdgeColor, 
                                            facecolor=currHlNodeColor)
                    else:
                        c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                            pos[n][1]-rec_height/2), 
                                            rec_width, 
                                            rec_height,
                                            boxstyle="round,pad=0.01, rounding_size=0.01",
                                            linewidth=self.nodeEdgelw, 
                                            edgecolor=self.nodeEdgeColor, 
                                            facecolor=self.reactionNodeColor)
                    if self.labelReactionIds:
                        ax.text(pos[n][0], 
                                pos[n][1], 
                                n, 
                                fontsize=self.fontsize, 
                                horizontalalignment='center', 
                                verticalalignment='center', 
                                color=self.labelColor)
                else:
                    if len(n) > 10:
                        rec_width = max(0.045*((len(n)/2)+1), 0.13)*(self.fontsize/20)
                        rec_height = 0.18*(self.fontsize/20)
                    else:
                        rec_width = max(0.04*(len(n)+1), 0.12)*(self.fontsize/20)
                        rec_height = 0.07*(self.fontsize/20)
                        
                    if self.drawInlineTimeCourse:
                        node_color = colorDict[n]
                    else:
                        if (n in self._Var.boundaryId) or (n == 'Input') or (n == 'Output'):
                            node_color = self.boundaryColor
                        else:
                            if self.analyzeRates:
                                norm = colors.Normalize(vmin=np.min(self._Var.reaction_rate),vmax=np.max(self._Var.reaction_rate))
                                colormap = cm.get_cmap(self.analyzeColorMap)
                                node_color = colormap(norm(self._Var.reaction_rate[0][self._Var.reaction_rate.colnames.index(n)]))
                            else:
                                node_color = self.nodeColor
                        
                    if n in self.highlight:
                        if type(self.hlNodeEdgeColor) == list:
                            currHlNodeEdgeColor = self.hlNodeEdgeColor[hlInd]
                            hlInd += 1
                        else:
                            currHlNodeEdgeColor = self.hlNodeEdgeColor
                        if type(self.hlNodeColor) == list:
                            currHlNodeColor = self.hlNodeColor[hlInd]
                            hlInd += 1
                        else:
                            currHlNodeColor = self.hlNodeColor
                            
                        c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                            pos[n][1]-rec_height/2),
                                            rec_width, 
                                            rec_height,
                                            boxstyle="round,pad=0.01, rounding_size=0.02",
                                            linewidth=self.nodeEdgelw, 
                                            edgecolor=currHlNodeEdgeColor, 
                                            facecolor=currHlNodeColor)
                    else:
                        c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                            pos[n][1]-rec_height/2), 
                                            rec_width, 
                                            rec_height,
                                            boxstyle="round,pad=0.01, rounding_size=0.02",
                                            linewidth=self.nodeEdgelw, 
                                            edgecolor=self.nodeEdgeColor, 
                                            facecolor=node_color)
                    if len(n) > 10:
                        ax.text(pos[n][0], pos[n][1], n[:int(len(n)/2)] + '\n' + n[int(len(n)/2):], 
                                fontsize=self.fontsize, horizontalalignment='center', 
                                verticalalignment='center', color=self.labelColor)
                    else:
                        ax.text(pos[n][0], pos[n][1], n, 
                                fontsize=self.fontsize, horizontalalignment='center', 
                                verticalalignment='center', color=self.labelColor)
                self._Var.G.nodes[n]['patch'] = c
            
            # add edges to the figure
            for i in range(len(self._Var.rid)):
                if (len(self._Var.rct[i]) == 1) or (len(self._Var.prd[i]) == 1): # UNI-involved
                    comb = list(itertools.combinations_with_replacement(self._Var.rct[i],len(self._Var.prd[i])))
                    for j in [list(zip(x,self._Var.prd[i])) for x in comb]:
                        for k in range(len(j)):
                            p1 = self._Var.G.nodes[j[k][0]]['patch']
                            p2 = self._Var.G.nodes[self._Var.rid[i]]['patch']
                            p3 = self._Var.G.nodes[j[k][1]]['patch']
                            
                            X1 = (p1.get_x()+p1.get_width()/2,p1.get_y()+p1.get_height()/2)
                            X2 = (p2.get_x()+p2.get_width()/2,p2.get_y()+p2.get_height()/2)
                            X3 = (p3.get_x()+p3.get_width()/2,p3.get_y()+p3.get_height()/2)
                            
                            if ((len(np.unique(self._Var.rct[i])) > len(self._Var.prd[i])) or 
                                (len(self._Var.rct[i]) < len(np.unique(self._Var.prd[i])))): # Uni-Bi or Bi-Uni
                                XY1 = np.vstack((X1, X2))
                                XY2 = np.vstack((X2, X3))
                                
                                tck1, u1 = interpolate.splprep([XY1[:,0], XY1[:,1]], 
                                                               k=1)
                                intX1, intY1 = interpolate.splev(np.linspace(0, 1, 100),
                                                                 tck1, 
                                                                 der=0)
                                stackXY1 = np.vstack((intX1, intY1))
                                tck2, u2 = interpolate.splprep([XY2[:,0], XY2[:,1]], 
                                                               k=1)
                                intX2, intY2 = interpolate.splev(np.linspace(0, 1, 100), 
                                                                 tck2, 
                                                                 der=0)
                                stackXY2 = np.vstack((intX2, intY2))
                                
                                if max(stackXY1[0]) > max_width[1]:
                                    max_width[1] = max(stackXY1[0])
                                if min(stackXY1[0]) < max_width[0]:
                                    max_width[0] = min(stackXY1[0])
                                if max(stackXY1[1]) > max_height[1]:
                                    max_height[1] = max(stackXY1[1])
                                if min(stackXY1[1]) < max_height[0]:
                                    max_height[0] = min(stackXY1[1])
                                    
                                if max(stackXY2[0]) > max_width[1]:
                                    max_width[1] = max(stackXY2[0])
                                if min(stackXY2[0]) < max_width[0]:
                                    max_width[0] = min(stackXY2[0])
                                if max(stackXY2[1]) > max_height[1]:
                                    max_height[1] = max(stackXY2[1])
                                if min(stackXY2[1]) < max_height[0]:
                                    max_height[0] = min(stackXY2[1])
                                
                                X3top = (p3.get_x()+p3.get_width()/2,
                                         p3.get_y()+p3.get_height())
                                X3bot = (p3.get_x()+p3.get_width()/2,
                                         p3.get_y())
                                X3left = (p3.get_x(),
                                          p3.get_y()+p3.get_height()/2)
                                X3right = (p3.get_x()+p3.get_width(),
                                           p3.get_y()+p3.get_height()/2)
                                
                                n_1 = -1
                                arrthres_h = .02
                                arrthres_v = .02
                                while (((stackXY2.T[n_1][0] > (X3left[0]-arrthres_h)) and
                                        (stackXY2.T[n_1][0] < (X3right[0]+arrthres_h)) and
                                        (stackXY2.T[n_1][1] > (X3bot[1]-arrthres_v)) and 
                                        (stackXY2.T[n_1][1] < (X3top[1]+arrthres_v))) and
                                        (np.abs(n_1) < np.shape(stackXY2)[1] - 75)):
                                    n_1 -= 1
                                
                                if self.edgeType == 'default':
                                    lpath1 = Path(stackXY1.T)
                                    lpath2 = Path(stackXY2.T[:n_1])
                                elif self.edgeType == 'bezier':
                                    bcp = toolbox.computeBezierControlPoints(stackXY1.T[0], X2, stackXY1.T[-1])
                                    lpath1 = Path([bcp[0], bcp[1], bcp[2], bcp[3]], 
                                                 [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                                    bcp = toolbox.computeBezierControlPoints(stackXY2.T[0], X2, stackXY2.T[n_1])
                                    lpath2 = Path([bcp[0], bcp[1], bcp[2], bcp[3]], 
                                                 [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                                
                                lw1 = (1+self.edgelw)
                                lw2 = (1+self.edgelw)
                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                e1color = self.reactionColor
                                e2color = self.reactionColor
                                
                                if self._Var.r_type[i] == 'reversible':
                                    X1top = (p1.get_x()+p1.get_width()/2,
                                             p1.get_y()+p1.get_height())
                                    X1bot = (p1.get_x()+p1.get_width()/2,
                                             p1.get_y())
                                    X1left = (p1.get_x(),
                                              p1.get_y()+p1.get_height()/2)
                                    X1right = (p1.get_x()+p1.get_width(),
                                               p1.get_y()+p1.get_height()/2)
                                    
                                    n_2 = 0
                                    while (((stackXY1.T[n_2][0] > (X1left[0]-arrthres_h)) and 
                                            (stackXY1.T[n_2][0] < (X1right[0]+arrthres_h)) and
                                            (stackXY1.T[n_2][1] > (X1bot[1]-arrthres_v)) and 
                                            (stackXY1.T[n_2][1] < (X1top[1]+arrthres_v))) and
                                            (np.abs(n_2) < np.shape(stackXY1)[1] - 75)):
                                        n_2 += 1
                                        
                                    if self.edgeType == 'default':
                                        lpath1 = Path(stackXY1.T[n_2:])
                                    elif self.edgeType == 'bezier':
                                        bcp = toolbox.computeBezierControlPoints(stackXY1.T[n_2], X2, stackXY1.T[-1])
                                        lpath1 = Path([bcp[0], bcp[1], bcp[2], bcp[3]], 
                                                     [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                                else:
                                    arrowstyle1 = ArrowStyle.Curve()
    
                                if self.analyzeFlux:
                                    if self._Var.flux[i] > 0:
                                        lw1 = (1+self.edgelw)
                                        lw2 = (4+self.edgelw)
                                        arrowstyle2 = ArrowStyle.CurveFilledB(head_length=1.2, head_width=0.8)
                                    elif self._Var.flux[i] < 0:
                                        lw1 = (4+self.edgelw)
                                        lw2 = (1+self.edgelw)
                                        arrowstyle1 = ArrowStyle.CurveFilledA(head_length=1.2, head_width=0.8)
                                    
                                    if self.analyzeColorScale:
                                        norm = colors.Normalize(vmin=np.min(self._Var.flux),vmax=np.max(self._Var.flux))
                                        colormap = cm.get_cmap(self.analyzeColorMap)
                                        e1color = colormap(norm(self._Var.flux[i]))#0.5-Var.flux[i]/(2*max(abs(self._Var.flux))))
                                        e2color = colormap(norm(self._Var.flux[i]))#0.5+Var.flux[i]/(2*max(abs(self._Var.flux))))
                                    else:
                                        e1color = self.analyzeColorLow
                                        e2color = self.analyzeColorHigh
                                
                                e1 = FancyArrowPatch(path=lpath1,
                                                    arrowstyle=arrowstyle1,
                                                    mutation_scale=10.0,
                                                    lw=lw1,
                                                    color=e1color)
                                
                                e2 = FancyArrowPatch(path=lpath2,
                                                    arrowstyle=arrowstyle2,
                                                    mutation_scale=10.0,
                                                    lw=lw2,
                                                    color=e2color)
                                    
                                ax.add_patch(e1)
                                ax.add_patch(e2)
                                
                                if j[k][0] in self._Var.floatingId:
                                    if (np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][0])][i]) > 1):
                                        # position calculation
                                        slope = ((lpath1.vertices[0][1] - lpath1.vertices[int(0.35*len(lpath1.vertices))][1])/
                                                 (lpath1.vertices[0][0] - lpath1.vertices[int(0.35*len(lpath1.vertices))][0]))
                                        x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                        y_prime = -slope*x_prime
                                        ax.text(x_prime+lpath1.vertices[int(0.35*len(lpath1.vertices))][0], 
                                                y_prime+lpath1.vertices[int(0.35*len(lpath1.vertices))][1], 
                                                int(np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][0])][i])), 
                                                fontsize=self.fontsize, 
                                                horizontalalignment='center', 
                                                verticalalignment='center', 
                                                color=self.reactionColor)
                                
                                if j[k][1] in self._Var.floatingId:
                                    if (np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][1])][i]) > 1):
                                        slope = ((lpath2.vertices[0][1] - lpath2.vertices[int(0.65*len(lpath2.vertices))][1])/
                                                 (lpath2.vertices[0][0] - lpath2.vertices[int(0.65*len(lpath2.vertices))][0]))
                                        x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                        y_prime = -slope*x_prime
                                        ax.text(x_prime+lpath2.vertices[int(0.65*len(lpath2.vertices))][0], 
                                                y_prime+lpath2.vertices[int(0.65*len(lpath2.vertices))][1], 
                                                int(np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][1])][i])), 
                                                fontsize=self.fontsize, 
                                                horizontalalignment='center', 
                                                verticalalignment='center', 
                                                color=self.reactionColor)
                                
                            else: 
                                if ((self._Var.rct[i] == self._Var.prd[i]) and (self._Var.rct[i] == self._Var.mod[i])): # Autoregulation
                                    lw1 = (1+self.edgelw)
                                    
                                    center = [(X1[0]+X2[0])/2, (X1[1]+X2[1])/2]
                                    radius = np.sqrt(np.square(X1[0]-X2[0])+np.square(X1[1]-X2[1]))
                                    angle = np.arctan((X1[1]-X2[1])/(X1[0]-X2[0]))/np.pi*180
                                    
                                    if (center[0]+radius)/2 > max_width[1]:
                                        max_width[1] = (center[0]+radius)/2
                                    if (center[0]-radius)/2 < max_width[0]:
                                        max_width[0] = (center[0]-radius)/2
                                    if (center[1]+radius)/2 > max_height[1]:
                                        max_height[1] = (center[1]+radius)/2
                                    if (center[1]-radius)/2 < max_height[0]:
                                        max_height[0] = (center[1]-radius)/2
                                    
                                    if self._Var.mod_type[i][0] == 'modifier':
                                        r = Arc(center, 
                                                radius,
                                                radius,
                                                theta1=angle+10, 
                                                theta2=360+angle-10,
                                                lw=lw1,
                                                linestyle=':',
                                                color=self.reactionColor)
                                    else:
                                        r = Arc(center, 
                                                radius,
                                                radius,
                                                theta1=angle+10, 
                                                theta2=360+angle-10,
                                                lw=lw1,
                                                color=self.reactionColor)
                                    ax.add_patch(r)
                                    
                                    endX1 = center[0]+(radius/2)*np.cos(np.radians(360+angle-10))
                                    endY1 = center[1]+(radius/2)*np.sin(np.radians(360+angle-10))
                                    
                                    if self._Var.mod_type[i][0] == 'inhibitor':
                                        slope = (endY1-center[1])/(endX1-center[0])
                                        endX2 = center[0]+(radius/2-0.07/2)*np.cos(np.radians(360+angle-10))
                                        endY2 = center[1]+(radius/2-0.07/2)*np.sin(np.radians(360+angle-10))
                                        
                                        ar = Rectangle((endX2, 
                                                        endY2),
                                                        0.07, 
                                                        0.001*lw1,
                                                        angle=360+angle-10,
                                                        linewidth=lw1,
                                                        color=self.reactionColor)
                                        ax.add_patch(ar)
                                    else:
                                        ar = RegularPolygon((endX1, endY1), 
                                                            3,
                                                            lw1/120,
                                                            np.radians(360+angle-10),
                                                            color=self.reactionColor)
                                    
                                        ax.add_patch(ar)
                                    
                                else: # Uni-Uni
                                    XY = np.vstack((X1, X2, X3))
                                    
                                    tck, u = interpolate.splprep([XY[:,0], XY[:,1]], k=2)
                                    intX, intY = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
                                    stackXY = np.vstack((intX, intY))
                                    
                                    if max(stackXY[0]) > max_width[1]:
                                        max_width[1] = max(stackXY[0])
                                    if min(stackXY[0]) < max_width[0]:
                                        max_width[0] = min(stackXY[0])
                                    if max(stackXY[1]) > max_height[1]:
                                        max_height[1] = max(stackXY[1])
                                    if min(stackXY[1]) < max_height[0]:
                                        max_height[0] = min(stackXY[1])
                                    
                                    X3top = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y()+p3.get_height())
                                    X3bot = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y())
                                    X3left = (p3.get_x(),
                                              p3.get_y()+p3.get_height()/2)
                                    X3right = (p3.get_x()+p3.get_width(),
                                               p3.get_y()+p3.get_height()/2)
                                    
                                    n_1 = -1
                                    arrthres_h = .02
                                    arrthres_v = .02
                                    while (((stackXY.T[n_1][0] > (X3left[0]-arrthres_h)) and 
                                            (stackXY.T[n_1][0] < (X3right[0]+arrthres_h)) and
                                            (stackXY.T[n_1][1] > (X3bot[1]-arrthres_v)) and 
                                            (stackXY.T[n_1][1] < (X3top[1]+arrthres_v))) and
                                            (np.abs(n_1) < np.shape(stackXY)[1] - 75)):
                                        n_1 -= 1
                                   
                                    if self._Var.r_type[i] == 'reversible':
                                        X1top = (p1.get_x()+p1.get_width()/2,
                                                 p1.get_y()+p1.get_height())
                                        X1bot = (p1.get_x()+p1.get_width()/2,
                                                 p1.get_y())
                                        X1left = (p1.get_x(),
                                                  p1.get_y()+p1.get_height()/2)
                                        X1right = (p1.get_x()+p1.get_width(),
                                                   p1.get_y()+p1.get_height()/2)
                                        
                                        n_2 = 0
                                        while (((stackXY.T[n_2][0] > (X1left[0]-arrthres_h)) and 
                                                (stackXY.T[n_2][0] < (X1right[0]+arrthres_h)) and
                                                (stackXY.T[n_2][1] > (X1bot[1]-arrthres_v)) and 
                                                (stackXY.T[n_2][1] < (X1top[1]+arrthres_v))) and
                                                (np.abs(n_2) < np.shape(stackXY)[1] - 75)):
                                            n_2 += 1
                                        
                                        if self.edgeType == 'default':
                                            lpath = Path(stackXY.T[n_2:n_1])
                                        elif self.edgeType == 'bezier':
                                            bcp = toolbox.computeBezierControlPoints(stackXY.T[n_2], X2, stackXY.T[n_1])
                                            lpath = Path([bcp[0], bcp[1], bcp[2], bcp[3]], 
                                                         [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                                        
                                        if self.analyzeFlux:
                                            if self._Var.flux[i] > 0:
                                                lw1 = (1+self.edgelw)
                                                lw2 = (4+self.edgelw)
                                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=1.2, head_width=0.8)
                                            elif self._Var.flux[i] < 0:
                                                lw1 = (4+self.edgelw)
                                                lw2 = (1+self.edgelw)
                                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=1.2, head_width=0.8)
                                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                            else:
                                                lw1 = (1+self.edgelw)
                                                lw2 = (1+self.edgelw)
                                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                            
                                            if self.analyzeColorScale:
                                                norm = colors.Normalize(vmin=np.min(self._Var.flux),vmax=np.max(self._Var.flux))
                                                colormap = cm.get_cmap(self.analyzeColorMap)
                                                e1color = colormap(norm(self._Var.flux[i]))
                                                e2color = colormap(norm(self._Var.flux[i]))
                                            else:
                                                e1color = self.analyzeColorLow
                                                e2color = self.analyzeColorHigh
                                                
                                            e1 = FancyArrowPatch(path=Path(stackXY.T[-n_1:50]),
                                                                arrowstyle=arrowstyle1,
                                                                mutation_scale=10.0,
                                                                lw=lw1,
                                                                color=e1color)
                                            e2 = FancyArrowPatch(path=Path(stackXY.T[50:n_1]),
                                                                arrowstyle=arrowstyle2,
                                                                mutation_scale=10.0,
                                                                lw=lw2,
                                                                color=e2color)
                                            ax.add_patch(e1)
                                            ax.add_patch(e2)
                                        else:
                                            arrowstyle = ArrowStyle.CurveFilledAB(head_length=0.8, head_width=0.4)
                                            e = FancyArrowPatch(path=lpath,
                                                            arrowstyle=arrowstyle,
                                                            mutation_scale=10.0,
                                                            lw=(1+self.edgelw),
                                                            color=self.reactionColor)
                                            ax.add_patch(e)
                                        
                                    else:
                                        e1color = self.reactionColor
                                        if self.analyzeFlux:
                                            if self.analyzeColorScale:
                                                norm = colors.Normalize(vmin=np.min(self._Var.flux),vmax=np.max(self._Var.flux))
                                                colormap = cm.get_cmap(self.analyzeColorMap)
                                                e1color = colormap(norm(self._Var.flux[i]))
                                            else:
                                                e1color = self.reactionColor
                                        
                                        if self.edgeType == 'default':
                                            lpath = Path(stackXY.T[:n_1])
                                        elif self.edgeType == 'bezier':
                                            bcp = toolbox.computeBezierControlPoints(stackXY.T[0], X2, stackXY.T[n_1])
                                            lpath = Path([bcp[0], bcp[1], bcp[2], bcp[3]], 
                                                         [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                                        
                                        arrowstyle1 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                        lw1 = (1+self.edgelw)
                                        e = FancyArrowPatch(path=lpath,
                                                            arrowstyle=arrowstyle1,
                                                            mutation_scale=10.0,
                                                            lw=lw1,
                                                            color=e1color)
                                        ax.add_patch(e)
                                
                                    if j[k][0] in self._Var.floatingId:
                                        if (np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][0])][i]) > 1):
                                            slope = ((lpath.vertices[0][1] - lpath.vertices[int(0.35*len(lpath.vertices))][1])/
                                                     (lpath.vertices[0][0] - lpath.vertices[int(0.35*len(lpath.vertices))][0]))
                                            x_prime = np.sqrt(0.01/(1 + np.square(slope)))*max(self.scale/2, 1)
                                            y_prime = -slope*x_prime
                                            ax.text(x_prime+lpath.vertices[int(0.35*len(lpath.vertices))][0], 
                                                    y_prime+lpath.vertices[int(0.35*len(lpath.vertices))][1], 
                                                    int(np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][0])][i])), 
                                                    fontsize=self.fontsize, 
                                                    horizontalalignment='center', 
                                                    verticalalignment='center', 
                                                    color=self.reactionColor)
                                    
                                    if j[k][1] in self._Var.floatingId:
                                        if (np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][1])][i]) > 1):
                                            slope = ((lpath.vertices[0][1] - lpath.vertices[int(0.75*len(lpath.vertices))][1])/
                                                     (lpath.vertices[0][0] - lpath.vertices[int(0.75*len(lpath.vertices))][0]))
                                            x_prime = np.sqrt(0.01/(1 + np.square(slope)))*max(self.scale/2, 1)
                                            y_prime = -slope*x_prime
                                            ax.text(x_prime+lpath.vertices[int(0.75*len(lpath.vertices))][0], 
                                                    y_prime+lpath.vertices[int(0.75*len(lpath.vertices))][1],
                                                    int(np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][1])][i])), 
                                                    fontsize=self.fontsize, 
                                                    horizontalalignment='center', 
                                                    verticalalignment='center',
                                                    color=self.reactionColor)
                        
                else: # BIBI or larger
                    if len(self._Var.rct[i]) < len(self._Var.prd[i]):
                        rVal = len(self._Var.rct[i])
                    else:
                        rVal = len(self._Var.prd[i])
                        
                    for j in [list(zip(x,self._Var.prd[i])) for x in itertools.combinations(self._Var.rct[i],rVal)][0]:
                        p1 = self._Var.G.nodes[j[0]]['patch']
                        p2 = self._Var.G.nodes[self._Var.rid[i]]['patch']
                        p3 = self._Var.G.nodes[j[1]]['patch']
                        
                        X1 = (p1.get_x()+p1.get_width()/2,p1.get_y()+p1.get_height()/2)
                        X2 = (p2.get_x()+p2.get_width()/2,p2.get_y()+p2.get_height()/2)
                        X3 = (p3.get_x()+p3.get_width()/2,p3.get_y()+p3.get_height()/2)
                        
                        XY = np.vstack((X1, X2, X3))
                        
                        tck, u = interpolate.splprep([XY[:,0], XY[:,1]], k=2)
                        intX, intY = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
                        stackXY = np.vstack((intX, intY))
                        
                        if max(stackXY[0]) > max_width[1]:
                            max_width[1] = max(stackXY[0])
                        if min(stackXY[0]) < max_width[0]:
                            max_width[0] = min(stackXY[0])
                        if max(stackXY[1]) > max_height[1]:
                            max_height[1] = max(stackXY[1])
                        if min(stackXY[1]) < max_height[0]:
                            max_height[0] = min(stackXY[1])
                        
                        X3top = (p3.get_x()+p3.get_width()/2,
                                 p3.get_y()+p3.get_height())
                        X3bot = (p3.get_x()+p3.get_width()/2,
                                 p3.get_y())
                        X3left = (p3.get_x(),
                                  p3.get_y()+p3.get_height()/2)
                        X3right = (p3.get_x()+p3.get_width(),
                                   p3.get_y()+p3.get_height()/2)
                        
                        n_1 = -1
                        arrthres_h = .02
                        arrthres_v = .02
                        while (((stackXY.T[n_1][0] > (X3left[0]-arrthres_h)) and
                                (stackXY.T[n_1][0] < (X3right[0]+arrthres_h)) and
                                (stackXY.T[n_1][1] > (X3bot[1]-arrthres_v)) and 
                                (stackXY.T[n_1][1] < (X3top[1]+arrthres_v))) and
                                (np.abs(n_1) < np.shape(stackXY)[1] - 75)):
                            n_1 -= 1
                        
                        if self._Var.r_type[i] == 'reversible':
                            X1top = (p1.get_x()+p1.get_width()/2,
                                     p1.get_y()+p1.get_height())
                            X1bot = (p1.get_x()+p1.get_width()/2,
                                     p1.get_y())
                            X1left = (p1.get_x(),
                                      p1.get_y()+p1.get_height()/2)
                            X1right = (p1.get_x()+p1.get_width(),
                                       p1.get_y()+p1.get_height()/2)
                            
                            n_2 = 0
                            while (((stackXY.T[n_2][0] > (X1left[0]-arrthres_h)) and 
                                    (stackXY.T[n_2][0] < (X1right[0]+arrthres_h)) and
                                    (stackXY.T[n_2][1] > (X1bot[1]-arrthres_v)) and 
                                    (stackXY.T[n_2][1] < (X1top[1]+arrthres_v))) and
                                    (np.abs(n_2) < np.shape(stackXY)[1] - 75)):
                                n_2 += 1
                            
                            if self.edgeType == 'default':
                                lpath = Path(stackXY.T[n_2:n_1])
                            elif self.edgeType == 'bezier':
                                bcp = toolbox.computeBezierControlPoints(stackXY.T[n_2], X2, stackXY.T[n_1])
                                lpath = Path([bcp[0], bcp[1], bcp[2], bcp[3]], 
                                             [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                            
                            if self.analyzeFlux:
                                if self._Var.flux[i] > 0:
                                    lw1 = (1+self.edgelw)
                                    lw2 = (4+self.edgelw)
                                    arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                    arrowstyle2 = ArrowStyle.CurveFilledB(head_length=1.2, head_width=0.8)
                                elif self._Var.flux[i] < 0:
                                    lw1 = (4+self.edgelw)
                                    lw2 = (1+self.edgelw)
                                    arrowstyle1 = ArrowStyle.CurveFilledA(head_length=1.2, head_width=0.8)
                                    arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                else:
                                    lw1 = (1+self.edgelw)
                                    lw2 = (1+self.edgelw)
                                    arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                    arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                
                                if self.analyzeColorScale:
                                    norm = colors.Normalize(vmin=np.min(self._Var.flux),vmax=np.max(self._Var.flux))
                                    colormap = cm.get_cmap(self.analyzeColorMap)
                                    e1color = colormap(norm(self._Var.flux[i]))
                                    e2color = colormap(norm(self._Var.flux[i]))
                                else:
                                    e1color = self.analyzeColorLow
                                    e2color = self.analyzeColorHigh
                                    
                                e1 = FancyArrowPatch(path=Path(stackXY.T[n_2:50]),
                                                    arrowstyle=arrowstyle1,
                                                    mutation_scale=10.0,
                                                    lw=lw1,
                                                    color=e1color)
                                e2 = FancyArrowPatch(path=Path(stackXY.T[50:n_1]),
                                                    arrowstyle=arrowstyle2,
                                                    mutation_scale=10.0,
                                                    lw=lw2,
                                                    color=e2color)
                                ax.add_patch(e1)
                                ax.add_patch(e2)
                            else:
                                arrowstyle = ArrowStyle.CurveFilledAB(head_length=0.8, head_width=0.4)
                                e = FancyArrowPatch(path=lpath,
                                                arrowstyle=arrowstyle,
                                                mutation_scale=10.0,
                                                lw=(1+self.edgelw),
                                                color=self.reactionColor)
                                ax.add_patch(e)
                            
                        else:
                            e1color = self.reactionColor
                            if self.analyzeFlux:
                                if self.analyzeColorScale:
                                    norm = colors.Normalize(vmin=np.min(self._Var.flux),vmax=np.max(self._Var.flux))
                                    colormap = cm.get_cmap(self.analyzeColorMap)
                                    e1color = colormap(norm(self._Var.flux[i]))
                                else:
                                    e1color = self.reactionColor
                            
                            if self.edgeType == 'default':
                                lpath = Path(stackXY.T[:n_1])
                            elif self.edgeType == 'bezier':
                                bcp = toolbox.computeBezierControlPoints(stackXY.T[0], X2, stackXY.T[n_1])
                                lpath = Path([bcp[0], bcp[1], bcp[2], bcp[3]], 
                                             [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                            
                            arrowstyle1 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                            lw1 = (1+self.edgelw)
                            e = FancyArrowPatch(path=lpath,
                                                arrowstyle=arrowstyle1,
                                                mutation_scale=10.0,
                                                lw=lw1,
                                                color=e1color)
                            ax.add_patch(e)
                        
                        if j[0] in self._Var.floatingId:
                            if (np.abs(self._Var.stoch[self._Var.stoch_row.index(j[0])][i]) > 1):
                                slope = ((lpath.vertices[0][1] - lpath.vertices[int(0.15*len(lpath.vertices))][1])/
                                         (lpath.vertices[0][0] - lpath.vertices[int(0.15*len(lpath.vertices))][0]))
                                x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                y_prime = -slope*x_prime
                                ax.text(x_prime+lpath.vertices[int(0.15*len(lpath.vertices))][0], 
                                        y_prime+lpath.vertices[int(0.15*len(lpath.vertices))][1], 
                                        int(np.abs(self._Var.stoch[self._Var.stoch_row.index(j[0])][i])), 
                                        fontsize=self.fontsize, 
                                        horizontalalignment='center', 
                                        verticalalignment='center', 
                                        color=self.reactionColor)
                        if j[1] in self._Var.floatingId:
                            if (np.abs(self._Var.stoch[self._Var.stoch_row.index(j[1])][i]) > 1):
                                slope = ((lpath.vertices[0][1] - lpath.vertices[int(0.8*len(lpath.vertices))][1])/
                                         (lpath.vertices[0][0] - lpath.vertices[int(0.8*len(lpath.vertices))][0]))
                                x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                y_prime = -slope*x_prime
                                ax.text(x_prime+lpath.vertices[int(0.8*len(lpath.vertices))][0], 
                                        y_prime+lpath.vertices[int(0.8*len(lpath.vertices))][1], 
                                        int(np.abs(self._Var.stoch[self._Var.stoch_row.index(j[1])][i])), 
                                        fontsize=self.fontsize,
                                        horizontalalignment='center', 
                                        verticalalignment='center', 
                                        color=self.reactionColor)
                        
            # Modifiers
            seen={}
            mc = 0
            for m in range(len(self._Var.mod)):
                if self._Var.rct[m] == self._Var.prd[m]:
                    pass
                else:
                    for j in range(len(self._Var.mod[m])):
                        s = self._Var.mod[m][j]
                        n1 = self._Var.G.nodes[s]['patch']
                        n2 = self._Var.G.nodes[self._Var.modtarget_flat[mc]]['patch']
                        rad = 0.1
                        shrinkB = 2.
                        
                        if (s,self._Var.modtarget_flat[mc]) in seen:
                            rad = seen.get((s,self._Var.modtarget_flat[mc])) # TODO: No curvature when there is just a single line between two nodes
                            rad = (rad+np.sign(rad)*0.1)*-1 # TODO: Change curvature
                            
                        X1 = (n1.get_x()+n1.get_width()/2,
                              n1.get_y()+n1.get_height()/2)
                        X2 = (n2.get_x()+n2.get_width()/2,
                              n2.get_y()+n2.get_height()/2)
                        
                        if self._Var.modtype_flat[mc] == 'inhibitor': # inhibition
                            color = self.modifierColor
                            arrowstyle = ArrowStyle.BarAB(widthA=0.0, angleA=None, widthB=1.0, angleB=None)
                            shrinkB = 10.
                            linestyle = '-'
                        elif self._Var.modtype_flat[mc] == 'activator': # activation
                            color = self.modifierColor
                            arrowstyle = arrowstyle = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                            linestyle = '-'
                        elif self._Var.modtype_flat[mc] == 'modifier': # Unknown modifier
                            color = self.modifierColor
                            arrowstyle = arrowstyle = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                            linestyle = ':'
                        e = FancyArrowPatch(X1,
                                            X2,
                                            patchA=n1,
                                            patchB=n2,
                                            shrinkB=shrinkB,
                                            arrowstyle=arrowstyle,
                                            connectionstyle='arc3,rad=%s'%rad,
                                            mutation_scale=10.0,
                                            lw=(1+self.edgelw),
                                            color=color,
                                            linestyle=linestyle)
                        seen[(s,self._Var.modtarget_flat[mc])]=rad
                        ax.add_patch(e)
                        ax.add_patch(n1)
                        mc += 1
        
            # Add reaction nodes at last to put it on top
            if self.drawReactionNode:
                allnodes = self._Var.speciesId + self._Var.rid
            else:
                allnodes = self._Var.speciesId
            
            if 'Input' in self._Var.G.nodes:
                allnodes += ['Input']
            if 'Output' in self._Var.G.nodes:
                allnodes += ['Output']
            for i in range(len(allnodes)):
                ax.add_patch(self._Var.G.nodes[allnodes[i]]['patch'])
        
            # Statistics
            if self.plotStatistics and (self.analyzeFlux or self.analyzeRates):
                pass
            
            # reset width and height
            ax.autoscale()
            fig.set_figwidth((abs(max_width[0] - max_width[1])+0.5)*5)
            fig.set_figheight((abs(max_height[0] - max_height[1])+0.5)*5)
            if self.plotColorbar:
                fig.set_figwidth((abs(max_width[0] - max_width[1])+0.5)*5 + 4)
                sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
                sm.set_array([])
                plt.colorbar(sm)
            plt.axis('off')
            
            if not self.drawInlineTimeCourse:
                plt.tight_layout()
            
            if savePath != None:
                try:
                    fig.savefig(savePath, bbox_inches='tight', dpi=dpi)
                except IOError as e:
                    raise Exception("Error saving diagram: " + str(e))
                    
            if show and self.customAxis == None:
                plt.show()
            plt.close()


    def savefig(self, path, dpi=150):
        """
        Save network diagram to specified location
        
        :param path: path to save the diagram
        :param dpi: dpi settings for the diagram
        """
        
        self.draw(show=False, savePath=path, dpi=dpi)
        

class NetworkEnsemble():
    
    def __init__(self, models):
        """
        Creates a new NetworkEnsemble object. 
        
        :param models: list of SBML or Antimony strings of models
        :type name: list
        """
        
        self._Var = _Variable()
        self._Var.models = models
        self._Var.boundaryIds = []
        self._Var.floatingIds = []
        self._Var.rids = []
        self.rrInstances = []
        
        for m in models:
            try:
                r = te.loadSBMLModel(m)
                self.rrInstances.append(r)
            except:
                try:
                    r = te.loadAntimonyModel(m)
                    self.rrInstances.append(r)
                except:
                    raise Exception("Input does not seem to be a valid list of SBML or Antimony string")
            
            try:
                self._Var.boundaryIds.append(r.getBoundarySpeciesIds())
                self._Var.floatingIds.append(r.getFloatingSpeciesIds())
                self._Var.rids.append(r.getReactionIds())
            except:
                raise Exception("Failed to analyze the file: Check the file is valid")
        
        self.reset()
        
    
    def reset(self):
        """
        Resets all properties
        """
        
        self.scale = 1.
        self.fontsize = 10
        self.edgelw = 10.
        self.nodeColor = 'tab:blue'
        self.reactionNodeColor = 'tab:gray'
        self.labelColor = 'w'
        self.labelReactionIds = False
        self.reactionColor = 'k'
        self.modifierColor = 'tab:red'
        self.boundaryColor = 'tab:green'
        self.nodeEdgeColor = 'k'
        self.nodeEdgelw = 0
        self.highlight = []
        self.hlNodeColor = 'tab:purple'
        self.hlNodeEdgeColor = 'tab:pink'
        self.edgeLabel = True
        self.disableReactionEdgeLabel = False
        self.disableModifierEdgeLabel = False
        self.edgeLabelFontSize = 10
        self.drawReactionNode = True
        self.breakBoundary = False
        self.tightLayout = False
        self.weights = []
        self.edgeTransparency = False
        self.plottingThreshold = 0.
        self.removeBelowThreshold = True
        self.analyzeFlux = False
        self.analyzeRates = False
        self.analyzeColorHigh = 'k'
        self.analyzeColorLow = 'k'
        self.analyzeColorMap = 'Reds'
        self.analyzeColorScale = False
        self.plotColorbar = False
        self.customAxis = None
        self.layoutAlgorithm = 'kamada-kawai'
        self._Var.pos = None
        
        
    def setLayout(self, pos):
        """
        Set custom layout and bypass whe is generated by the layout algorothm
        
        :param pos: Dictionary of all nodes and corresponding coordinates
        :type name: Dict
        """
        
        self._Var.pos = pos
        
    
    def getLayout(self):
        """
        Return the layout of the model
        """
        
        if self.layoutAlgorithm not in toolbox.getListOfAlgorithms():
            raise Exception("Unsupported layout algorithm: '" + str(self.layoutAlgorithm) + "'")
        
        avoid = ['C', 'CC', 'Ci', 'E1', 'EX', 'Ei', 'FF', 'GF', 'Ge', 'Gt', 'I', 'LC',
                 'LM', 'LT', 'Le', 'Li', 'Lt', 'N', 'Ne', 'O', 'Q', 'QQ', 'RR', 'S',
                 'Si', 'ZZ', 'ff', 'fu', 'im', 'jn', 'li', 'ln', 'oo', 'pi', 're',
                 'rf', 'yn']
        
        # extract reactant, product, modifiers, and kinetic laws
        allRxn = []
        allMod = []
        count = []
        count_mod = []
        self._Var.rid = []
        self._Var.r_type = []
        self._Var.mod_type = []
        allBoundary = []
        rid_ind = 0
        
        if len(self.weights) > 0 and len(self.weights) != len(self.rrInstances):
            raise Exception("The dimension of weights provides does not match "
                            "the number of models given")
    
        for rind, r in enumerate(self.rrInstances):
            numBnd = r.getNumBoundarySpecies()
            numFlt = r.getNumFloatingSpecies()
            boundaryId = self._Var.boundaryIds[rind]
            floatingId = self._Var.floatingIds[rind]
            
            allBoundary.append(boundaryId)
            
            # prepare symbols for sympy
            boundaryId_sympy = [] 
            floatingId_sympy = []
            
            # Fix issues with reserved characters
            for i in range(numBnd):
                if boundaryId[i] in avoid:
                    boundaryId_sympy.append('_' + boundaryId[i])
                else:
                    boundaryId_sympy.append(boundaryId[i])
            
            for i in range(numFlt):
                if floatingId[i] in avoid:
                    floatingId_sympy.append('_' + floatingId[i])
                else:
                    floatingId_sympy.append(floatingId[i])
                    
            paramIdsStr = ' '.join(r.getGlobalParameterIds())
            floatingIdsStr = ' '.join(floatingId_sympy)
            boundaryIdsStr = ' '.join(boundaryId_sympy)
            comparmentIdsStr = ' '.join(r.getCompartmentIds())
            
            allIds = paramIdsStr + ' ' + floatingIdsStr + ' ' + boundaryIdsStr + ' ' + comparmentIdsStr
            
            avsym = sympy.symbols(allIds)
            
            # extract reactant, product, modifiers, and kinetic laws
            rct = []
            prd = []
            mod_m = []
            mod_target = []
            kineticLaw = []
            r_type_m = []
            rr_type_m = []
            mod_type_m = []
            
            doc = libsbml.readSBMLFromString(r.getSBML())
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
                
                if len(temprct) == 0:
                    rct.append(['Input'])
                else:
                    rct.append(sorted(temprct, key=lambda v: (v.upper(), v[0].islower())))
                if len(tempprd) == 0:
                    prd.append(['Output'])
                else:
                    prd.append(sorted(tempprd, key=lambda v: (v.upper(), v[0].islower())))
                mod_m.append(sorted(tempmod, key=lambda v: (v.upper(), v[0].islower())))
                
                if sbmlreaction.getReversible():
                    rr_type_m.append('reversible')
                else:
                    rr_type_m.append('irreversible')
                
                if kl == None:
                    kineticLaw.append(None)
                else:
                    # Update kinetic law according to change in species name
                    kl_split = re.split('( |\(|\))', kl.getFormula())
                    for i in range(len(kl_split)):
                        if kl_split[i] in avoid:
                            kl_split[i] = '_' + kl_split[i]
                    
                    kineticLaw.append(''.join(kl_split))
            
            nkl = 0
            # use sympy for analyzing modifiers weSmart
            for ml in range(len(mod_m)):
                mod_type_temp = []
                expression = kineticLaw[ml]
                if expression == None:
                    for ml_i in range(len(mod_m[ml])):
                        mod_type_temp.append('modifier')
                    r_type_m.append(rr_type_m[nkl])
                    nkl += 1
                else:
                    n,d = sympy.fraction(expression)
                    for ml_i in range(len(mod_m[ml])):
                        if n.has(sympy.symbols(mod_m[ml][ml_i])):
                            mod_type_temp.append('activator')
                        elif d.has(sympy.symbols(mod_m[ml][ml_i])):
                            mod_type_temp.append('inhibitor')
                        else:
                            mod_type_temp.append('modifier')
                    n = '(' + str(n) + ')'
                    
                    # In case all products are in rate law, assume it is a reversible reaction
                    if (all(ext in str(n) for ext in [s + '/' for s in prd[ml]]) or
                        all(ext in str(n) for ext in [s + ')' for s in prd[ml]]) or
                        all(ext in str(n) for ext in [s + '*' for s in prd[ml]]) or
                        all(ext in str(n) for ext in [s + ';' for s in prd[ml]]) or
                        all(ext in str(n) for ext in [s + '+' for s in prd[ml]]) or
                        all(ext in str(n) for ext in [s + '-' for s in prd[ml]]) or
                        all(ext in str(n) for ext in [s + ' ' for s in prd[ml]])):
                        r_type_m.append('reversible')
                    else:
                        r_type_m.append('irreversible')
                mod_type_m.append(mod_type_temp)
            
            for i in range(len(mod_m)):
                if len(mod_m[i]) > 0:
                    for j in range(len(mod_m[i])):
                        mod_target.append([mod_m[i][j],[rct[i]],[prd[i]],mod_type_m[i][j]])

            speciesId = list(rct + prd)
            speciesId = [item for sublist in speciesId for item in sublist]
            self._Var.speciesId = list(set(speciesId))
            
            for t in range(sbmlmodel.getNumReactions()):
                if r_type_m[t] == 'irreversible':
                    if ([rct[t], prd[t]] not in allRxn):
                        allRxn.append([rct[t], prd[t]])
                        if len(self.weights) > 0:
                            count.append(1*self.weights[rind])
                        else:
                            count.append(1)
                        self._Var.rid.append("J" + str(rid_ind))
                        rid_ind += 1
                        self._Var.r_type.append(r_type_m[t])
                    else:
                        if len(self.weights) > 0:
                            count[allRxn.index([rct[t], prd[t]])] += 1*self.weights[rind]
                        else:
                            count[allRxn.index([rct[t], prd[t]])] += 1
                else:
                    if ([rct[t],prd[t]] not in allRxn) and ([prd[t],rct[t]] not in allRxn):
                        allRxn.append([rct[t], prd[t]])
                        if len(self.weights) > 0:
                            count.append(1*self.weights[rind])
                        else:
                            count.append(1)
                        self._Var.rid.append("J" + str(rid_ind))
                        rid_ind += 1
                        self._Var.r_type.append(r_type_m[t])
                    elif ([rct[t],prd[t]] not in allRxn) and ([prd[t],rct[t]] in allRxn):
                        if len(self.weights) > 0:
                            count[allRxn.index([prd[t], rct[t]])] += 1*self.weights[rind]
                        else:
                            count[allRxn.index([prd[t], rct[t]])] += 1
                    elif ([rct[t],prd[t]] in allRxn) and ([prd[t],rct[t]] not in allRxn):
                        if len(self.weights) > 0:
                            count[allRxn.index([rct[t], prd[t]])] += 1*self.weights[rind]
                        else:
                            count[allRxn.index([rct[t], prd[t]])] += 1                    

            for mi in range(len(mod_target)):
                # Reversible reactions
                try:
                    mod_ridx = allRxn.index([mod_target[mi][1][0], mod_target[mi][2][0]])
                except:
                    mod_ridx = allRxn.index([mod_target[mi][2][0], mod_target[mi][1][0]])

                if ([[mod_target[mi][0]], ['J'+str(mod_ridx)]] in allMod):
                    rxnind = [i for i, x in enumerate(allMod) if x == [[mod_target[mi][0]], ['J'+str(mod_ridx)]]]
                    if ([mod_target[mi][3]] in np.array(self._Var.mod_type)[rxnind]):
                        modind = [i for i, x in enumerate(self._Var.mod_type) if x == mod_target[mi][3]]
                        if len(self.weights) > 0:
                            count_mod[list(set(rxnind) & set(modind))[0]] += 1*self.weights[rind]
                        else:
                            count_mod[list(set(rxnind) & set(modind))[0]] += 1
                    else:
                        allMod.append([[mod_target[mi][0]], ['J'+str(mod_ridx)]])
                        if len(self.weights) > 0:
                            count_mod.append(1*self.weights[rind])
                        else:
                            count_mod.append(1)
                        if len(mod_target[mi][3]) > 0:
                            self._Var.mod_type.append(mod_target[mi][3])
                else:
                    allMod.append([[mod_target[mi][0]], ['J'+str(mod_ridx)]])
                    if len(self.weights) > 0:
                        count_mod.append(1*self.weights[rind])
                    else:
                        count_mod.append(1)
                    if len(mod_target[mi][3]) > 0:
                        self._Var.mod_type.append(mod_target[mi][3])
    
        self._Var.allRxn = allRxn + allMod
        count = count + count_mod
        
        # Break boundary
        self._Var.allBoundary = np.unique(allBoundary).tolist()
        if self.breakBoundary:
            speciesId_temp = []
            for i in range(len(self._Var.speciesId)):
                if self._Var.speciesId[i] not in allBoundary + ['Input', 'Output']:
                    speciesId_temp.append(self._Var.speciesId[i])
                
            self._Var.speciesId = speciesId_temp
            
            boundaryId_temp = []
            bc = 0
            for i in range(len(self._Var.allRxn)):
                for j in range(len(self._Var.allRxn[i][0])):
                    if self._Var.allRxn[i][0][j] in allBoundary + ['Input', 'Output']:
                        self._Var.allRxn[i][0][j] = self._Var.allRxn[i][0][j] + '_' + str(bc)
                        self._Var.speciesId.append(self._Var.allRxn[i][0][j])
                        boundaryId_temp.append(self._Var.allRxn[i][0][j])
                        bc += 1
                for k in range(len(self._Var.allRxn[i][1])):
                    if self._Var.allRxn[i][1][k] in allBoundary + ['Input', 'Output']:
                        self._Var.allRxn[i][1][k] = self._Var.allRxn[i][1][k] + '_' + str(bc)
                        self._Var.speciesId.append(self._Var.allRxn[i][1][k])
                        boundaryId_temp.append(self._Var.allRxn[i][1][k])
                        bc += 1
            self._Var.allBoundary = boundaryId_temp
            
        self._Var.count = np.divide(count, len(self.rrInstances))
        
        # initialize directional graph
        self._Var.G = nx.DiGraph()
    
        # add edges
        sid_used = []
        self._Var.rid_used = []
        rid_idx = 0
        
        for i in range(len(self._Var.allRxn)):
            if self._Var.allRxn[i][1][0] not in self._Var.rid:
                if count[i] > self.plottingThreshold:
                    for k in range(len(self._Var.allRxn[i][0])):
                        self._Var.G.add_edges_from([(self._Var.allRxn[i][0][k], self._Var.rid[rid_idx])])
                        
                    for j in range(len(self._Var.allRxn[i][1])):
                        self._Var.G.add_edges_from([(self._Var.rid[rid_idx], self._Var.allRxn[i][1][j])])
                    
                    sid_used.append(self._Var.allRxn[i][0][k])
                    sid_used.append(self._Var.allRxn[i][1][j])
                    self._Var.rid_used.append(self._Var.rid[rid_idx])
                    
                rid_idx += 1
            else:
                if count[i] > self.plottingThreshold:
                    self._Var.G.add_edges_from([(self._Var.allRxn[i][0][0], self._Var.allRxn[i][1][0])])
                    sid_used.append(self._Var.allRxn[i][0][0])
        
        self._Var.sid_used = np.unique(sid_used).tolist()
        
        # calcutate positions
        thres = 0.3
        if self.layoutAlgorithm == 'spring':
            pos = nx.spring_layout(self._Var.G, scale=self.scale, seed=1)
        elif self.layoutAlgorithm == 'kamada-kawai':
            shortest_dist = dict(nx.shortest_path_length(self._Var.G, weight='weight'))
            pos = nx.kamada_kawai_layout(self._Var.G, dist=shortest_dist, scale=self.scale)
            
            maxIter = 5
            maxIter_n = 0
            dist_flag = True
            
            if self.tightLayout:
                comId = self._Var.speciesId
            else:
                comId = self._Var.speciesId + self._Var.rid
            
            while dist_flag and (maxIter_n < maxIter):
                dist_flag = False
                for i in itertools.combinations(pos.keys(), 2):
                    pos_dist = np.sqrt((pos[i[0]][0] - pos[i[1]][0])**2 + (pos[i[0]][1] - pos[i[1]][1])**2)
                    if pos_dist < thres:
                        dist_flag = True
                        shortest_dist[i[0]][i[1]] = 2
                        shortest_dist[i[1]][i[0]] = 2
                pos = nx.kamada_kawai_layout(self._Var.G, dist=shortest_dist, scale=self.scale)
                maxIter_n += 1
        elif self.layoutAlgorithm == 'twopi' or self.layoutAlgorithm == 'neato' or self.layoutAlgorithm == 'dot':
            from networkx.drawing.nx_pydot import graphviz_layout
            try:
                pos = graphviz_layout(self._Var.G, prog=self.layoutAlgorithm)
            except:
                raise Exception("Error running graphviz: Please check graphviz is properly configured.")
            keylist = np.array(list(pos.keys()))
            poslist = np.array(list(pos.values()))
            poslist /= np.max(np.abs(poslist),axis=0)
            pos = dict(zip(keylist, poslist))
        else:
            raise Exception("Unsupported layout algorithm.")
            
        return pos
    
    
    def drawWeightedDiagram(self, show=True, savePath=None, dpi=150):
        """     
        Draw weighted reaction network based on frequency of reactions
        
        :param show: flag to show the diagram
        :param savePath: path to save the diagram
        :param dpi: dpi settings for the diagram
        :returns allRxn: list of all reactions in the list of models presented as a pair of reactants and products
        :returns count: normalized count of reactions in allRxn throughout the list of models
        """
        
        if self._Var.pos == None:
            pos = self.getLayout()
        else:
            pos = self.getLayout()
            pos = self._Var.pos
            assert(len(self._Var.boundaryIds[0])+len(self._Var.floatingIds[0])+len(self._Var.rids[0]) == len(pos))
        
        if not self.removeBelowThreshold:
            rid_idx = 0
            for i in range(len(self._Var.allRxn)):
                if self._Var.allRxn[i][1][0] not in self._Var.rid:
                    if self._Var.count[i] <= self.plottingThreshold:
                        for k in range(len(self._Var.allRxn[i][0])):
                            self._Var.G.add_edges_from([(self._Var.allRxn[i][0][k], self._Var.rid[rid_idx])])
                            
                        for j in range(len(self._Var.allRxn[i][1])):
                            self._Var.G.add_edges_from([(self._Var.rid[rid_idx], self._Var.allRxn[i][1][j])])
                        
                    rid_idx += 1
                else:
                    if self._Var.count[i] <= self.plottingThreshold:
                        self._Var.G.add_edges_from([(self._Var.allRxn[i][0][0], self._Var.allRxn[i][1][0])])
            
            pos = nx.spring_layout(self._Var.G, pos=pos, fixed=self._Var.sid_used+self._Var.rid_used, scale=self.scale, seed=1)
            self._Var.sid_used = self._Var.speciesId
            self._Var.rid_used = self._Var.rid
        
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
        if self.customAxis != None:
            ax = self.customAxis
        else:
            ax = plt.gca()
        
        # add nodes to the figure
        for n in self._Var.G:
            if n in self._Var.rid:
                rec_width = 0.05*(self.fontsize/20)
                rec_height = 0.05*(self.fontsize/20)
                if n in self.highlight:
                    c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                        pos[n][1]-rec_height/2),
                                        rec_width, 
                                        rec_height,
                                        boxstyle="round,pad=0.01, rounding_size=0.01",
                                        linewidth=self.nodeEdgelw, 
                                        edgecolor=self.hlNodeEdgeColor, 
                                        facecolor=self.hlNodeColor)
                else:
                    c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                        pos[n][1]-rec_height/2),
                                        rec_width, 
                                        rec_height,
                                        boxstyle="round,pad=0.01, rounding_size=0.01",
                                        linewidth=self.nodeEdgelw, 
                                        edgecolor=self.nodeEdgeColor, 
                                        facecolor=self.reactionNodeColor)
                if self.labelReactionIds:
                    ax.text(pos[n][0], 
                            pos[n][1], 
                            n, 
                            fontsize=self.fontsize, 
                            horizontalalignment='center', 
                            verticalalignment='center', 
                            color=self.labelColor)
            else:
                if len(n) > 10:
                    rec_width = max(0.045*((len(n)/2)+1), 0.13)*(self.fontsize/20)
                    rec_height = 0.20*(self.fontsize/20)
                else:
                    rec_width = max(0.045*(len(n)+1), 0.13)*(self.fontsize/20)
                    rec_height = 0.11*(self.fontsize/20)
                    
                if (n in self._Var.allBoundary) or (n == 'Input') or (n == 'Output'):
                    node_color = self.boundaryColor
                else:
                    node_color = self.nodeColor
                    
                if n in self.highlight:
                    c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                        pos[n][1]-rec_height/2),
                                        rec_width, 
                                        rec_height,
                                        boxstyle="round,pad=0.01, rounding_size=0.02",
                                        linewidth=self.nodeEdgelw, 
                                        edgecolor=self.hlNodeEdgeColor, 
                                        facecolor=self.hlNodeColor)
                else:
                    c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                        pos[n][1]-rec_height/2),
                                        rec_width, 
                                        rec_height,
                                        boxstyle="round,pad=0.01, rounding_size=0.02",
                                        linewidth=self.nodeEdgelw, 
                                        edgecolor=self.nodeEdgeColor, 
                                        facecolor=node_color)
                if len(n) > 10:
                    ax.text(pos[n][0], pos[n][1], n[:int(len(n)/2)] + '\n' + n[int(len(n)/2):], 
                            fontsize=self.fontsize, horizontalalignment='center', 
                            verticalalignment='center', color=self.labelColor)
                else:
                    ax.text(pos[n][0], pos[n][1], n, 
                            fontsize=self.fontsize, horizontalalignment='center', 
                            verticalalignment='center', color=self.labelColor)
            self._Var.G.nodes[n]['patch'] = c
        
        # add edges to the figure
        mod_idx = 0
        rid_idx = 0
        rad_track = True
        
        seen={}
        for i in range(len(self._Var.allRxn)):
            if self._Var.allRxn[i][1][0] not in self._Var.rid:
                if self._Var.count[i] > self.plottingThreshold or not self.removeBelowThreshold:
                    if (len(self._Var.allRxn[i][0]) == 1) or (len(self._Var.allRxn[i][1]) == 1): # UNI-involved
                        comb = list(itertools.combinations_with_replacement(self._Var.allRxn[i][0],len(self._Var.allRxn[i][1])))
                        for j in [list(zip(x,self._Var.allRxn[i][1])) for x in comb]:
                            for k in range(len(j)):
                                p1 = self._Var.G.nodes[j[k][0]]['patch']
                                p2 = self._Var.G.nodes[self._Var.rid[rid_idx]]['patch']
                                p3 = self._Var.G.nodes[j[k][1]]['patch']
                    
                                X1 = (p1.get_x()+p1.get_width()/2,p1.get_y()+p1.get_height()/2)
                                X2 = (p2.get_x()+p2.get_width()/2,p2.get_y()+p2.get_height()/2)
                                X3 = (p3.get_x()+p3.get_width()/2,p3.get_y()+p3.get_height()/2)
                                
                                if ((len(np.unique(self._Var.allRxn[i][0])) > len(self._Var.allRxn[i][1])) or 
                                    (len(self._Var.allRxn[i][0]) < len(np.unique(self._Var.allRxn[i][1])))): # Uni-Bi or Bi-Uni
                                    XY1 = np.vstack((X1, X2))
                                    XY2 = np.vstack((X2, X3))
                                    
                                    tck1, u1 = interpolate.splprep([XY1[:,0], XY1[:,1]], 
                                                                   k=1)
                                    intX1, intY1 = interpolate.splev(np.linspace(0, 1, 100),
                                                                     tck1, 
                                                                     der=0)
                                    stackXY = np.vstack((intX1, intY1))
                                    tck2, u2 = interpolate.splprep([XY2[:,0], XY2[:,1]], 
                                                                   k=1)
                                    intX2, intY2 = interpolate.splev(np.linspace(0, 1, 100), 
                                                                     tck2, 
                                                                     der=0)
                                    stackXY2 = np.vstack((intX2, intY2))
                                    
                                    if max(stackXY[0]) > max_width[1]:
                                        max_width[1] = max(stackXY[0])
                                    if min(stackXY[0]) < max_width[0]:
                                        max_width[0] = min(stackXY[0])
                                    if max(stackXY[1]) > max_height[1]:
                                        max_height[1] = max(stackXY[1])
                                    if min(stackXY[1]) < max_height[0]:
                                        max_height[0] = min(stackXY[1])
                                        
                                    if max(stackXY2[0]) > max_width[1]:
                                        max_width[1] = max(stackXY2[0])
                                    if min(stackXY2[0]) < max_width[0]:
                                        max_width[0] = min(stackXY2[0])
                                    if max(stackXY2[1]) > max_height[1]:
                                        max_height[1] = max(stackXY2[1])
                                    if min(stackXY2[1]) < max_height[0]:
                                        max_height[0] = min(stackXY2[1])
                                    
                                    X3top = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y()+p3.get_height())
                                    X3bot = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y())
                                    X3left = (p3.get_x(),
                                              p3.get_y()+p3.get_height()/2)
                                    X3right = (p3.get_x()+p3.get_width(),
                                               p3.get_y()+p3.get_height()/2)
                                    
                                    n = -1
                                    arrthres_h = .02
                                    arrthres_v = .01
                                    while (((stackXY2.T[n][0] > (X3left[0]-arrthres_h)) and
                                            (stackXY2.T[n][0] < (X3right[0]+arrthres_h)) and
                                            (stackXY2.T[n][1] > (X3bot[1]-arrthres_v)) and 
                                            (stackXY2.T[n][1] < (X3top[1]+arrthres_v))) and
                                            (np.abs(n) < np.shape(stackXY2)[1] - 10)):
                                        n -= 1
                                   
                                    lpath1 = Path(stackXY.T)
                                    lpath2 = Path(stackXY2.T[:n])
                                    
                                    if self.edgeTransparency:
                                        alpha = self._Var.count[i]
                                    else:
                                        alpha = None
                                    
                                    arrowstyle1 = ArrowStyle.Curve()
                                    arrowstyle2 = ArrowStyle.CurveFilledB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                          head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                                    
                                    if self._Var.r_type[i] == 'reversible':
                                        lpath1 = Path(stackXY.T[-n:])
                                        arrowstyle1 = ArrowStyle.CurveFilledA(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                              head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                                        arrowstyle2 = ArrowStyle.CurveFilledB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                              head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                                    
                                    e1 = FancyArrowPatch(path=lpath1,
                                                        arrowstyle=arrowstyle1,
                                                        mutation_scale=10.0,
                                                        lw=(self._Var.count[i]*self.edgelw),
                                                        alpha=alpha,
                                                        color=self.reactionColor)
                                    
                                    e2 = FancyArrowPatch(path=lpath2,
                                                        arrowstyle=arrowstyle2,
                                                        mutation_scale=10.0,
                                                        lw=(self._Var.count[i]*self.edgelw),
                                                        alpha=alpha,
                                                        color=self.reactionColor)
                                    
                                    ax.add_patch(e1)
                                    ax.add_patch(e2)
                                    
                                else: # Uni-Uni    
                                    XY = np.vstack((X1, X2, X3))
                                    
                                    tck, u = interpolate.splprep([XY[:,0], XY[:,1]], k=2)
                                    intX, intY = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
                                    stackXY = np.vstack((intX, intY))
                                    
                                    if max(stackXY[0]) > max_width[1]:
                                        max_width[1] = max(stackXY[0])
                                    if min(stackXY[0]) < max_width[0]:
                                        max_width[0] = min(stackXY[0])
                                    if max(stackXY[1]) > max_height[1]:
                                        max_height[1] = max(stackXY[1])
                                    if min(stackXY[1]) < max_height[0]:
                                        max_height[0] = min(stackXY[1])
                                    
                                    X3top = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y()+p3.get_height())
                                    X3bot = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y())
                                    X3left = (p3.get_x(),
                                              p3.get_y()+p3.get_height()/2)
                                    X3right = (p3.get_x()+p3.get_width(),
                                               p3.get_y()+p3.get_height()/2)
                                    
                                    n = -1
                                    arrthres_h = .02
                                    arrthres_v = .01
                                    while (((stackXY.T[n][0] > (X3left[0]-arrthres_h)) and 
                                            (stackXY.T[n][0] < (X3right[0]+arrthres_h)) and
                                            (stackXY.T[n][1] > (X3bot[1]-arrthres_v)) and 
                                            (stackXY.T[n][1] < (X3top[1]+arrthres_v))) and
                                            (np.abs(n) < np.shape(stackXY)[1] - 10)):
                                        n -= 1
                                   
                                    lpath = Path(stackXY.T[:n])
                                    
                                    if self.edgeTransparency:
                                        alpha = self._Var.count[i]
                                    else:
                                        alpha = None
                                    
                                    arrowstyle = ArrowStyle.CurveFilledB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                          head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                                    
                                    if self._Var.r_type[i] == 'reversible':
                                        lpath = Path(stackXY.T[-n:n])
                                        arrowstyle = ArrowStyle.CurveFilledAB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                              head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                                    
                                    e = FancyArrowPatch(path=lpath,
                                                        arrowstyle=arrowstyle,
                                                        mutation_scale=10.0,
                                                        lw=(self._Var.count[i]*self.edgelw),
                                                        alpha=alpha,
                                                        color=self.reactionColor)
                                    ax.add_patch(e)
                    else: # BIBI or larger
                        if len(self._Var.allRxn[i][0]) < len(self._Var.allRxn[i][1]):
                            rVal = len(self._Var.allRxn[i][0])
                        else:
                            rVal = len(self._Var.allRxn[i][1])
                            
                        for j in [list(zip(x,self._Var.allRxn[i][1])) for x in itertools.combinations(self._Var.allRxn[i][0],rVal)][0]:
                            p1 = self._Var.G.nodes[j[0]]['patch']
                            p2 = self._Var.G.nodes[self._Var.rid[rid_idx]]['patch']
                            p3 = self._Var.G.nodes[j[1]]['patch']
                            
                            X1 = (p1.get_x()+p1.get_width()/2,p1.get_y()+p1.get_height()/2)
                            X2 = (p2.get_x()+p2.get_width()/2,p2.get_y()+p2.get_height()/2)
                            X3 = (p3.get_x()+p3.get_width()/2,p3.get_y()+p3.get_height()/2)
                            
                            XY = np.vstack((X1, X2, X3))
                            
                            tck, u = interpolate.splprep([XY[:,0], XY[:,1]], k=2)
                            intX, intY = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
                            stackXY = np.vstack((intX, intY))
                            
                            if max(stackXY[0]) > max_width[1]:
                                max_width[1] = max(stackXY[0])
                            if min(stackXY[0]) < max_width[0]:
                                max_width[0] = min(stackXY[0])
                            if max(stackXY[1]) > max_height[1]:
                                max_height[1] = max(stackXY[1])
                            if min(stackXY[1]) < max_height[0]:
                                max_height[0] = min(stackXY[1])
                            
                            X3top = (p3.get_x()+p3.get_width()/2,
                                     p3.get_y()+p3.get_height())
                            X3bot = (p3.get_x()+p3.get_width()/2,
                                     p3.get_y())
                            X3left = (p3.get_x(),
                                      p3.get_y()+p3.get_height()/2)
                            X3right = (p3.get_x()+p3.get_width(),
                                       p3.get_y()+p3.get_height()/2)
                            
                            n = -1
                            arrthres_h = .02
                            arrthres_v = .01
                            while (((stackXY.T[n][0] > (X3left[0]-arrthres_h)) and
                                    (stackXY.T[n][0] < (X3right[0]+arrthres_h)) and
                                    (stackXY.T[n][1] > (X3bot[1]-arrthres_v)) and 
                                    (stackXY.T[n][1] < (X3top[1]+arrthres_v))) and
                                    (np.abs(n) < np.shape(stackXY)[1] - 10)):
                                n -= 1
                           
                            lpath = Path(stackXY.T[:n])
                            
                            if self.edgeTransparency:
                                alpha = self._Var.count[i]
                            else:
                                alpha = None
                            
                            arrowstyle = ArrowStyle.CurveFilledB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                 head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                            
                            if self._Var.r_type[i] == 'reversible':
                                lpath = Path(stackXY.T[-n:n])
                                arrowstyle = ArrowStyle.CurveFilledAB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                      head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                            e = FancyArrowPatch(path=lpath,
                                                arrowstyle=arrowstyle,
                                                mutation_scale=10.0,
                                                lw=(self._Var.count[i]*self.edgelw),
                                                alpha=alpha,
                                                color=self.reactionColor)
                            ax.add_patch(e)
                    # Edge labels
                    if self.edgeLabel:
                        if not self.disableReactionEdgeLabel:
                            c = FancyBboxPatch((stackXY.T[50,0]-0.0325, stackXY.T[50,1]+0.005),
                                               0.125, 
                                               0.05,
                                               boxstyle="round,pad=0.01, rounding_size=0.01",
                                               color='w')
                            ax.add_patch(c)
                            ax.text(stackXY.T[50,0]+0.03, stackXY.T[50,1]+0.03, round(self._Var.count[i], 3), 
                                    fontsize=self.edgeLabelFontSize, horizontalalignment='center', 
                                    verticalalignment='center')
                rid_idx += 1
            else:
                # Modifiers
                if self._Var.count[i] > self.plottingThreshold or not self.removeBelowThreshold:
                    for m, e in enumerate(self._Var.allRxn[i][0]):
                        n1 = self._Var.G.nodes[e]['patch']
                        n2 = self._Var.G.nodes[self._Var.allRxn[i][1][0]]['patch']
                        rad = 0.15
                        shrinkB = 2.
                        
                        if (self._Var.allRxn[i][0][0],self._Var.allRxn[i][1][0]) in seen:
                            rad = seen.get((e,self._Var.allRxn[i][1][0])) # TODO: No curvature when there is just a single line between two nodes
                            rad = (rad+np.sign(rad)*0.15)*-1 # TODO: Change curvature
                            
                        X1 = (n1.get_x()+n1.get_width()/2,
                              n1.get_y()+n1.get_height()/2)
                        X2 = (n2.get_x()+n2.get_width()/2,
                              n2.get_y()+n2.get_height()/2)
                        
                        XY = np.vstack((X1, X2))
                        
                        tck, u = interpolate.splprep([XY[:,0], XY[:,1]], k=1)
                        intX, intY = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
                        stackXY = np.vstack((intX, intY))
                        
                        if self._Var.mod_type[mod_idx] == 'inhibitor': # inhibition
                            color = self.modifierColor
                            arrowstyle = ArrowStyle.BarAB(widthA=0.0, angleA=None, widthB=1.0, angleB=None)
                            shrinkB = 10.
                            linestyle = '-'
                        elif self._Var.mod_type[mod_idx] == 'activator': # activation
                            color = self.modifierColor
                            arrowstyle = ArrowStyle.CurveFilledB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                 head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                            linestyle = '-'
                        elif self._Var.mod_type[mod_idx] == 'modifier': # Unknown modifier
                            color = self.modifierColor
                            arrowstyle = ArrowStyle.CurveFilledB(head_length=(0.8 + 0.01*self._Var.count[i]*self.edgelw), 
                                                                 head_width=(0.4 + 0.01*self._Var.count[i]*self.edgelw))
                            linestyle = ':'
                            
                        if self.edgeTransparency:
                            alpha = self._Var.count[i]
                        else:
                            alpha = None
                        
                        e = FancyArrowPatch(X1,
                                            X2,
                                            patchA=n1,
                                            patchB=n2,
                                            shrinkB=shrinkB,
                                            arrowstyle=arrowstyle,
                                            connectionstyle='arc3,rad=%s'%rad,
                                            mutation_scale=10.0,
                                            lw=(self._Var.count[i]*self.edgelw),
                                            color=color,
                                            alpha=alpha,
                                            linestyle=linestyle)
                        
                        seen[(self._Var.allRxn[i][0][0],self._Var.allRxn[i][1][0])]=rad
                        ax.add_patch(e)
                        
                        # Edge labels
                        if self.edgeLabel:
                            if not self.disableModifierEdgeLabel:
                                verts = e.get_path().vertices
                                trans = e.get_patch_transform()
                                points = trans.transform(verts)
                                if self._Var.mod_type[mod_idx] == 'inhibitor':
                                    c = FancyBboxPatch((stackXY.T[50,0]+((points[5,0]-stackXY.T[50,0])/2)-0.0625,
                                                        stackXY.T[50,1]+((points[5,1]-stackXY.T[50,1])/2)-0.025),
                                                       0.125, 
                                                       0.05,
                                                       boxstyle="round,pad=0.01, rounding_size=0.01",
                                                       color='w')
                                    ax.add_patch(c)
                                    ax.text(stackXY.T[50,0]+((points[5,0]-stackXY.T[50,0])/2),
                                            stackXY.T[50,1]+((points[5,1]-stackXY.T[50,1])/2),
                                            round(self._Var.count[i], 3), 
                                            fontsize=self.edgeLabelFontSize, 
                                            horizontalalignment='center', 
                                            verticalalignment='center', 
                                            color='r')
                                elif self._Var.mod_type[mod_idx] == 'activator':
                                    plt.plot(points[:1,0],points[:1,1])
                                    c = FancyBboxPatch((stackXY.T[50,0]+((points[1,0]-stackXY.T[50,0])/2)-0.0625,
                                                        stackXY.T[50,1]+((points[1,1]-stackXY.T[50,1])/2)-0.025),
                                                       0.125, 
                                                       0.05,
                                                       boxstyle="round,pad=0.01, rounding_size=0.01",
                                                       color='w')
                                    ax.add_patch(c)
                                    ax.text(stackXY.T[50,0]+((points[1,0]-stackXY.T[50,0])/2),
                                            stackXY.T[50,1]+((points[1,1]-stackXY.T[50,1])/2),
                                            round(self._Var.count[i], 3), 
                                            fontsize=self.edgeLabelFontSize, 
                                            horizontalalignment='center', 
                                            verticalalignment='center', 
                                            color='r')
                    
                mod_idx += 1
            
        # Add nodes at last to put it on top
        if self.drawReactionNode:
            allnodes = self._Var.sid_used + self._Var.rid_used
        else:
            allnodes = self._Var.sid_used
            
        if 'Input' in self._Var.G.nodes:
            allnodes += ['Input']
        if 'Output' in self._Var.G.nodes:
            allnodes += ['Output']
        for i in range(len(allnodes)):
            ax.add_patch(self._Var.G.nodes[allnodes[i]]['patch'])
        
        # reset width and height
        ax.autoscale()
        fig.set_figwidth((abs(max_width[0] - max_width[1])+0.5)*5)
        fig.set_figheight((abs(max_height[0] - max_height[1])+0.5)*5)
        plt.axis('off')
        plt.axis('equal')
        
        if savePath != None:
            try:
                fig.savefig(savePath, bbox_inches='tight', dpi=dpi)
            except IOError as e:
                raise Exception("Error saving diagram: " + str(e))
            return self._Var.allRxn, self._Var.count
        else:
            if show and self.customAxis == None:
                plt.show()
            plt.close()
            return self._Var.allRxn, self._Var.count
    
    
    def drawNetworkGrid(self, nrows, ncols, auto=False, show=True, savePath=None, dpi=150):
        """
        Plot a grid of network diagrams
        
        :param nrows: number of rows
        :param ncols: number of columns
        :param auto: Automatically configure nrows and ncols based on the number of models. Overrides nrows and ncols.
        :param show: flag to show the diagram
        :param savePath: path to save the diagram
        :param dpi: dpi settings for the diagram
        """
        
        if self.layoutAlgorithm not in toolbox.getListOfAlgorithms():
            raise Exception("Unsupported layout algorithm: '" + str(self.layoutAlgorithm) + "'")
        
        edgelw_backup = self.edgelw
        
        fig, ax = plt.subplots(nrows, ncols, squeeze=False, sharex=True, sharey=True)
        fig.set_figheight(7*nrows)
        fig.set_figwidth(7*ncols)
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        
        for mdl in range(nrows*ncols):
            plt.sca(fig.axes[mdl])
            fig.axes[mdl].axis('off')
            
            if mdl < len(self._Var.models):
                net = Network(self._Var.models[mdl])
                
                if self._Var.pos == None:
                    pos, iVar = net.getLayout(returnState=True)
                else:
                    pos = self._Var.pos
                    assert(len(self._Var.boundaryId)+len(self._Var.floatingId)+len(self._Var.rid) == len(pos))
        
                # check the range of x and y positions
                max_width = []
                max_height = []
                for key, value in pos.items():
                    max_width.append(value[0])
                    max_height.append(value[1])
                
                max_width = [min(max_width), max(max_width)]
                max_height = [min(max_height), max(max_height)]
                
                # add nodes to the figure
                for n in iVar.G:
                    if n in iVar.rid:
                        rec_width = 0.05*(self.fontsize/20)
                        rec_height = 0.05*(self.fontsize/20)
                        if n in self.highlight:
                            c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                                pos[n][1]-rec_height/2),
                                                rec_width, 
                                                rec_height,
                                                boxstyle="round,pad=0.01, rounding_size=0.01",
                                                linewidth=self.nodeEdgelw, 
                                                edgecolor=self.hlNodeEdgeColor, 
                                                facecolor=self.hlNodeColor)
                        else:
                            c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                                pos[n][1]-rec_height/2), 
                                                rec_width, 
                                                rec_height,
                                                boxstyle="round,pad=0.01, rounding_size=0.01",
                                                linewidth=self.nodeEdgelw, 
                                                edgecolor=self.nodeEdgeColor, 
                                                facecolor=self.reactionNodeColor)
                        if self.labelReactionIds:
                            plt.text(pos[n][0], 
                                     pos[n][1], 
                                     n, 
                                     fontsize=self.fontsize, 
                                     horizontalalignment='center', 
                                     verticalalignment='center', 
                                     color=self.labelColor)
                    else:
                        if len(n) > 10:
                            rec_width = max(0.045*((len(n)/2)+1), 0.15)*(self.fontsize/20)
                            rec_height = 0.27*(self.fontsize/20)
                        else:
                            rec_width = max(0.045*(len(n)+1), 0.15)*(self.fontsize/20)
                            rec_height = 0.13*(self.fontsize/20)
                            
                        if (n in iVar.boundaryId) or (n == 'Input') or (n == 'Output'):
                            node_color = self.boundaryColor
                        else:
                            node_color = self.nodeColor
                            
                        if n in self.highlight:
                            c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                                pos[n][1]-rec_height/2),
                                                rec_width, 
                                                rec_height,
                                                boxstyle="round,pad=0.01, rounding_size=0.02",
                                                linewidth=self.nodeEdgelw, 
                                                edgecolor=self.hlNodeEdgeColor, 
                                                facecolor=self.hlNodeColor)
                        else:
                            c = FancyBboxPatch((pos[n][0]-rec_width/2, 
                                                pos[n][1]-rec_height/2), 
                                                rec_width, 
                                                rec_height,
                                                boxstyle="round,pad=0.01, rounding_size=0.02",
                                                linewidth=self.nodeEdgelw, 
                                                edgecolor=self.nodeEdgeColor, 
                                                facecolor=node_color)
                        if len(n) > 10:
                            plt.text(pos[n][0], pos[n][1], n[:int(len(n)/2)] + '\n' + n[int(len(n)/2):], 
                                     fontsize=self.fontsize, horizontalalignment='center', 
                                     verticalalignment='center', color=self.labelColor)
                        else:
                            plt.text(pos[n][0], pos[n][1], n, 
                                     fontsize=self.fontsize, horizontalalignment='center', 
                                     verticalalignment='center', color=self.labelColor)
                    iVar.G.nodes[n]['patch'] = c
                
                # add edges to the figure
                for i in range(len(iVar.rid)):
                    if (len(iVar.rct[i]) == 1) or (len(iVar.prd[i]) == 1): # UNI-involved
                        comb = list(itertools.combinations_with_replacement(iVar.rct[i],len(iVar.prd[i])))
                        for j in [list(zip(x,iVar.prd[i])) for x in comb]:
                            for k in range(len(j)):
                                p1 = iVar.G.nodes[j[k][0]]['patch']
                                p2 = iVar.G.nodes[iVar.rid[i]]['patch']
                                p3 = iVar.G.nodes[j[k][1]]['patch']
                                
                                X1 = (p1.get_x()+p1.get_width()/2,p1.get_y()+p1.get_height()/2)
                                X2 = (p2.get_x()+p2.get_width()/2,p2.get_y()+p2.get_height()/2)
                                X3 = (p3.get_x()+p3.get_width()/2,p3.get_y()+p3.get_height()/2)
                                
                                if ((len(np.unique(iVar.rct[i])) > len(iVar.prd[i])) or 
                                    (len(iVar.rct[i]) < len(np.unique(iVar.prd[i])))): # Uni-Bi or Bi-Uni
                                    XY1 = np.vstack((X1, X2))
                                    XY2 = np.vstack((X2, X3))
                                    
                                    tck1, u1 = interpolate.splprep([XY1[:,0], XY1[:,1]], 
                                                                   k=1)
                                    intX1, intY1 = interpolate.splev(np.linspace(0, 1, 100),
                                                                     tck1, 
                                                                     der=0)
                                    stackXY1 = np.vstack((intX1, intY1))
                                    tck2, u2 = interpolate.splprep([XY2[:,0], XY2[:,1]], 
                                                                   k=1)
                                    intX2, intY2 = interpolate.splev(np.linspace(0, 1, 100), 
                                                                     tck2, 
                                                                     der=0)
                                    stackXY2 = np.vstack((intX2, intY2))
                                    
                                    X3top = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y()+p3.get_height())
                                    X3bot = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y())
                                    X3left = (p3.get_x(),
                                              p3.get_y()+p3.get_height()/2)
                                    X3right = (p3.get_x()+p3.get_width(),
                                               p3.get_y()+p3.get_height()/2)
                                    
                                    n_1 = -1
                                    arrthres_h = .02
                                    arrthres_v = .02
                                    while (((stackXY2.T[n_1][0] > (X3left[0]-arrthres_h)) and
                                            (stackXY2.T[n_1][0] < (X3right[0]+arrthres_h)) and
                                            (stackXY2.T[n_1][1] > (X3bot[1]-arrthres_v)) and 
                                            (stackXY2.T[n_1][1] < (X3top[1]+arrthres_v))) and
                                            (np.abs(n_1) < np.shape(stackXY2)[1] - 75)):
                                        n_1 -= 1
                                    
                                    lpath1 = Path(stackXY1.T)
                                    lpath2 = Path(stackXY2.T[:n_1])
                                    lw1 = (1+self.edgelw)
                                    lw2 = (1+self.edgelw)
                                    arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                    arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                    
                                    if iVar.r_type[i] == 'reversible':
                                        X1top = (p1.get_x()+p1.get_width()/2,
                                                 p1.get_y()+p1.get_height())
                                        X1bot = (p1.get_x()+p1.get_width()/2,
                                                 p1.get_y())
                                        X1left = (p1.get_x(),
                                                  p1.get_y()+p1.get_height()/2)
                                        X1right = (p1.get_x()+p1.get_width(),
                                                   p1.get_y()+p1.get_height()/2)
                                        
                                        n_2 = 0
                                        while (((stackXY1.T[n_2][0] > (X1left[0]-arrthres_h)) and 
                                                (stackXY1.T[n_2][0] < (X1right[0]+arrthres_h)) and
                                                (stackXY1.T[n_2][1] > (X1bot[1]-arrthres_v)) and 
                                                (stackXY1.T[n_2][1] < (X1top[1]+arrthres_v))) and
                                                (np.abs(n_2) < np.shape(stackXY1)[1] - 75)):
                                            n_2 += 1
                                        
                                        lpath1 = Path(stackXY1.T[n_2:])
                                        
                                        if self.analyzeFlux:
                                            if iVar.flux[i] > 0:
                                                lw1 = (1+self.edgelw)
                                                lw2 = (4+self.edgelw)
                                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=1.2, head_width=0.8)
                                            elif self._Var.flux[i] < 0:
                                                lw1 = (4+self.edgelw)
                                                lw2 = (1+self.edgelw)
                                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=1.2, head_width=0.8)
                                    else:
                                        arrowstyle1 = ArrowStyle.Curve()
                                        
                                    e1 = FancyArrowPatch(path=lpath1,
                                                        arrowstyle=arrowstyle1,
                                                        mutation_scale=10.0,
                                                        lw=lw1,
                                                        color=self.reactionColor)
                                    
                                    e2 = FancyArrowPatch(path=lpath2,
                                                        arrowstyle=arrowstyle2,
                                                        mutation_scale=10.0,
                                                        lw=lw2,
                                                        color=self.reactionColor)
                                        
                                    fig.axes[mdl].add_patch(e1)
                                    fig.axes[mdl].add_patch(e2)
                                    
                                    if j[k][0] in iVar.floatingId:
                                        if (np.abs(iVar.stoch[iVar.stoch_row.index(j[k][0])][i]) > 1):
                                            # position calculation
                                            slope = ((lpath1.vertices[0][1] - lpath1.vertices[10][1])/
                                                     (lpath1.vertices[0][0] - lpath1.vertices[10][0]))
                                            x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                            y_prime = -slope*x_prime
                                            plt.text(x_prime+lpath1.vertices[10][0], 
                                                     y_prime+lpath1.vertices[10][1], 
                                                     int(np.abs(iVar.stoch[iVar.stoch_row.index(j[k][0])][i])), 
                                                     fontsize=self.fontsize, 
                                                     horizontalalignment='center', 
                                                     verticalalignment='center', 
                                                     color=self.reactionColor)
                                    
                                    if j[k][1] in iVar.floatingId:
                                        if (np.abs(iVar.stoch[iVar.stoch_row.index(j[k][1])][i]) > 1):
                                            slope = ((lpath2.vertices[0][1] - lpath2.vertices[-20][1])/
                                                     (lpath2.vertices[0][0] - lpath2.vertices[-20][0]))
                                            x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                            y_prime = -slope*x_prime
                                            plt.text(x_prime+lpath2.vertices[-20][0], 
                                                     y_prime+lpath2.vertices[-20][1], 
                                                     int(np.abs(self._Var.stoch[self._Var.stoch_row.index(j[k][1])][i])), 
                                                     fontsize=self.fontsize, 
                                                     horizontalalignment='center', 
                                                     verticalalignment='center', 
                                                     color=self.reactionColor)
                                    
                                else: # Uni-Uni
                                    XY = np.vstack((X1, X2, X3))
                                    
                                    tck, u = interpolate.splprep([XY[:,0], XY[:,1]], k=2)
                                    intX, intY = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
                                    stackXY = np.vstack((intX, intY))
                                    
                                    X3top = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y()+p3.get_height())
                                    X3bot = (p3.get_x()+p3.get_width()/2,
                                             p3.get_y())
                                    X3left = (p3.get_x(),
                                              p3.get_y()+p3.get_height()/2)
                                    X3right = (p3.get_x()+p3.get_width(),
                                               p3.get_y()+p3.get_height()/2)
                                    
                                    n_1 = -1
                                    arrthres_h = .02
                                    arrthres_v = .02
                                    while (((stackXY.T[n_1][0] > (X3left[0]-arrthres_h)) and 
                                            (stackXY.T[n_1][0] < (X3right[0]+arrthres_h)) and
                                            (stackXY.T[n_1][1] > (X3bot[1]-arrthres_v)) and 
                                            (stackXY.T[n_1][1] < (X3top[1]+arrthres_v))) and
                                            (np.abs(n_1) < np.shape(stackXY)[1] - 75)):
                                        n_1 -= 1
                                   
                                    if iVar.r_type[i] == 'reversible':
                                        X1top = (p1.get_x()+p1.get_width()/2,
                                                 p1.get_y()+p1.get_height())
                                        X1bot = (p1.get_x()+p1.get_width()/2,
                                                 p1.get_y())
                                        X1left = (p1.get_x(),
                                                  p1.get_y()+p1.get_height()/2)
                                        X1right = (p1.get_x()+p1.get_width(),
                                                   p1.get_y()+p1.get_height()/2)
                                        
                                        n_2 = 0
                                        while (((stackXY.T[n_2][0] > (X1left[0]-arrthres_h)) and 
                                                (stackXY.T[n_2][0] < (X1right[0]+arrthres_h)) and
                                                (stackXY.T[n_2][1] > (X1bot[1]-arrthres_v)) and 
                                                (stackXY.T[n_2][1] < (X1top[1]+arrthres_v))) and
                                                (np.abs(n_2) < np.shape(stackXY)[1] - 75)):
                                            n_2 += 1
                                        
                                        lpath = Path(stackXY.T[n_2:n_1])
                                        
                                        if self.analyzeFlux:
                                            if iVar.flux[i] > 0:
                                                lw1 = (1+self.edgelw)
                                                lw2 = (4+self.edgelw)
                                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=1.2, head_width=0.8)
                                            elif iVar.flux[i] < 0:
                                                lw1 = (4+self.edgelw)
                                                lw2 = (1+self.edgelw)
                                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=1.2, head_width=0.8)
                                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                            else:
                                                lw1 = (1+self.edgelw)
                                                lw2 = (1+self.edgelw)
                                                arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                                arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                                
                                            e1 = FancyArrowPatch(path=Path(stackXY.T[-n_1:50]),
                                                                arrowstyle=arrowstyle1,
                                                                mutation_scale=10.0,
                                                                lw=lw1,
                                                                color=self.reactionColor)
                                            e2 = FancyArrowPatch(path=Path(stackXY.T[50:n_1]),
                                                                arrowstyle=arrowstyle2,
                                                                mutation_scale=10.0,
                                                                lw=lw2,
                                                                color=self.reactionColor)
                                            fig.axes[mdl].add_patch(e1)
                                            fig.axes[mdl].add_patch(e2)
                                        else:
                                            arrowstyle = ArrowStyle.CurveFilledAB(head_length=0.8, head_width=0.4)
                                            e = FancyArrowPatch(path=lpath,
                                                            arrowstyle=arrowstyle,
                                                            mutation_scale=10.0,
                                                            lw=(1+self.edgelw),
                                                            color=self.reactionColor)
                                            fig.axes[mdl].add_patch(e)
                                        
                                    else:
                                        lpath = Path(stackXY.T[:n_1])
                                        arrowstyle = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                        e = FancyArrowPatch(path=lpath,
                                                            arrowstyle=arrowstyle,
                                                            mutation_scale=10.0,
                                                            lw=(1+self.edgelw),
                                                            color=self.reactionColor)
                                        fig.axes[mdl].add_patch(e)
                                
                                    if j[k][0] in iVar.floatingId:
                                        if (np.abs(iVar.stoch[iVar.stoch_row.index(j[k][0])][i]) > 1):
                                            slope = ((lpath.vertices[0][1] - lpath.vertices[10][1])/
                                                     (lpath.vertices[0][0] - lpath.vertices[10][0]))
                                            x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                            y_prime = -slope*x_prime
                                            plt.text(x_prime+lpath.vertices[10][0], 
                                                     y_prime+lpath.vertices[10][1], 
                                                     int(np.abs(iVar.stoch[iVar.stoch_row.index(j[k][0])][i])), 
                                                     fontsize=self.fontsize, 
                                                     horizontalalignment='center', 
                                                     verticalalignment='center', 
                                                     color=self.reactionColor)
                                    
                                    if j[k][1] in iVar.floatingId:
                                        if (np.abs(iVar.stoch[iVar.stoch_row.index(j[k][1])][i]) > 1):
                                            slope = ((lpath.vertices[0][1] - lpath.vertices[-20][1])/
                                                     (lpath.vertices[0][0] - lpath.vertices[-20][0]))
                                            x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                            y_prime = -slope*x_prime
                                            plt.text(x_prime+lpath.vertices[-20][0], 
                                                     y_prime+lpath.vertices[-20][1],
                                                     int(np.abs(iVar.stoch[iVar.stoch_row.index(j[k][1])][i])), 
                                                     fontsize=self.fontsize, 
                                                     horizontalalignment='center', 
                                                     verticalalignment='center',
                                                     color=self.reactionColor)
                            
                    else: # BIBI or larger
                        if len(iVar.rct[i]) < len(iVar.prd[i]):
                            rVal = len(iVar.rct[i])
                        else:
                            rVal = len(iVar.prd[i])
                            
                        for j in [list(zip(x,iVar.prd[i])) for x in itertools.combinations(iVar.rct[i],rVal)][0]:
                            p1 = iVar.G.nodes[j[0]]['patch']
                            p2 = iVar.G.nodes[iVar.rid[i]]['patch']
                            p3 = iVar.G.nodes[j[1]]['patch']
                            
                            X1 = (p1.get_x()+p1.get_width()/2,p1.get_y()+p1.get_height()/2)
                            X2 = (p2.get_x()+p2.get_width()/2,p2.get_y()+p2.get_height()/2)
                            X3 = (p3.get_x()+p3.get_width()/2,p3.get_y()+p3.get_height()/2)
                            
                            XY = np.vstack((X1, X2, X3))
                            
                            tck, u = interpolate.splprep([XY[:,0], XY[:,1]], k=2)
                            intX, intY = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
                            stackXY = np.vstack((intX, intY))
                            
                            X3top = (p3.get_x()+p3.get_width()/2,
                                     p3.get_y()+p3.get_height())
                            X3bot = (p3.get_x()+p3.get_width()/2,
                                     p3.get_y())
                            X3left = (p3.get_x(),
                                      p3.get_y()+p3.get_height()/2)
                            X3right = (p3.get_x()+p3.get_width(),
                                       p3.get_y()+p3.get_height()/2)
                            
                            n_1 = -1
                            arrthres_h = .02
                            arrthres_v = .02
                            while (((stackXY.T[n_1][0] > (X3left[0]-arrthres_h)) and
                                    (stackXY.T[n_1][0] < (X3right[0]+arrthres_h)) and
                                    (stackXY.T[n_1][1] > (X3bot[1]-arrthres_v)) and 
                                    (stackXY.T[n_1][1] < (X3top[1]+arrthres_v))) and
                                    (np.abs(n_1) < np.shape(stackXY)[1] - 75)):
                                n_1 -= 1
                            
                            if iVar.r_type[i] == 'reversible':
                                X1top = (p1.get_x()+p1.get_width()/2,
                                         p1.get_y()+p1.get_height())
                                X1bot = (p1.get_x()+p1.get_width()/2,
                                         p1.get_y())
                                X1left = (p1.get_x(),
                                          p1.get_y()+p1.get_height()/2)
                                X1right = (p1.get_x()+p1.get_width(),
                                           p1.get_y()+p1.get_height()/2)
                                
                                n_2 = 0
                                while (((stackXY.T[n_2][0] > (X1left[0]-arrthres_h)) and 
                                        (stackXY.T[n_2][0] < (X1right[0]+arrthres_h)) and
                                        (stackXY.T[n_2][1] > (X1bot[1]-arrthres_v)) and 
                                        (stackXY.T[n_2][1] < (X1top[1]+arrthres_v))) and
                                        (np.abs(n_2) < np.shape(stackXY)[1] - 75)):
                                    n_2 += 1
                                
                                lpath = Path(stackXY.T[n_2:n_1])
                                
                                if self.analyzeFlux:
                                    if iVar.flux[i] > 0:
                                        lw1 = (1+self.edgelw)
                                        lw2 = (4+self.edgelw)
                                        arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                        arrowstyle2 = ArrowStyle.CurveFilledB(head_length=1.2, head_width=0.8)
                                    elif iVar.flux[i] < 0:
                                        lw1 = (4+self.edgelw)
                                        lw2 = (1+self.edgelw)
                                        arrowstyle1 = ArrowStyle.CurveFilledA(head_length=1.2, head_width=0.8)
                                        arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                    else:
                                        lw1 = (1+self.edgelw)
                                        lw2 = (1+self.edgelw)
                                        arrowstyle1 = ArrowStyle.CurveFilledA(head_length=0.8, head_width=0.4)
                                        arrowstyle2 = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                        
                                    e1 = FancyArrowPatch(path=Path(stackXY.T[n_2:50]),
                                                        arrowstyle=arrowstyle1,
                                                        mutation_scale=10.0,
                                                        lw=lw1,
                                                        color=self.reactionColor)
                                    e2 = FancyArrowPatch(path=Path(stackXY.T[50:n_1]),
                                                        arrowstyle=arrowstyle2,
                                                        mutation_scale=10.0,
                                                        lw=lw2,
                                                        color=self.reactionColor)
                                    fig.axes[mdl].add_patch(e1)
                                    fig.axes[mdl].add_patch(e2)
                                else:
                                    arrowstyle = ArrowStyle.CurveFilledAB(head_length=0.8, head_width=0.4)
                                    e = FancyArrowPatch(path=lpath,
                                                    arrowstyle=arrowstyle,
                                                    mutation_scale=10.0,
                                                    lw=(1+self.edgelw),
                                                    color=self.reactionColor)
                                    fig.axes[mdl].add_patch(e)
                                
                            else:
                                lpath = Path(stackXY.T[:n_1])
                                arrowstyle = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                                e = FancyArrowPatch(path=lpath,
                                                    arrowstyle=arrowstyle,
                                                    mutation_scale=10.0,
                                                    lw=(1+self.edgelw),
                                                    color=self.reactionColor)
                                fig.axes[mdl].add_patch(e)
                            
                            if j[0] in iVar.floatingId:
                                if (np.abs(iVar.stoch[iVar.stoch_row.index(j[0])][i]) > 1):
                                    slope = ((lpath.vertices[0][1] - lpath.vertices[15][1])/
                                             (lpath.vertices[0][0] - lpath.vertices[15][0]))
                                    x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                    y_prime = -slope*x_prime
                                    plt.text(x_prime+lpath.vertices[15][0], 
                                             y_prime+lpath.vertices[15][1], 
                                             int(np.abs(iVar.stoch[iVar.stoch_row.index(j[0])][i])), 
                                             fontsize=self.fontsize, 
                                             horizontalalignment='center', 
                                             verticalalignment='center', 
                                             color=self.reactionColor)
                            if j[1] in iVar.floatingId:
                                if (np.abs(iVar.stoch[iVar.stoch_row.index(j[1])][i]) > 1):
                                    slope = ((lpath.vertices[0][1] - lpath.vertices[-20][1])/
                                             (lpath.vertices[0][0] - lpath.vertices[-20][0]))
                                    x_prime = np.sqrt(0.01/(1 + np.square(slope)))*(self.fontsize/20)*max(self.scale/2, 1)
                                    y_prime = -slope*x_prime
                                    plt.text(x_prime+lpath.vertices[-20][0], 
                                             y_prime+lpath.vertices[-20][1], 
                                             int(np.abs(iVar.stoch[iVar.stoch_row.index(j[1])][i])), 
                                             fontsize=self.fontsize,
                                             horizontalalignment='center', 
                                             verticalalignment='center', 
                                             color=self.reactionColor)
                            
                # Modifiers
                seen={}
                for i, e in enumerate(iVar.mod_flat):
                    n1 = iVar.G.nodes[e]['patch']
                    n2 = iVar.G.nodes[iVar.modtarget_flat[i]]['patch']
                    rad = 0.1
                    shrinkB = 2.
                    
                    if (e,iVar.modtarget_flat[i]) in seen:
                        rad = seen.get((e,iVar.modtarget_flat[i])) # TODO: No curvature when there is just a single line between two nodes
                        rad = (rad+np.sign(rad)*0.1)*-1 # TODO: Change curvature
                        
                    X1 = (n1.get_x()+n1.get_width()/2,
                          n1.get_y()+n1.get_height()/2)
                    X2 = (n2.get_x()+n2.get_width()/2,
                          n2.get_y()+n2.get_height()/2)
                    
                    if iVar.modtype_flat[i] == 'inhibitor': # inhibition
                        color = self.modifierColor
                        arrowstyle = ArrowStyle.BarAB(widthA=0.0, angleA=None, widthB=1.0, angleB=None)
                        shrinkB = 10.
                        linestyle = '-'
                    elif iVar.modtype_flat[i] == 'activator': # activation
                        color = self.modifierColor
                        arrowstyle = arrowstyle = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                        linestyle = '-'
                    elif iVar.modtype_flat[i] == 'modifier': # Unknown modifier
                        color = self.modifierColor
                        arrowstyle = arrowstyle = ArrowStyle.CurveFilledB(head_length=0.8, head_width=0.4)
                        linestyle = ':'
                    e = FancyArrowPatch(X1,
                                        X2,
                                        patchA=n1,
                                        patchB=n2,
                                        shrinkB=shrinkB,
                                        arrowstyle=arrowstyle,
                                        connectionstyle='arc3,rad=%s'%rad,
                                        mutation_scale=10.0,
                                        lw=(1+self.edgelw),
                                        color=color,
                                        linestyle=linestyle)
                    seen[(e,iVar.modtarget_flat[i])]=rad
                    fig.axes[mdl].add_patch(e)
                    fig.axes[mdl].add_patch(n1)
                
                # Add reaction nodes at last to put it on top
                if self.drawReactionNode:
                    allnodes = iVar.speciesId + iVar.rid
                else:
                    allnodes = iVar.speciesId
                
                if 'Input' in iVar.G.nodes:
                    allnodes += ['Input']
                if 'Output' in iVar.G.nodes:
                    allnodes += ['Output']
                for i in range(len(allnodes)):
                    fig.axes[mdl].add_patch(iVar.G.nodes[allnodes[i]]['patch'])
                
                fig.axes[mdl].autoscale()
        
        if savePath != None:
            try:
                fig.savefig(savePath, bbox_inches='tight', dpi=dpi)
            except IOError as e:
                raise Exception("Error saving diagram: " + str(e))
        else:
            if show and self.customAxis == None:
                plt.show()
            plt.close()
        
        self.edgelw = edgelw_backup


    def savefig(self, path, dpi=150):
        """
        Save weighted network diagram to specified location
        
        :param path: path to save the diagram
        :param dpi: dpi settings for the diagram
        """
        
        allRxn, count, fig = self.drawWeightedDiagram(savePath=path, dpi=dpi)
        
