#!/usr/bin/env python
import gtk
import gobject

import matplotlib.pyplot as plt   

from matplotlib.figure import Figure
# uncomment to select /GTK/GTKAgg/GTKCairo
# from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
# from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import sys, os
cur = os.path.dirname(os.path.realpath(__file__))

sys.path.append(cur + '/../UtilityLib')
from moduleXML import loadNetworkFromXML 
from processing import linearWaveSplitting
from processing import nonLinearWaveSplitting
from processing import minMaxFunction

from moduleStartUp import parseOptions

from modulePickle import loadSolutionDataFile
from modulePickle import loadExternalDataSet


sys.path.append(cur + '/../NetworkLib')
from classVascularNetwork import VascularNetwork 

from matplotlib.font_manager import FontProperties

import itertools
import subprocess

from optparse import OptionParser
import cPickle
import numpy as np

from pprint import pprint as pp

from numpy import arange, sin, pi

from copy import deepcopy as copy

class Visualisation2DPlotWindowAdjustValues(gtk.Window):
    
    def __init__(self, parent, variableName, actualValuesDict, inititalValuesDict):
        super(Visualisation2DPlotWindowAdjustValues, self).__init__()
        
        self.set_resizable(False)
        #self.set_position(gtk.WIN_POS_CENTER)
                                
        self.parentWindow       = parent
        self.variableName       = variableName
        self.inititalValuesDict = copy(inititalValuesDict)
        self.actualValuesDict   = actualValuesDict
        
        vbox = gtk.VBox(False, 1)
               
        buttonUpdate = gtk.Button("update")
        buttonUpdate.connect("clicked", self.on_clickedUpdate)
        buttonUpdate.set_size_request(120,30)
                
        buttonReset = gtk.Button("reset")
        buttonReset.connect("clicked", self.on_clickedReset)
        buttonReset.set_size_request(120,30)
        
        hbox = gtk.HBox(False, 10) 
        hbox.pack_start(buttonUpdate, fill=False, expand=True)  
        hbox.pack_start(buttonReset, fill=False, expand=True)   
        vbox.pack_start(hbox, expand=True, fill=False)
                
        self.entries = {}
        
        countX = 1
        countY = 1
        
        for key,values in actualValuesDict.iteritems():
            hbox = gtk.HBox(False, 10) 
            label = gtk.Label(key)
            hbox.pack_start(label, fill=False, expand=True)     
            self.entries[key] = []
            
            countX = countX +1
            
            for value in values:
                entry = gtk.Entry()
                entry.add_events(gtk.gdk.KEY_RELEASE_MASK)
                entry.set_text(str(value))
                hbox.pack_start(entry, fill=False, expand=False)  
                self.entries[key].append(entry)   
                
                countY = countY+1
            vbox.pack_start(hbox, expand=True, fill=False)
            
        self.set_size_request(30*countX, 15*countY)
            
        self.connect("destroy", self.on_destroy)   
        
        self.add(vbox)
        self.show_all()
    
    def on_destroy(self, widget):
        self.parentWindow.update({'limitsWindow':None})
    
    def on_clickedUpdate(self, widget):
        for key, entries in self.entries.iteritems():
            for entry in entries:
                text = entry.get_text()
                index = entries.index(entry)
                if text is not '':
                    try: self.actualValuesDict[key][index] = float(text)
                    except: print "WARNING: insert float please and not :",text
        self.parentWindow.update({self.variableName:self.actualValuesDict})
        self.parentWindow.plot()

    def on_clickedReset(self, widget):
        for key, entries in self.entries.iteritems():
            for entry,i in zip(entries,xrange(len(entries))):
                value = self.inititalValuesDict[key][i]
                entry.set_text(str(value))
        self.parentWindow.update({self.variableName: copy(self.inititalValuesDict)})
        self.parentWindow.plot()

class Visualisation2DPlotWindowGui(gtk.Window):
    def __init__(self):
        # self.set_title(' '.join(['2d Visualisation - ',self.vascularNetwork.name,' : ',self.vessel.name,'- dataNumber ',dataNumber]))
        super(Visualisation2DPlotWindowGui, self).__init__()
        
        self.axisX = 'Time'
        self.scale = gtk.HScale()
        self.scale.set_range(0, 10)
        self.scale.set_increments(1, 1)
        self.scale.set_digits(0)
        self.scale.set_size_request(400, 30)
        self.scale.connect('value-changed', self.on_changedSlider)
        self.scale.set_draw_value(False)        
        
        self.nodeLabel = gtk.Label('Node 0')
        
        # # image legend
        self.legend = False
        cBlegend = gtk.CheckButton("show legend")
        # cBlegend.connect("toggled", self.on_changedLegendBool)
        cBlegend.set_size_request(120, 30)
        
        # plot title
        self.plotTitle = False
        cBtitle = gtk.CheckButton("show title")
        # cBtitle.connect("toggled", self.on_changedTitleBool)
        cBtitle.set_size_request(120, 30)
        
        # medical units
#         self.medicalUnit = True
#         cBmedical = gtk.CheckButton("use medical units")
#         #cBmedical.connect("toggled", self.on_changedMedicalBool)
#         cBmedical.set_size_request(120,30)
#         cBmedical.set_active(1) 

        # show decriptions
        self.showDescription = False
        cBdescription = gtk.CheckButton('show descriptions')
        # cBdescription.connect("toggled", self.on_changedDescription)
        cBdescription.set_size_request(120, 30)
        # render movie button
        self.buttonRenderMovie = gtk.Button("render movie")
        # self.buttonRenderMovie.connect("clicked", self.on_clickedRenderMovie)
        self.buttonRenderMovie.set_size_request(120, 30)
        self.buttonRenderMovie.set_sensitive(False)

        #### second row of buttons
        # min max plots
        self.plotMinMaxPoints = False
        cBplotMinMaxPoints = gtk.CheckButton('plot MinMax-points')
        # cBplotMinMaxPoints.connect("toggled", self.on_changedPlotMinMax)
        cBplotMinMaxPoints.set_size_request(120, 30)
        
#         cBchangeMinMaxDelta = gtk.Button('change MinMax-deltas')
#         cBchangeMinMaxDelta.connect("clicked", self.on_clickedMinMaxDelta)
#         cBchangeMinMaxDelta.set_size_request(120,30)
#         cBchangeMinMaxDelta.set_sensitive(False)

        buttonLimits = gtk.Button("update limits")
        buttonLimits.connect("clicked", self.on_changeLimits)
        buttonLimits.set_size_request(120, 30)

        
        self.fig = Figure(figsize=(5, 4), dpi=100)
        # self.plot(0)
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.canvas.set_size_request(640, 690)
        
        toolbar = NavigationToolbar(self.canvas, self)
        
        cBtitle.set_active(0) 
        cBdescription.set_active(0) 
        
        cbType = gtk.combo_box_new_text()
        cbType.append_text('Plot P,Q')
        cbType.append_text('Plot P,Q with linear wavesplitting')
        cbType.append_text('Plot P,Q with non-linear wavesplitting')
        cbType.append_text('Plot CFL, wave speed')
        cbType.append_text('Plot Area, Compliance')
        cbType.set_active(0) 
        cbType.connect("changed", self.on_changePlotType)
        
        self.cbXaxis = gtk.combo_box_new_text()
        self.cbXaxis.connect("changed", self.on_changePlotXaxis)
        self.cbXaxis.append_text('Plot over Time')
        self.cbXaxis.append_text('Plot over Space')
        self.cbXaxis.set_active(0) 
        self.cbXaxis.set_sensitive(False)

        vbox = gtk.VBox(False, 1)
        hboxLegend = gtk.HBox(False, 10) 
        hboxCheckboxes2 = gtk.HBox(False, 10) 
        hbox = gtk.HBox(False, 1)
        hbox2 = gtk.HBox(False, 1)
                
        # Checkbutton series
        hboxLegend.pack_start(cBlegend)
        hboxLegend.pack_start(cBdescription)
        hboxLegend.pack_start(buttonLimits)
        # hboxLegend.pack_start(cBmedical)
        hboxLegend.pack_start(self.buttonRenderMovie)
        hboxCheckboxes2.pack_start(cBplotMinMaxPoints)
        # hboxCheckboxes2.pack_start(cBchangeMinMaxDelta)

        # align pictures canvas
        alignIm = gtk.Alignment(0, 1 , 1, 0)
        alignIm.add(self.canvas)
        # align node switcher scale
        hbox.pack_start(self.nodeLabel)
        hbox.pack_start(self.scale)
        alignHbox = gtk.Alignment(0, 1, 1, 0)
        alignHbox.add(hbox)
        # align combobox
        hbox2.pack_start(cbType)
        hbox2.pack_start(self.cbXaxis)
        alignCB = gtk.Alignment(0, 1, 1, 0)  
        alignCB.add(hbox2)
        # align navigation toolbox
        alignNT = gtk.Alignment(0, 1, 1, 0)  
        alignNT.add(toolbar)
        
        # put all together
        vbox.pack_start(hboxLegend)
        vbox.pack_start(hboxCheckboxes2)
        vbox.pack_start(alignCB)
        vbox.pack_start(alignNT)
        vbox.pack_start(alignIm)
        vbox.pack_start(alignHbox)
        
        self.add(vbox)         
        self.show_all()
            
    def on_changePlotXaxis(self, widget):
        pass
        
    def on_changedSlider(self, widget):
        pass
            
    def on_changeLimits(self, widget):
        pass
            
class Visualisation2DPlotWindow(Visualisation2DPlotWindowGui):
    def __init__(self, selectedNetworks, selectedVesselIds, selectedExternalData):
        
        self.plot = lambda :''
        
        self.limitsWindow = None
        
        # # variables for the slider
        self.sliderValue = 0  # grid node / time point
        self.limits = {'Time':      [0, 10],
                       'Space':     [0, 10],
                       'gridNodes': [0, 5]}
        
        
        # # initialize super class
        super(Visualisation2DPlotWindow, self).__init__()
        
        self.selectedNetworks = selectedNetworks
        self.selectedVesselIds = selectedVesselIds
        self.selectedExternalData = selectedExternalData
        
        # # activate space/time changer if only one network is loaded
        if len(self.selectedNetworks) == 1:
            self.cbXaxis.set_sensitive(True)
                
        self.linewidth = 1.5  
        self.fontSizeLabel = 14
        
        self.unitPtext = '$mmHg$'
        self.unitFtext = '$ml/s$'
        
        self.createGraph()
        self.estimatePlotLimits()
        self.plot = self.updateLinesPQ
        
        # # update slider and so --> update the plot
        self.scale.set_range(*self.limits['gridNodes'])
        self.plot()
    
    def update(self, visualisationData):
        '''
        updates the vascularNetwork data using a dictionary in form of 
        visualisationData = {'variableName': value}
        '''
        for key,value in visualisationData.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except: 
                print 'WARNING vascularNetwork.update(): wrong key: %s, could not update vascularNetwork' %key 
    
       
    def createGraph(self):
        '''
        create graph with 2 subplots and all necessary lines 
        '''
        self.fig = plt.figure(figsize=(6, 4), dpi=100, edgecolor='k')
        self.fig.subplots_adjust(right=0.86)
        self.fig.subplots_adjust(left=0.17)
        self.fig.subplots_adjust(top=0.95)
        self.fig.subplots_adjust(bottom=0.15)
        self.fig.subplots_adjust(hspace=0.18)
        
        fontLegend = FontProperties()
        fontLegend.set_size(self.fontSizeLabel)
        
        from matplotlib import rc
        from matplotlib import rcParams
        
        rcParams['text.usetex'] = True
        rcParams['text.latex.unicode'] = True
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.size'] = self.fontSizeLabel
        rcParams['savefig.dpi'] = 300.
        
        
        ax1 = plt.subplot(2, 1, 1, frameon=True)
        
        ax12 = ax1.twinx()
        
        ax12.set_visible(False)
        
        plt.xticks(np.linspace(ax1.get_xlim()[0], ax1.get_xlim()[1], 2), ['', '']) 
        ax1.tick_params(axis='x', top='off', bottom='off')
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.tick_params(axis='y', right='off')
        
        #plt.yticks(np.linspace(ax12.get_xlim()[0], ax12.get_xlim()[1], 2), ['', '']) 
                
        ax2 = plt.subplot(2, 1, 2, frameon=True)
        ax2.spines['top'].set_visible(False)
        ax2.tick_params(axis='x', top='off')
        ax2.spines['right'].set_visible(False)
        ax2.tick_params(axis='y', right='off')
           
        ax2.set_xlabel('Time $s$}', fontsize=self.fontSizeLabel)
           
        ax22 = ax2.twinx()
        ax22.set_visible(False)
        
        #plt.yticks(np.linspace(ax22.get_xlim()[0], ax22.get_xlim()[1], 2), ['', '']) 
                           
        self.axis = [ax1, ax12 ,ax2, ax22]
        self.lines = [[], []]
        # # add 3 lines for each network case to each subplot
        colors = ['b', 'r', 'm', 'g', 'c']
        for i in xrange(len(self.selectedVesselIds)):
            
            self.lines[0].append([self.axis[0].plot(-1, 0, color=colors[i], linestyle='-', linewidth=self.linewidth)[0]])
            self.lines[0][i].append(self.axis[1].plot(-1, 0, color=colors[i], linestyle='--', linewidth=self.linewidth)[0])
            self.lines[0][i].append(self.axis[1].plot(-1, 0, color=colors[i], linestyle=':', linewidth=self.linewidth)[0])
            
            self.lines[1].append([self.axis[2].plot(-1, 0, color=colors[i], linestyle='-', linewidth=self.linewidth)[0]])
            self.lines[1][i].append(self.axis[3].plot(-1, 0, color=colors[i], linestyle='--', linewidth=self.linewidth)[0])
            self.lines[1][i].append(self.axis[3].plot(-1, 0, color=colors[i], linestyle=':', linewidth=self.linewidth)[0])   
        
    def clearLines(self):
        for lineAxisPack in self.lines:
            for linePack in lineAxisPack:
                for line in linePack:
                    line.set_data(-1, 0)
                    
        self.axis[1].set_visible(False)
        self.axis[3].set_visible(False)   
        self.updatePlotWindow()
        
        
        
    def updatePlotWindow(self):
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()        
                
    def updateLinesPQ(self):
        
        gridNode = self.sliderValue
        # 1. set axis lable
        self.axis[0].set_ylabel('Pressure ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis[2].set_ylabel('Flow ' + self.unitFtext, fontsize=self.fontSizeLabel)
        # 2. update lines for P and Q over time for grid node 0
                
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                try:
                    yData00 = vascularNetwork.vessels[vesselId].Psol[:, [gridNode]] / 133.32
                    ydata10 = vascularNetwork.vessels[vesselId].Qsol[:, [gridNode]] * 1e6    
                    
                    xData = vascularNetwork.simulationTime
                    
                    self.lines[0][i][0].set_data(xData, yData00)
                    self.lines[1][i][0].set_data(xData, ydata10)
                    
                    self.axis[2].set_xlabel('Time $s$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Time'])
                    self.axis[2].set_xlim(self.limits['Time'])
                except:
                    self.lines[0][i][0].set_data(-1, 0)
                    self.lines[1][i][0].set_data(-1, 0)
            
            elif self.axisX == "Space":
                try: 
                    yData00 = vascularNetwork.vessels[vesselId].Psol[gridNode] / 133.32
                    ydata10 = vascularNetwork.vessels[vesselId].Qsol[gridNode] * 1e6    
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    
                    self.lines[0][i][0].set_data(xData, yData00)
                    self.lines[1][i][0].set_data(xData, ydata10)      
                    
                    self.axis[1].set_xlabel('Space $cm$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Space'])
                    self.axis[1].set_xlim(self.limits['Space'])              
                except:
                    self.lines[0][i][0].set_data(-1, 0)
                    self.lines[1][i][0].set_data(-1, 0)
        
        self.axis[0].set_ylim(self.limits['P'])
        self.axis[2].set_ylim(self.limits['Q'])
        
        self.updatePlotWindow()
          
    def updateLinesPQsplitLin(self):
        gridNode = self.sliderValue
        # 1. set axis lable
        self.axis[0].set_ylabel('Pressure ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis[1].set_ylabel('Contribution ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis[2].set_ylabel('Flow ' + self.unitFtext, fontsize=self.fontSizeLabel)
        self.axis[3].set_ylabel('Contribution ' + self.unitFtext, fontsize=self.fontSizeLabel)
        # 2. update lines for P and Q over time for grid node 0
        
        self.axis[1].set_visible(True)
        self.axis[3].set_visible(True) 
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                 try:
                    Psol = vascularNetwork.vessels[vesselId].Psol[:, [gridNode]]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[:, [gridNode]]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[:, [gridNode]]   
                    csol = vascularNetwork.vessels[vesselId].csol[:, [gridNode]]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = linearWaveSplitting(Psol, Qsol, Asol, csol, vascularNetwork.vessels[vesselId].rho)
                    
                    yData00 = Psol / 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32
                    yData10 = Qsol * 1.e6   
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    xData = vascularNetwork.simulationTime
                    
                    self.lines[0][i][0].set_data(xData,     yData00)
                    self.lines[0][i][1].set_data(xData[1:], yData01)
                    self.lines[0][i][2].set_data(xData[1:], yData02)
                    
                    self.lines[1][i][0].set_data(xData,     yData10)
                    self.lines[1][i][1].set_data(xData[1:], yData11)
                    self.lines[1][i][2].set_data(xData[1:], yData12)
                                        
                    self.axis[2].set_xlabel('Time $s$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Time'])
                    self.axis[2].set_xlim(self.limits['Time'])
                 except:
                    self.lines[0][i][0].set_data(-1,0)
                    self.lines[1][i][0].set_data(-1,0)
            
            elif self.axisX == "Space":
                try:                     
                    Psol = vascularNetwork.vessels[vesselId].Psol[gridNode]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[gridNode]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[gridNode]   
                    csol = vascularNetwork.vessels[vesselId].csol[gridNode]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = linearWaveSplitting(Psol, Qsol, Asol, csol, vascularNetwork.vessels[vesselId].rho)
                    
                    yData00 = Psol / 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32
                    yData10 = Qsol * 1.e6   
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    print len(xData)
                    
                    self.lines[0][i][0].set_data(xData,     yData00)
                    self.lines[0][i][1].set_data(xData[1:], yData01)
                    self.lines[0][i][2].set_data(xData[1:], yData02)
                    
                    self.lines[1][i][0].set_data(xData,     yData10)
                    self.lines[1][i][1].set_data(xData[1:], yData11)
                    self.lines[1][i][2].set_data(xData[1:], yData12)
                                        
                    self.axis[2].set_xlabel('Space $cm$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Space'])
                    self.axis[2].set_xlim(self.limits['Space']) 
                                 
                except:
                    self.lines[0][i][0].set_data(-1, 0)
                    self.lines[1][i][0].set_data(-1, 0)
                    
        
        self.axis[0].set_ylim(self.limits['P'])
        self.axis[2].set_ylim(self.limits['Q'])
        
        self.axis[1].set_ylim(self.limits['Pfb'])
        self.axis[3].set_ylim(self.limits['Qfb'])
        
        
        self.updatePlotWindow()
    
    def updateLinesPQspliNonLin(self):
        pass
    
    def updateLinesWaveCFL(self):
        gridNode = self.sliderValue
        # 1. set axis lable
        self.axis[0].set_ylabel('Wave Speed $m/s$', fontsize=self.fontSizeLabel)
        self.axis[2].set_ylabel('CFL $-$' + self.unitFtext, fontsize=self.fontSizeLabel)
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                try:
                    # wave speed
                    yData00 = vascularNetwork.vessels[vesselId].csol[:, [gridNode]]
                    # CFL 
                    dz = vascularNetwork.vessels[vesselId].dz
                    ydata10 = yData00 * vascularNetwork.dt / sum(dz) * len(dz)
                    
                    xData = vascularNetwork.simulationTime
                    
                    self.lines[0][i][0].set_data(xData, yData00)
                    self.lines[1][i][0].set_data(xData, ydata10)
                    
                    self.axis[2].set_xlabel('Time $s$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Time'])
                    self.axis[2].set_xlim(self.limits['Time'])
                except:
                    self.lines[0][i][0].set_data(-1, 0)
                    self.lines[1][i][0].set_data(-1, 0)
            
            elif self.axisX == "Space":
                try: 
                    # wave speed
                    yData00 = vascularNetwork.vessels[vesselId].csol[gridNode]
                    # CFL 
                    dz = vascularNetwork.vessels[vesselId].dz
                    ydata10 = yData00 * vascularNetwork.dt / sum(dz) * len(dz)
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    print np.shape(xData)
                    
                    print yData00, xData
                    
                    self.lines[0][i][0].set_data(xData, yData00)
                    self.lines[1][i][0].set_data(xData, ydata10)      
                    
                    self.axis[2].set_xlabel('Space $cm$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Space'])
                    self.axis[2].set_xlim(self.limits['Space'])     
                except:
                    self.lines[0][i][0].set_data(-1, 0)
                    self.lines[1][i][0].set_data(-1, 0)
        
        self.axis[0].set_ylim(self.limits['c'])
        self.axis[2].set_ylim(self.limits['CFL'])
        
        self.updatePlotWindow()
    
    def updateLinesAreaComp(self):
        gridNode = self.sliderValue
        # 1. set axis lable
        self.axis[0].set_ylabel('Area $mm^2$', fontsize=self.fontSizeLabel)
        self.axis[2].set_ylabel('Compliance $mm^2/mmHg$', fontsize=self.fontSizeLabel)
        # 2. update lines for P and Q over time for grid node 0
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                try:
                    yData00 = vascularNetwork.vessels[vesselId].Asol[:, gridNode] * 1000 * 1000
                    Psol = vascularNetwork.vessels[vesselId].Psol   
                    ydata10 = vascularNetwork.vessels[vesselId].C(Psol)[:, gridNode] * 1000 * 1000 / 133.32 
                    
                    xData = vascularNetwork.simulationTime
                    
                    self.lines[0][i][0].set_data(xData, yData00)
                    self.lines[1][i][0].set_data(xData, ydata10)
                    
                    self.axis[2].set_xlabel('Time $s$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Time'])
                    self.axis[2].set_xlim(self.limits['Time'])
                except:
                    self.lines[0][i][0].set_data(-1, 0)
                    self.lines[1][i][0].set_data(-1, 0)
            
            elif self.axisX == "Space":
                try: 
                    yData00 = vascularNetwork.vessels[vesselId].Asol[gridNode] * 1000 * 1000
                    Psol = vascularNetwork.vessels[vesselId].Psol[gridNode]   
                    ydata10 = vascularNetwork.vessels[vesselId].C(Psol) * 1000 * 1000 / 133.32   
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    
                    self.lines[0][i][0].set_data(xData, yData00)
                    self.lines[1][i][0].set_data(xData, ydata10)      
                    
                    self.axis[2].set_xlabel('Space $cm$}', fontsize=self.fontSizeLabel)
                    self.axis[0].set_xlim(self.limits['Space'])
                    self.axis[2].set_xlim(self.limits['Space'])              
                except:
                    self.lines[0][i][0].set_data(-1, 0)
                    self.lines[1][i][0].set_data(-1, 0)
        
        self.axis[0].set_ylim(self.limits['A'])
        self.axis[2].set_ylim(self.limits['C'])
        
        self.updatePlotWindow()
        
    def updateGravity(self):
        '''
        creates a plot of gravity in the first sub plot
        '''
        
    def on_changeLimits(self, widget):
        '''
        create window to change plot limits
        '''
        if self.limitsWindow == None:
            self.limitsWindow = Visualisation2DPlotWindowAdjustValues(self,'limits',self.limits,self.limitsInit)    
        else:
            self.limitsWindow.destroy()
        self.plot()
                  
    def estimatePlotLimits(self):
        '''
        This function evaluates all limits for the plots
        It may take a while to calculate all
        '''     
        self.limits = {'P':  [1e50, -1e50],
                       'Pf': [1e50, -1e50],
                       'Pb': [1e50, -1e50],
                       'Pfb': [1e50, -1e50],
                       'Q':  [1e50, -1e50],
                       'Qf': [1e50, -1e50],
                       'Qb': [1e50, -1e50],
                       'Qfb': [1e50, -1e50],
                       'Time':   [0, -1e50],
                       'Space':  [0, -1e50],
                       'c':  [1e50, -1e50],
                       'CFL':[ 0, 1.1],
                       'A':  [1e50, -1e50],
                       'C':  [1e50, -1e50],
                       'gridNodes': [0, -1e50]}
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            # pressure
            limit = 'P'
            Psol = vascularNetwork.vessels[vesselId].Psol 
            sol = Psol / 133.32
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
            # flow
            limit = 'Q'
            Qsol = vascularNetwork.vessels[vesselId].Qsol 
            sol = Qsol * 1.e6
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]  
            # time
            limit = 'Time'
            sol = vascularNetwork.simulationTime
            self.limits[limit] = [0, max(self.limits[limit][1], np.max(sol))]
            # space
            limit = 'Space'
            sol = vascularNetwork.vessels[vesselId].length * 100.
            self.limits[limit] = [0, max(self.limits[limit][1], np.max(sol))]
            # wave speed
            limit = 'c'
            csol = vascularNetwork.vessels[vesselId].csol
            sol = csol
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
            # CFL
            limit = 'CFL'
            sol = csol * vascularNetwork.dt / vascularNetwork.vessels[vesselId].dz[0] 
            self.limits[limit] = [0, max([self.limits[limit][1], np.max(sol)])]
            # area
            limit = 'A'
            Asol = vascularNetwork.vessels[vesselId].Asol 
            sol = Asol * 1000 * 1000
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
            # compliance
            limit = 'C'
            sol = vascularNetwork.vessels[vesselId].C(Psol) * 1000 * 1000 / 133.32
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
            # gridNodes
            limit = 'gridNodes'
            sol = vascularNetwork.vessels[vesselId].N - 1
            self.limits[limit] = [0, max([self.limits[limit][1], np.max(sol)])]
            # pressure / flow,  forward backward
            for n in xrange(int(vascularNetwork.vessels[vesselId].N)):
                pf,pb,qf,qb =  linearWaveSplitting(Psol[:,n],Qsol[:,n],Asol[:,n],csol[:,n],vascularNetwork.vessels[vesselId].rho)
                
                limit = 'Pf'
                sol = pf/133.32
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
                
                limit = 'Pb'
                sol = pb/133.32
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
                
                limit = 'Qf'
                sol = qf*1.e6
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
                
                limit = 'Qb'
                sol = qb*1.e6
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
                
                limit = 'Pfb'
                self.limits[limit] = [min([self.limits[limit][0], np.min(pf/133.32), np.min(pb/133.32)]), max([self.limits[limit][1], np.max(pf/133.32), np.max(pb/133.32)])]
                
                limit = 'Qfb'
                self.limits[limit] = [min([self.limits[limit][0], np.min(qf*1.e6), np.min(qb*1.e6)]), max([self.limits[limit][1], np.max(qf*1.e6), np.max(qb*1.e6)])]
                
        self.limitsInit = copy(self.limits)
            
            
            
        
    # callback function
    def on_changePlotXaxis(self, widget):   
        '''
        callback function for changed plot axis menu time/space
        
        where space plot is only possible if one network is loaded
        '''    
        
        cbIndex = widget.get_active()
        if cbIndex == 0:
            self.axisX = 'Time'
            self.nodeLabel.set_text('Node ' + str(self.sliderValue))
            self.buttonRenderMovie.set_sensitive(False)
            self.scale.set_range(*self.limits['gridNodes']) 
        
        elif cbIndex == 1:
            self.axisX = 'Space'
            NodeValue = str(self.selectedNetworks[0].simulationTime[self.sliderValue][0])
            self.nodeLabel.set_text(' '.join(['Time ', str(NodeValue)[:6]]))
            self.buttonRenderMovie.set_sensitive(True)
            self.scale.set_range(0, len(self.selectedNetworks[0].simulationTime) - 1)   
            
        self.plot()
        

    def on_changedSlider(self, widget):
        '''
        callback function for changed node/time slider
        '''
        self.sliderValue = int(widget.get_value())
        if self.axisX == 'Time':
            self.nodeLabel.set_text('Node ' + str(self.sliderValue))
        elif self.axisX == 'Space':
            NodeValue = str(self.selectedNetworks[0].simulationTime[self.sliderValue][0])
            self.nodeLabel.set_text(' '.join(['Time ', str(NodeValue)[:6]]))
        self.plot()
        
    def on_changePlotType(self, widget):
        '''
        call function of the combobox to choose type of plot
        '''
        self.clearLines()
        cbIndex = widget.get_active()
        if   cbIndex == 0:
            self.plot = self.updateLinesPQ
        elif cbIndex == 1:
            self.plot = self.updateLinesPQsplitLin
        elif cbIndex == 2:
            self.plot = self.updateLinesPQspliNonLin
        elif cbIndex == 3:
            self.plot = self.updateLinesWaveCFL
        elif cbIndex == 4:
            self.plot = self.updateLinesAreaComp
        
        self.plot()
        
class Visualisation2DMainCase(object):
    def __init__(self, number):
        
        self.networkInfo = {"choose simulation case" : [[], '-'] }
        self.networkCases = ["choose simulation case"]
        self.currentNetwork = self.networkCases[0]
        self.currentVesselId = None
                
        # description of the dataSet
        self.networkDescription = gtk.Label(' here is the description what can be rather long')
        self.networkDescription.set_size_request(400, 35)
        # # Combo boxes
        # vessel chooser
        self.comboBoxVessels = gtk.combo_box_new_text()
        self.comboBoxVessels.append_text("choose vessel")
        self.comboBoxVessels.connect("changed", self.on_changedCBvessels) 
        self.comboBoxVessels.set_active(0)  
        # Network SolutionData chooser
        self.comboBoxNetworks = gtk.combo_box_new_text()
        self.comboBoxNetworks.append_text("choose simulation case")
        self.comboBoxNetworks.connect("changed", self.on_changedNetworkComboBox) 
        self.comboBoxNetworks.set_size_request(250, 35)
        self.comboBoxNetworks.set_active(0)    
        textComboBoxNetworks = gtk.Label(' '.join(['Case', str(number)]))
        textComboBoxNetworks.set_size_request(100, 35)             
        
        # Case
        self.vBoxCase = gtk.VBox(False, 10)
        hBoxButtons = gtk.HBox(False, 10)
        hBoxDecription = gtk.HBox(False, 10)

        spacingText25Box = gtk.Label('')
        spacingText25Box.set_size_request(5, 35)
        spacingText25BoxDesc = gtk.Label('')
        spacingText25BoxDesc.set_size_request(25, 35)

        hBoxButtons.pack_start(spacingText25Box, fill=False, expand=False)
        hBoxButtons.pack_start(textComboBoxNetworks, fill=False, expand=False)
        hBoxButtons.pack_start(self.comboBoxNetworks, fill=False, expand=False)
        hBoxButtons.pack_start(self.comboBoxVessels, fill=False , expand=False)
         
        hBoxDecription.pack_start(spacingText25BoxDesc, fill=False, expand=False)
        hBoxDecription.pack_start(self.networkDescription, fill=False, expand=False)

        separator = gtk.HSeparator()
        
        self.vBoxCase.pack_start(hBoxButtons, expand=True, fill=False)
        self.vBoxCase.pack_start(hBoxDecription, expand=True, fill=False)
        self.vBoxCase.pack_start(separator, expand=True, fill=False)
        
    def updateVesselComboBox(self):
        '''
        Update the comboBox entries with vessel ids
        '''
        self.comboBoxVessels.get_model().clear()
        self.comboBoxVessels.append_text("choose vessel")
        for vesselName in self.networkInfo[self.currentNetwork][0]:
            self.comboBoxVessels.append_text(vesselName)
        self.comboBoxVessels.set_active(0)
                        
    def updateNetworkComboBox(self, networkCases, networkInfo):
        '''
        update the comboBox entries if available network names 
        have changed
        '''
        self.networkCases = networkCases
        self.networkInfo = networkInfo
        
        self.comboBoxNetworks.get_model().clear()
        for casName in networkCases:
            self.comboBoxNetworks.append_text(casName)     
        
        if len(networkCases) > 1:
            self.comboBoxNetworks.set_active(1)
            self.currentNetwork = self.networkCases[1]
        else:
            self.comboBoxNetworks.set_active(0)
            self.currentNetwork = self.networkCases[0]       
        
    def on_changedCBvessels(self, widget):
        '''
        call function of the combobox to choose vessel of the network
        '''
        cbIndex = widget.get_active()
        if cbIndex <= 0 : self.currentVesselId = None
        else: self.currentVesselId = self.networkInfo[self.currentNetwork][0][cbIndex - 1                                                                     ]
        
    def on_changedNetworkComboBox(self, widget):
        '''
        call function of the combobox to choose solution data set number
        '''
        cbIndex = widget.get_active()
        try: self.currentNetwork = self.networkCases[cbIndex]
        except: pass
        self.updateVesselComboBox()
        
        self.networkDescription.set_text(self.networkInfo[self.currentNetwork][1])

class Visualisation2DMainGUI(gtk.Window):
    '''
    Class defining the GUI of the visualisation main window
    '''
    def __init__(self):
        '''
        Initialize
        '''
        self.networkInfo = {"choose simulation case" : [[], '-'] }
        self.networkCases = ["choose simulation case"]
        
        super(Visualisation2DMainGUI, self).__init__()
        # variables for windows                
        self.set_size_request(550, 290)
        self.set_position(gtk.WIN_POS_CENTER)
        self.connect("destroy", gtk.main_quit)
        self.mainTitle = "2D Visualisation - "
        self.set_title(self.mainTitle)
        self.set_resizable(False)
                       
        # open Button
        self.buttonOpenSolutionData = gtk.Button("Open solutionData")
        self.buttonOpenSolutionData.connect("clicked", self.on_clickedLoad)
        self.buttonOpenSolutionData.set_size_request(200, 35)
        
        # open plot window
        self.buttonOpenPlotWindow = gtk.Button("Show plots")
        self.buttonOpenPlotWindow.connect("clicked", self.on_clickedPlots)
        self.buttonOpenPlotWindow.set_size_request(120, 35)
        
         # add new Case
        self.cases = []
        self.buttonAddCase = gtk.Button("Add case")
        self.buttonAddCase.connect("clicked", self.on_clickedAddCase)
        self.buttonAddCase.set_size_request(120, 35)
        
        # ExternalData
        self.extDataDescription = gtk.Label("-")
        self.extDataDescription.set_size_request(400, 35)
        # open Button
        self.buttonOpenExternalData = gtk.Button('Open external data for comparison')
        self.buttonOpenExternalData.connect("clicked", self.on_clickedLoadExternal)
        self.buttonOpenExternalData.set_size_request(250, 35)
        # enable check box
        self.buttonEnableExternalData = gtk.CheckButton("plot external data")
        self.buttonEnableExternalData.set_size_request(150, 35)
        self.buttonEnableExternalData.set_active(0)
        
        # alignment of the boxes
        self.vBox = gtk.VBox(False, 10)
        
        # Load And Plot buttons
        hBox1 = gtk.HBox(False, 10)
        
        spacingText25Box1 = gtk.Label('')
        spacingText25Box1.set_size_request(25, 35)
        
        hBox1.pack_start(spacingText25Box1, fill=False, expand=False)
        hBox1.pack_start(self.buttonOpenSolutionData, fill=False, expand=False)
        hBox1.pack_start(self.buttonOpenPlotWindow, fill=False, expand=False)
        separator1 = gtk.HSeparator()
        
        self.vBox.pack_start(hBox1, expand=False, fill=False)
        self.vBox.pack_start(separator1, expand=False, fill=False)
        
        # # add first case button
        newCase = Visualisation2DMainCase(len(self.cases) + 1)
        self.cases.append(newCase)
        self.vBox.pack_start(newCase.vBoxCase, expand=False, fill=False)
                
        # add more button
        hBoxAdd = gtk.HBox(False, 10)
        spacingText25Add = gtk.Label('')
        spacingText25Add.set_size_request(25, 35)
        hBoxAdd.pack_start(spacingText25Add, fill=False, expand=False)
        hBoxAdd.pack_start(self.buttonAddCase, fill=False, expand=False)
        separator2 = gtk.HSeparator()
        
        self.vBox.pack_start(hBoxAdd, expand=False, fill=False)
        self.vBox.pack_start(separator2, expand=False, fill=False)
                                
        # external Data Set
        
        spacingText25Box4text = gtk.Label('')
        spacingText25Box4text.set_size_request(15, 25)
        spacingText25Box4Desc = gtk.Label('')
        spacingText25Box4Desc.set_size_request(25, 35)
        
        hBoxExternalData = gtk.HBox(False, 10)
        hBoxExternalDecription = gtk.HBox(False, 10)
        hBoxExternalData.pack_start(spacingText25Box4text, fill=False, expand=False)
        hBoxExternalData.pack_start(self.buttonOpenExternalData, fill=False, expand=False)
        hBoxExternalData.pack_start(self.buttonEnableExternalData, fill=False, expand=False)        
        hBoxExternalDecription.pack_start(spacingText25Box4Desc, fill=False, expand=False)
        hBoxExternalDecription.pack_start(self.extDataDescription, fill=False, expand=False)
        
        # external Data Set
        self.vBox.pack_start(hBoxExternalData, expand=False, fill=False)
        self.vBox.pack_start(hBoxExternalDecription, expand=False, fill=False)       
                     
        self.add(self.vBox)
        self.show_all()
        
    def on_clickedAddCase(self, widget):
        '''
        add new simulation case to compare plots
        bounded to 5 cases now ..
        '''
        if len(self.cases) < 5:
            newCase = Visualisation2DMainCase(len(self.cases) + 1)
            newCase.updateNetworkComboBox(self.networkCases, self.networkInfo)
            self.cases.append(newCase)
            self.vBox.pack_start(newCase.vBoxCase, expand=False, fill=False)
            self.vBox.reorder_child(newCase.vBoxCase, len(self.cases) + 1)
            
            width, height = self.get_size()
            self.set_size_request(width, height + 102)
            
            self.show_all()
                
    def on_clickedLoad(self, widget):
        pass
    
    def on_clickedPlots(self, widget):
        pass
        
    def on_clickedLoadExternal(self):
        pass
        
    def on_clickedLoad(self):
        pass
    
    def on_clickedLoadExternal(self):
        pass
    
class Visualisation2DMain(Visualisation2DMainGUI):
    '''
    Class for the Main GUI window
    '''
    def __init__(self, networkName=None, dataNumber=None, connect=None):
        super(Visualisation2DMain, self).__init__()

        self.externalData = None
        
        self.networks     = { "choose simulation case"   : None}
        self.networkCases = [ "choose simulation case"]
        self.networkInfo  = { "choose simulation case"   : [[], '-'] }
        
        self.loadVascularNetwork(networkName, dataNumber)

    def on_clickedPlots(self, widget):
        '''
        Open new plot window to plot information of all selected vessels of the selected cases  
        '''
        # # check out selected networks and ids
        selectedNetworks = []
        selectedVesselIds = []
        for case in self.cases:
            currentNetwork = case.currentNetwork
            currentVesselId = case.currentVesselId
            if currentNetwork != "choose simulation case":
                if currentVesselId:
                    currentVesselId = int(currentVesselId.split('-')[0])
                    if currentVesselId not in selectedVesselIds:
                        selectedNetworks.append(self.networks[currentNetwork])
                        selectedVesselIds.append(currentVesselId)
        # # check selected external data
        selectedExternalData = None
        if self.buttonEnableExternalData.get_active() == True:
            print " oh yes and i will plot external data if given"
            selectedExternalData = self.externalData
        # # open plot window        
        if selectedNetworks != []:
            Visualisation2DPlotWindow(selectedNetworks, selectedVesselIds, selectedExternalData)

    def loadVascularNetwork(self, networkName, dataNumber):
        '''
        This function actually loads the vascular networks with the given names and ids
            networkName
        The network is added to the existing ones if is not existing
        '''
        # load vascular network
        vascularNetwork = loadSolutionDataFile(networkName, dataNumber) 
        # # save it and refresh GUi setup
        networkSolutionName = '_'.join([networkName, dataNumber])  
        # # add data name and corresponding network
        self.networks[networkSolutionName] = vascularNetwork
        # # names of the solution data sets available sorted keys of the stuff above
        self.networkCases.append(networkSolutionName)
        # # get vessel names
        vesselNames = []
        for vessel in vascularNetwork.vessels.itervalues():
            vesselNames.append('{:3}-{}'.format(vessel.Id, vessel.name))
        # # reference network showing network case name and [corresponding number of vessels, vesselNames, description]
        self.networkInfo[networkSolutionName] = [ vesselNames, vascularNetwork.description]
        
        # # update comboBoxes
        for case in self.cases:
            case.updateNetworkComboBox(self.networkCases, self.networkInfo)

    def on_clickedLoad(self, widget):
        '''
        Call function of the Open Solution Data Button
        '''
        fileFilter = gtk.FileFilter()
        fileFilter.set_name("SolutionData")
        fileFilter.add_pattern("*.pickle")
        
        filenames = self.LoadDialog(fileFilter)
        
        for filename in filenames:
            networkName = filename.split('/')[-1].split('_SolutionData_')[0]
            dataNumber = filename.split('/')[-1].split('_SolutionData_')[1].split('.')[0]
            if '_'.join([networkName, dataNumber]) not in self.networkCases:
                self.loadVascularNetwork(networkName, dataNumber)
        
                            
    def LoadDialog(self, fileFilter):
        '''
        Dialog window function: 
        displaying the file browser where solution data files can be selected
        to be loaded into the program.
        
        return: [filename,..] (abs-path to the files)
        '''
        dialog = gtk.FileChooserDialog("Open Solution Data File..",
                                     None,
                                     gtk.FILE_CHOOSER_ACTION_OPEN,
                                     (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                      gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        directory = ''.join([cur, '/../NetworkFiles/'])
        dialog.set_current_folder(directory)
        
        dialog.set_select_multiple(True)
        
        dialog.add_filter(fileFilter)
        
        filenames = []
        
        response = dialog.run()
        
        if response == gtk.RESPONSE_OK:
            filenames = dialog.get_filenames()
        elif response == gtk.RESPONSE_CANCEL:
            print 'Closed, no files selected'
        dialog.destroy()    
        return filenames
    

    def on_clickedLoadExternal(self, widget):
        fileFilter = gtk.FileFilter()
        fileFilter.add_pattern("*.v1dfExD")
        filenames = self.LoadDialog(fileFilter)
        
        try:
            fileName = filenames[0]
            self.externalData = loadExternalDataSet(fileName)
            self.extDataDescription.set_text(self.externalData['Description'])
        except: 
            pass


if __name__ == '__main__':
               
    optionsDict = parseOptions(['f', 'n', 'c'], visualisationOnly=True)
    
    networkName = optionsDict['networkName']
    dataNumber = optionsDict['dataNumber']
    connect = optionsDict['connect']
             
    Visualisation2DMain(networkName=networkName, dataNumber=dataNumber, connect=connect)
    gtk.main()
