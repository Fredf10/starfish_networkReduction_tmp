#!/usr/bin/env python
import gtk
import gobject

import matplotlib.pyplot as plt   

from matplotlib.figure import Figure
# uncomment to select /GTK/GTKAgg/GTKCairo
#from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append(cur+'/../UtilityLib')
from moduleXML import loadNetworkFromXML 
from processing import linearWaveSplitting
from processing import nonLinearWaveSplitting
from processing import minMaxFunction

from moduleStartUp import parseOptions

from modulePickle import loadSolutionDataFile
from modulePickle import loadExternalDataSet


sys.path.append(cur+'/../NetworkLib')
from classVascularNetwork import VascularNetwork 

from matplotlib.font_manager import FontProperties

import itertools
import subprocess

from optparse import OptionParser
import cPickle
import numpy as np

from pprint import pprint as pp

from numpy import arange, sin, pi

class Visualisation2DMainCase(object):
    def __init__(self, number):
        
        self.networkInfo = {"choose simulation case" : [[],'-'] }
        self.networkCases = ["choose simulation case"]
        self.currentNetwork  = self.networkCases[0]
        self.currentVesselId = None
                
        # description of the dataSet
        self.networkDescription = gtk.Label(' here is the description what can be rather long')
        self.networkDescription.set_size_request(400,35)
        ## Combo boxes
        # vessel chooser
        self.comboBoxVessels = gtk.combo_box_new_text()
        self.comboBoxVessels.append_text("choose vessel")
        self.comboBoxVessels.connect("changed", self.on_changedCBvessels) 
        self.comboBoxVessels.set_active(0)  
        # Network SolutionData chooser
        self.comboBoxNetworks = gtk.combo_box_new_text()
        self.comboBoxNetworks.append_text("choose simulation case")
        self.comboBoxNetworks.connect("changed", self.on_changedNetworkComboBox) 
        self.comboBoxNetworks.set_size_request(250,35)
        self.comboBoxNetworks.set_active(0)    
        textComboBoxNetworks = gtk.Label(' '.join(['Case',str(number)]))
        textComboBoxNetworks.set_size_request(100,35)             
        
        # Case
        self.vBoxCase = gtk.VBox(False,10)
        hBoxButtons = gtk.HBox(False,10)
        hBoxDecription = gtk.HBox(False,10)

        spacingText25Box = gtk.Label('')
        spacingText25Box.set_size_request(5,35)
        spacingText25BoxDesc = gtk.Label('')
        spacingText25BoxDesc.set_size_request(25,35)

        hBoxButtons.pack_start(spacingText25Box,fill = False, expand = False)
        hBoxButtons.pack_start(textComboBoxNetworks,fill = False, expand = False)
        hBoxButtons.pack_start(self.comboBoxNetworks,fill = False, expand = False)
        hBoxButtons.pack_start(self.comboBoxVessels, fill = False , expand = False)
         
        hBoxDecription.pack_start(spacingText25BoxDesc,fill = False, expand = False)
        hBoxDecription.pack_start(self.networkDescription,fill = False, expand = False)

        separator = gtk.HSeparator()
        
        self.vBoxCase.pack_start(hBoxButtons, expand = True, fill = False)
        self.vBoxCase.pack_start(hBoxDecription, expand = True, fill = False)
        self.vBoxCase.pack_start(separator, expand = True, fill = False)
        
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
            self.currentNetwork  = self.networkCases[1]
        else:
            self.comboBoxNetworks.set_active(0)
            self.currentNetwork  = self.networkCases[0]       
        
    def on_changedCBvessels(self,widget):
        '''
        call function of the combobox to choose vessel of the network
        '''
        cbIndex = widget.get_active()
        if cbIndex <= 0 : self.currentVesselId = None
        else: self.currentVesselId = self.networkInfo[self.currentNetwork][0][cbIndex-1
                                                                              ]
        
    def on_changedNetworkComboBox(self,widget):
        '''
        call function of the combobox to choose solution data set number
        '''
        cbIndex = widget.get_active()
        try: self.currentNetwork  = self.networkCases[cbIndex]
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
        self.networkInfo = {"choose simulation case" : [[],'-'] }
        self.networkCases = ["choose simulation case"]
        
        super(Visualisation2DMainGUI, self).__init__()
        # variables for windows                
        self.set_size_request(550,290)
        self.set_position(gtk.WIN_POS_CENTER)
        self.connect("destroy", gtk.main_quit)
        self.mainTitle = "2D Visualisation - "
        self.set_title(self.mainTitle)
        self.set_resizable(False)
                       
        # open Button
        self.buttonOpenSolutionData = gtk.Button("Open solutionData")
        self.buttonOpenSolutionData.connect("clicked", self.on_clickedLoad)
        self.buttonOpenSolutionData.set_size_request(200,35)
        
        # open plot window
        self.buttonOpenPlotWindow = gtk.Button("Show plots")
        self.buttonOpenPlotWindow.connect("clicked", self.on_clickedPlots)
        self.buttonOpenPlotWindow.set_size_request(120,35)
        
         # add new Case
        self.cases = []
        self.buttonAddCase = gtk.Button("Add case")
        self.buttonAddCase.connect("clicked", self.on_clickedAddCase)
        self.buttonAddCase.set_size_request(120,35)
        
        # ExternalData
        self.extDataDescription = gtk.Label("-")
        self.extDataDescription.set_size_request(400,35)
        # open Button
        self.buttonOpenExternalData = gtk.Button('Open external data for comparison')
        self.buttonOpenExternalData.connect("clicked", self.on_clickedLoadExternal)
        self.buttonOpenExternalData.set_size_request(250,35)
        # enable check box
        self.buttonEnableExternalData = gtk.CheckButton("plot external data")
        self.buttonEnableExternalData.set_size_request(150,35)
        self.buttonEnableExternalData.set_active(0)
        
        #alignment of the boxes
        self.vBox = gtk.VBox(False,10)
        
        # Load And Plot buttons
        hBox1 = gtk.HBox(False,10)
        
        spacingText25Box1 = gtk.Label('')
        spacingText25Box1.set_size_request(25,35)
        
        hBox1.pack_start(spacingText25Box1,fill = False, expand = False)
        hBox1.pack_start(self.buttonOpenSolutionData,fill = False, expand = False)
        hBox1.pack_start(self.buttonOpenPlotWindow,fill = False, expand = False)
        separator1 = gtk.HSeparator()
        
        self.vBox.pack_start(hBox1, expand = False, fill = False)
        self.vBox.pack_start(separator1, expand = False, fill = False)
        
        ## add first case button
        newCase = Visualisation2DMainCase(len(self.cases)+1)
        self.cases.append(newCase)
        self.vBox.pack_start(newCase.vBoxCase, expand = False, fill = False)
                
        # add more button
        hBoxAdd = gtk.HBox(False,10)
        spacingText25Add = gtk.Label('')
        spacingText25Add.set_size_request(25,35)
        hBoxAdd.pack_start(spacingText25Add,fill = False, expand = False)
        hBoxAdd.pack_start(self.buttonAddCase,fill = False, expand = False)
        separator2 = gtk.HSeparator()
        
        self.vBox.pack_start(hBoxAdd, expand = False, fill = False)
        self.vBox.pack_start(separator2, expand = False, fill = False)
                                
        # external Data Set
        
        spacingText25Box4text = gtk.Label('')
        spacingText25Box4text.set_size_request(15,25)
        spacingText25Box4Desc = gtk.Label('')
        spacingText25Box4Desc.set_size_request(25,35)
        
        hBoxExternalData = gtk.HBox(False,10)
        hBoxExternalDecription = gtk.HBox(False,10)
        hBoxExternalData.pack_start(spacingText25Box4text, fill = False, expand = False)
        hBoxExternalData.pack_start(self.buttonOpenExternalData, fill = False, expand = False)
        hBoxExternalData.pack_start(self.buttonEnableExternalData, fill = False, expand = False)        
        hBoxExternalDecription.pack_start(spacingText25Box4Desc,fill = False, expand = False)
        hBoxExternalDecription.pack_start(self.extDataDescription,fill = False, expand = False)
        
        # external Data Set
        self.vBox.pack_start(hBoxExternalData, expand = False, fill = False)
        self.vBox.pack_start(hBoxExternalDecription, expand = False, fill = False)       
                     
        self.add(self.vBox)
        self.show_all()
        
    def on_clickedAddCase(self, widget):
        '''
        add new simulation case to compare plots
        '''
        newCase = Visualisation2DMainCase(len(self.cases)+1)
        newCase.updateNetworkComboBox(self.networkCases, self.networkInfo)
        self.cases.append(newCase)
        self.vBox.pack_start(newCase.vBoxCase, expand = False, fill = False)
        self.vBox.reorder_child(newCase.vBoxCase,len(self.cases)+1)
        
        width,height = self.get_size()
        self.set_size_request(width,height+102)
        
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
    def __init__(self, networkName = None, dataNumber = None, connect = None):
        super(Visualisation2D, self).__init__()

        self.externalData = None
        
        self.networks     = { "choose simulation case"   : None}
        self.networkCases = [ "choose simulation case"]
        self.networkInfo  = { "choose simulation case"   : [[],'-'] }
        
        self.loadVascularNetwork(networkName, dataNumber)

    def on_clickedPlots(self, widget):
        '''
        Open new plot window to plot information of all selected vessels of the selected cases  
        '''
        selectedNetworks = []
        selectedVesselIds = []
        for case in self.cases:
            currentNetwork  = case.currentNetwork
            currentVesselId = case.currentVesselId
            if currentNetwork != "choose simulation case":
                if currentVesselId and currentVesselId not in selectedVesselIds:
                    selectedCases.append(self.networks[currentNetwork])
                    selectedVesselIds.append(currentVesselId)

        print "Follwing networks selected:", selectedNetworks
        print "with following network ids", selectedVesselIds

        selectedExternalData = None
        if self.buttonEnableExternalData.get_active() == True:
            print " oh yes and i will plot external data if given"
            selectedExternalData = self.externalData
        
        #Visualisation2DPlotWindow(selectedNetworks, selectedVesselIds, selectedExternalData)

    def loadVascularNetwork(self, networkName, dataNumber):
        '''
        This function actually loads the vascular networks with the given names and ids
            networkName
        The network is added to the existing ones if is not existing
        '''
        # load vascular network
        vascularNetwork = loadSolutionDataFile(networkName, dataNumber) 
        ## save it and refresh GUi setup
        networkSolutionName = '_'.join([networkName,dataNumber])  
        ## add data name and corresponding network
        self.networks[networkSolutionName] = vascularNetwork
        ## names of the solution data sets available sorted keys of the stuff above
        self.networkCases.append(networkSolutionName)
        ## get vessel names
        vesselNames = []
        for vessel in vascularNetwork.vessels.itervalues():
            vesselNames.append('{:3} {}'.format(vessel.Id, vessel.name))
        ## reference network showing network case name and [corresponding number of vessels, vesselNames, description]
        self.networkInfo[networkSolutionName] = [ vesselNames, vascularNetwork.description]
        
        ## update comboBoxes
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
            dataNumber  = filename.split('/')[-1].split('_SolutionData_')[1].split('.')[0]
            if '_'.join([networkName,dataNumber]) not in self.networkCases:
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
        
        directory = ''.join([cur,'/../NetworkFiles/'])
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
               
    optionsDict = parseOptions(['f','n','c'], visualisationOnly = True)
    
    networkName  = optionsDict['networkName']
    dataNumber   = optionsDict['dataNumber']
    connect      = optionsDict['connect']
             
    Visualisation2DMain(networkName = networkName, dataNumber = dataNumber, connect = connect)
    gtk.main()