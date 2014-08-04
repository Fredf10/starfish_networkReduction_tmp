import pydot
import xdot
import gtk
import cPickle

class MyDotWindow(xdot.DotWindow):

    def __init__(self):
        xdot.DotWindow.__init__(self)


# class for the graph

class Graph(object):
    
    def __init__(self):
        '''
        Constructor
        '''
        # create the graph
        self.graph = pydot.Dot(graph_type ='graph')
    
    
    def addNodeToGraph(self, node):

        # create end node if not existing
        if self.graph.get_node(node) == None or self.graph.get_node(node) == []:
            self.graph.add_node(pydot.Node(node))
        try:
            self.graph.get_node(node)[0].set_shape('circle')
            self.graph.get_node(node)[0].set_fontsize('9')
        except:
            self.graph.get_node(node).set_shape('circle')
            self.graph.get_node(node).set_fontsize('9')


    def addVesselToGraph(self, name, startNode, endNode, typeOfDaughter, visible):

        node1 = str(startNode)
        node2 = str(endNode)
    
        # create start node if not existing
        self.addNodeToGraph(node1)
        self.addNodeToGraph(node2)
        
        # create vessel
        if typeOfDaughter == 'leftDaughter':      currentColor = 'red3'
        elif typeOfDaughter == 'rightDaughter' : currentColor = 'steelblue'
        elif typeOfDaughter == 'bif': currentColor = 'green'
        elif typeOfDaughter == 'anas': currentColor = 'chocolate1'
        else: currentColor = 'black'
        self.graph.add_edge(pydot.Edge(node1,node2, label= name, color=currentColor, style = visible))
        
        
    def add_bc(self,networkNode, bcNode, name, rootCON = False):
        ## creates boundary node
        
        # check if networkNode exists
        if self.graph.get_node(networkNode) == None or self.graph.get_node(networkNode) == []:
            print 'Error: something went wrong, not such network nodes'
        
        # creates boundary condition node
        if self.graph.get_node(bcNode) == None or self.graph.get_node(bcNode) == []:
            self.graph.add_node(pydot.Node(bcNode, label= name))
            if type(self.graph.get_node(bcNode)) == list:
                self.graph.get_node(bcNode)[0].set_shape('box')
                self.graph.get_node(bcNode)[0].set_color('darkgreen')
            else:
                self.graph.get_node(bcNode).set_shape('box')
                self.graph.get_node(bcNode).set_color('darkgreen')
                
        # creates edge connection
        if rootCON: 
            self.graph.add_edge(pydot.Edge(bcNode,networkNode,color = 'darkgreen'))
        else: 
            self.graph.add_edge(pydot.Edge(networkNode,bcNode,color = 'darkgreen'))


    def getGraph(self):
        return self.graph.to_string()
    
    
    def update_graph(self, vascularNetwork, window):
     
        self.resetGraph()
        
        if vascularNetwork != None: 
        
            vascularNetwork.evaluateConnections()
            vascularNetwork.findStartAndEndNodes()
            
            if vascularNetwork.vessels.keys() != []:
                ## go through network as a binary tree and create vessel graph
                
                root = vascularNetwork.root
                self.addVesselToGraph(' '.join([str(root),vascularNetwork.vessels[root].name]), str(vascularNetwork.vessels[root].startNode), str(vascularNetwork.vessels[root].endNode), '', ' ')
                        
                for LM,RM,LD,RD in vascularNetwork.treeTraverseConnections:
                    ## link
                    if RM == None and RD == None:
                        #add leftDaughter
                        self.addVesselToGraph(' '.join([str(LD),vascularNetwork.vessels[LD].name]), str(vascularNetwork.vessels[LD].startNode), str(vascularNetwork.vessels[LD].endNode), 'leftDaughter', ' ')
                        
                    
                    ## bifurcation
                    elif RM == None:
                        # add leftDaughter
                        self.addVesselToGraph(' '.join([str(LD),vascularNetwork.vessels[LD].name]), str(vascularNetwork.vessels[LD].startNode), str(vascularNetwork.vessels[LD].endNode), 'leftDaughter', ' ')
                        # add rightDaughter
                        self.addVesselToGraph( ' '.join([str(RD),vascularNetwork.vessels[RD].name]), str(vascularNetwork.vessels[RD].startNode), str(vascularNetwork.vessels[RD].endNode), 'rightDaughter', ' ')
                        # connect leftDaughter to right daughter
                        #addVesselToGraph(' ', str(vascularNetwork.vessels[RD].endNode), str(vascularNetwork.vessels[LD].endNode), 'bif', 'invis')
                        #addVesselToGraph(' ', str(vascularNetwork.vessels[LD].endNode), str(vascularNetwork.vessels[RD].endNode), 'bif', 'invis')
                                  
                    ## anastomosis
                    elif RD == None:
                        print LM,RM,LD,RD 
                        # connect LM and RM
                        #addVesselToGraph(' ', str(vascularNetwork.vessels[RM].startNode), str(vascularNetwork.vessels[LM].startNode), 'anas', ' ')
                        # add rightDaughter
                        self.addVesselToGraph(' '.join([str(LD),vascularNetwork.vessels[LD].name]), str(vascularNetwork.vessels[LD].startNode), str(vascularNetwork.vessels[LD].endNode), 'leftDaughter', ' ')
                   
                ## create boundary condition nodes
            
                ## root node bc condition
                # find applied bc conditions
                bcNames = ['?']
                bcTrue = False
                if root in vascularNetwork.boundaryConditions.keys() and vascularNetwork.boundaryConditions[root] != {}:
                        bcNames.pop(0)
                        for boundaryInstance in vascularNetwork.boundaryConditions[root]:
                            bcNames.append(boundaryInstance.name)
                for bcName in bcNames:
                    if '_' not in bcName:
                        self.add_bc(str(vascularNetwork.vessels[root].startNode), ''.join(['bcR',str(root),str(bcNames.index(bcName))]), bcName, rootCON = True) 
                        bcTrue = True
                if bcTrue == False:
                    self.add_bc(str(vascularNetwork.vessels[root].startNode), ''.join(['bcR',str(root),'0']), '?', rootCON = True)
                # end node bc conditions
                for endNode in vascularNetwork.boundaryVessels:
                    bcNames = ['?']
                    if endNode in vascularNetwork.boundaryConditions.keys() and vascularNetwork.boundaryConditions[endNode] != {}:
                            bcNames.pop(0)
                            for boundaryInstance in vascularNetwork.boundaryConditions[endNode]:
                                bcNames.append(boundaryInstance.name)
                    if endNode == root:
                        bcTrue = False
                        for bcName in bcNames:
                            if '_' in bcName: 
                                self.add_bc(str(vascularNetwork.vessels[endNode].endNode), ''.join(['bcEq',str(endNode),str(bcNames.index(bcName))]), bcName) 
                                bcTrue = True
                        if bcTrue == False:
                            self.add_bc(str(vascularNetwork.vessels[endNode].endNode), ''.join(['bcEq',str(endNode),'0']), '?')
                    else:
                        for bcName in bcNames:
                            self.add_bc(str(vascularNetwork.vessels[endNode].endNode), ''.join(['bcEq',str(endNode),str(bcNames.index(bcName))]), bcName) 
            # save the vascularNetwork temporary                          
            tempSave(vascularNetwork)
        window.set_dotcode(self.graph.to_string())
        window.show()
    
    def resetGraph(self):
        self.graph = pydot.Dot(graph_type ='graph')
    

def bifTest(vessels,testNode):
    # bifTest returns true if less then 2 nodes start at testNode
    # bifTest retruns false if already 2 nodes start here
    bifCounter = 0
    for vessel_i in vessels.itervalues():
        if vessel_i.startNode == testNode:
            bifCounter = bifCounter + 1
    if bifCounter <2:
        return True
    else:
        return False
    
def enterFilename(filename, endString,recentFilenames = None):
    print "     current filename: ",filename
    filenameT = str(raw_input("     enter/change filename (only! ENTER to use current filename):\n     "))
    
    if filenameT == "":
        if filename != None:
            filenameT = filename+endString
        else:
            filenameT = None
            
    if filenameT != None:       
        if filenameT in [str(int(i)) for i in range(0,9)] and len(filenameT)==1 and recentFilenames!= None:
            number = int(filenameT)
            if len(recentFilenames)>= number:
                filenameT = (recentFilenames[number-1]+endString)
                
        if endString not in filenameT:
            filenameT = (filenameT+endString)
        
    return filenameT
    
def tempSave(vascularNetwork):
    FILE = open('.tempData.pickle',"w")
    vascularNetwork.quiet = True
    # store pickle
    vascularNetwork.prepareToSave()
    cPickle.dump(vascularNetwork, FILE, protocol=2)
    FILE.close()
    
def tempLoad():
    try:
        FILE = open('.tempData.pickle',"rb")
        vascularNetwork = cPickle.load(FILE)
        FILE.close()    
    except:
        print "WARNING vnc-classes.py - could nof find tempData.pickle"
    vascularNetwork.evaluateConnections()
    vascularNetwork.quiet = True
    
    return vascularNetwork


def findAllDaughters(vascularNetwork, motherVesselID):
    daughters = []
    viz = []
    root = motherVesselID
    if vascularNetwork.vessels[root].leftDaughter != None:
        viz.append(root)
    # loop through tree until all daughters are added to the graph
    while len(viz) != 0:
        # get the mother vessel (already added) and add its daughters
        motherVessel = viz.pop(0)
        # find left daughter
        leftDaughter = vascularNetwork.vessels[motherVessel].leftDaughter
        # add left daughter
        daughters.append(leftDaughter)
        # check if leftDaughter has also daughters 
        if vascularNetwork.vessels[leftDaughter].leftDaughter != None:
            viz.append(leftDaughter)
        # find right daughter
        if vascularNetwork.vessels[motherVessel].rightDaughter != None:
            rightDaughter = vascularNetwork.vessels[motherVessel].rightDaughter
            # add right daughter
            daughters.append(rightDaughter)
            # check if rightDaughter has also daughters
            if vascularNetwork.vessels[rightDaughter].leftDaughter != None:
                viz.append(rightDaughter)
    return daughters
    
