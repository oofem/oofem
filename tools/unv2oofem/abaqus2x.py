'''
Created on Aug 5, 2013

@author: Carl Sandstrom

This Python module is used to parse an Abaqus .inp file and store the data in the same way
as the unv2x module. In order to do so, imaginary elements on the faces of the 3D element are 
created from surfaces. The resulting element type depends on the input type and are defined 
in the method setupElement (for example, the c3d10 element are convertyed to type 118 in UNV 
and has surface element of type 42) 
'''
from FEM import *
from oofemctrlreader import *

class ElementProperties:
    def __init__(self, name, unvType=-1, edgeMask=None, faceMask=None, nnodes=-1, surfaceElementType=-1):
        if (edgeMask != None):
            self.name = name.lower() # string - name as Abaqus elements
            self.unvType = unvType
            self.edgeMask = edgeMask # 2D array expressing node masks for OOFEM's edge 1,2,..., original Abaqus node numbering
            self.faceMask = faceMask # 2D array expressing node masks for OOFEM's face 1,2,..., original Abaqus node numbering
            self.nnodes = nnodes
            self.surfaceElementType = surfaceElementType # Element type on surface
            
class AbaqusParser:
    '''
    classdocs
    '''
    
    def __init__(self,filename):
        self.file = None
        self.filename = filename
        self.FEM = FEM()
        self.ElementConfiguration = self.setupElements()
        
    def mapping(self):
        oofem_elemProp = []
        oofem_elemProp.append(oofem_elementProperties("None", [0], [], []))#leave this line [0] as it is
        oofem_elemProp.append(oofem_elementProperties("RepresentsBoundaryLoad", [],[],[]))#special element representing boundary load
        oofem_elemProp.append(oofem_elementProperties("tet21stokes", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [], [[0,1,2,4,5,6],[0,1,3,4,7,8],[1,2,3,5,8,9],[0,2,3,6,7,9]]))
#        oofem_elemProp.append(oofem_elementProperties("PlaneStress2DXFEM", [0,1,2,3], [], [[0,1],[1,2],[2,3],[3,0]]))
        oofem_elemProp.append(oofem_elementProperties("TrPlaneStress2DXFEM", [0,2,1], [], [[0,2],[2,1],[1,0]]))
        oofem_elemProp.append(oofem_elementProperties("intelline1", [0,1,2,3], [], []))
        return oofem_elemProp
    
    def setupElements(self):
        es = [];
        es.append(ElementProperties("c3d10", 118, [], [[0,1,2,4,5,6],[0,1,3,4,7,8],[1,2,3,5,8,9],[0,2,3,6,7,9]], 10, 42))
        # es.append(ElementProperties("c3d20R", 118, [], [[0,1,2,4,5,6],[0,1,3,4,7,8],[1,2,3,5,8,9],[0,2,3,6,7,9]], 20, 42))

        # 3-node triangle
        es.append(ElementProperties("cpe3", 1, [], [[0,1],[1,2],[2,3],[3,0]], 3, 42))

        # 2d interface element  
        es.append(ElementProperties("coh1d4", 2, [], [], 4))
        return es
        
    def _read_multiple_lines(self):
        ''' Read a line from input file distributed across several rows. If a row end with a comma, it resumes on the next line '''
        thisLine = ""
        doContinue = True

        while True:
            start_of_line = self.file.tell()
            newLine = self.file.readline().strip()

            if len(newLine) > 0:
                if newLine[0]=='*':
                    self.file.seek(start_of_line)
                    if thisLine=="":
                        thisLine = newLine
                    break
                thisLine = thisLine + newLine
                if thisLine[-1] != ',':
                    thisLine = thisLine + ','
            
        return thisLine

    def _read_nodes(self):
        self.file.readline()
        while True:
            start_of_line = self.file.tell()
            nodeData = self.file.readline().strip().split(",")
            
            #print 'nodeData: ',nodeData
            if len(nodeData) > 0:
                if len(nodeData[0]) > 0:
                    if nodeData[0][0] == '*':
                        self.file.seek(start_of_line)
                        break
            
                    id = int(nodeData[0])
                    coords = map(float, nodeData[1:4])
            
                    self.FEM.nodes.append(Node(id, coords))
                    self.FEM.nnodes = self.FEM.nnodes+1            

    def _read_elements(self):
         elementData = self.file.readline().strip().split(",")
         elementType = elementData[1].split("=")[1]
         
         nnodes = -1;
         
         for es in self.ElementConfiguration:
             if elementType.lower() == es.name:
                 nnodes = es.nnodes
                 elType = es.unvType
                 
         if nnodes == -1:
             print ("Element not supported")
             exit(0)
         
         while True:
             while True:
                 start_of_line = self.file.tell()
                 data = self._read_multiple_lines()
                 
                 if data[0]=='*': # Reached en of element section
                     break
                 
                 elementSection=[]
                 elementSection.extend(data.split(","))
                 
                 while '' in elementSection:
                    elementSection.remove('')                 
                 
                 # Add element to FEM object (id,type,material,color,nnodes,cntvt)
                 for i in range (0, len(elementSection), nnodes+1):
                     thisElement = elementSection[i:(i+nnodes+1)]
                     self.FEM.elems.append(Element(int ( thisElement[0] ), elType, 0, 0, nnodes, map(int, thisElement[1:])))

             if data[0]=='*':
                 self.file.seek(start_of_line)
                 break
    
    def _read_generate(self):
        outList = []
        start_of_line = self.file.tell()
        generateData = self._read_multiple_lines().split(",")
        doContinue = generateData[0][0] != '*'
        
        while doContinue:
            start = int ( generateData[0] )
            end = int ( generateData[1] )
            if len(generateData) > 2:
                step = int ( generateData[2] )
            else:
                step = 1
            outList.extend( range(start, end+1, step) )
            start_of_line = self.file.tell()
            generateData = self._read_multiple_lines().split(",")
            doContinue = generateData[0][0] != '*'

        self.file.seek(start_of_line)

        return outList

    def _read_set(self):
        setData = self.file.readline().strip().split(",")
        data = []
        setName = setData[1].split("=")[1]
        
        doGenerate = False;
        
        for s in setData:
            if s.lower().strip()=='generate':
                doGenerate = True
                break 
        
        if doGenerate:
            gData = self._read_generate()
            data.extend( gData )
        else:
            dataLine = self._read_multiple_lines().split(',')
            while '' in dataLine:
                dataLine.remove('')
            data = map(int, dataLine)
            
        return setName, data
            
    def _read_element_set(self):
        elementSetName, elements = self._read_set()
        
        ID = len(self.FEM.elemsets)+1;
        eset = Group(ID, elementSetName)
        eset.type = 8
        eset.items.extend(elements)
        eset.nitems = len(eset.items)
        self.FEM.elemsets.append(eset)
        self.FEM.nelemsets = len(self.FEM.elemsets)
        pass
    
    def _read_node_set(self):
        nodeSetName, nodes = self._read_set()
        
        ID = len(self.FEM.nodesets)+1
        nset = Group(ID, nodeSetName)
        nset.type = 7
        nset.items.extend ( nodes )
        nset.nitems = len(nset.items)
        self.FEM.nodesets.append(nset)
        self.FEM.nnodesets = len(self.FEM.nodesets)
        pass

    def _read_surface(self):
        surfaceData = self.file.readline().strip().split(",")

        header = []
        elements = []
        
        for s in surfaceData:
            header.append(s.strip().split('='))
            if header[-1][0].lower() == 'name':
                elementSetName =  header[-1][1]
            
        while True:
            start_of_line = self.file.tell()
            lineData = self.file.readline().strip().split(',')
            
            if lineData[0][0]=='*':
                self.file.seek(start_of_line)
                break;
            
            face = int( lineData[1].strip()[1] )
            setName = lineData[0].strip()
            
            # Create virtual element from the set 
            for es in self.FEM.elemsets:
                if es.name == setName:
                    # es is the element set containing elements having face "face" facing the surface
                    for e in es.items:
                        # e is the index of the current element
                        eNodes = self.FEM.elems[e-1].cntvt
                        eProp = []
                        for t in self.ElementConfiguration:
                            if t.unvType == self.FEM.elems[e-1].type:
                                eProp = t
                                break;
                        
                        surfaceNodes = []
                        
                        for id in eProp.faceMask[face-1]:
                            surfaceNodes.append(eNodes[id])
                            
                        # Create element
                        eID = len(self.FEM.elems)+1
                        elements.append(eID)
                        self.FEM.elems.append(Element(eID, eProp.surfaceElementType, 0, 0, len( surfaceNodes ), surfaceNodes))
                    
                    #print es 
                    break
                
        # Create element set
        ID = len(self.FEM.elemsets)+1;
        eset = Group(ID, elementSetName)
        eset.type = 8
        eset.items.extend(elements)
        eset.nitems = len(eset.items)
        self.FEM.elemsets.append(eset)
        self.FEM.nelemsets = len(self.FEM.elemsets)
    def parse(self):
        """ parse Abaqus file to fill the FEM data structure"""
        self.file=open(self.filename,'r')
        
        while True:
            start_of_line = self.file.tell()
            line = self.file.readline()
            keyword = line.strip().split(",")[0].lower()
            if keyword == "*heading":
                self.heading = self.file.readline()
            elif keyword == "*node":
                self.file.seek(start_of_line)
                self._read_nodes()
            elif keyword == "*element":
                self.file.seek(start_of_line)
                self._read_elements()
            elif keyword == "*elset":
                self.file.seek(start_of_line)
                self._read_element_set()
            elif keyword == "*nset":
                self.file.seek(start_of_line)
                self._read_node_set()
            elif keyword == "*surface":
                self.file.seek(start_of_line)
                self._read_surface()                
            elif not line:
                break
                        
        self.file.close()
        
        # Determine which sets that represents surfaces/edges

        return self.FEM
            
            
