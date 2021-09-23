# -*- coding: utf-8 -*-
# UNV2X base code, J.Cugnoni, www.caelinux.com,9.2006
# Licence: GPL
# Purpose: this is a set of objects & functions to read a Universal File
#          into a simple FEM object structure in order to simplify the
#          conversion of Mesh definitions from UNV to any other format
# How it works:
#          This code is based on two main objects:
#           1) a simple FEM object structure to store nodes, elements and groups
#           2) a UNVParser which provides a simple & modular solution to read
#              some datasets from UNV file and store them in a FEM object structure
# To Do:   add UNV data handler functions for other datasets (to be defined)
#          add your own code to write the model into your own file format
# For documentation see http://sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-file-datasets-summary

import os
import os.path
import sys
from oofemctrlreader import *
from FEM import *

class UNVParser:
    """ Universal file parser class"""
    def __init__(self,filename):
        self.file=None
        self.filename=filename
        self.FEM=FEM()
        self.startFlag='    -1'
        self.endFlag='    -1'
        # list of supported datasets and corresponding dataset handler functions
        self.datasetsIds=[2411,2412,2467, 2477]
        self.datasetsHandlers=[self.UNV2411Reader, self.UNV2412Reader, self.UNV2467Reader, self.UNV2467Reader]
        self.sections=[]

    def mapping(self):
        """Returns mapping for .unv elements"""
        #Table of element properties. It contains mapping of nodes, edges and faces between unv and OOFEM element.
    
        # TODO: Use a linked list where each oofem element is linked to the type of element and use the linked list when mapping occurs. In that way, we only need to specify each type of element (discretization) once.
    
        oofem_elemProp = []
        oofem_elemProp.append(oofem_elementProperties("None", [0], [], []))#leave this line [0] as it is
        oofem_elemProp.append(oofem_elementProperties("RepresentsBoundaryLoad", [],[],[]))#special element representing boundary load (edge or surface)
        oofem_elemProp.append(oofem_elementProperties("Truss1D", [0,1], [], []))
        oofem_elemProp.append(oofem_elementProperties("Interface1d", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Truss2D", [0,1], [0,1],[]))
        oofem_elemProp.append(oofem_elementProperties("Truss3D",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Beam2D",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LIBeam2D",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LIBeam2Dnl",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Beam3D",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LIBeam3D",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LIBeam3D2",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LIBeam3Dnl",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LIBeam3Dnl2",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("IntELPoint",oofem_elemProp[-1]))
        #oofem_elemProp.append(oofem_elementProperties("TrPlaneStress2D", [0,2,1], [[0,2],[2,1],[1,0]],[])) #checked - current numbering of triangle nodes is anti-clockwise, the same orientation as in OOFEM.
        oofem_elemProp.append(oofem_elementProperties("TrPlaneStress2D", [0,1,2], [[0,1],[1,2],[2,0]],[])) #old version of UNV export in SALOME, nodes on triangular elements are numbered clockwise.
        oofem_elemProp.append(oofem_elementProperties("TrPlaneStress2DXFEM", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("TrPlaneStrain",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Axisymm3D",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Tr1ht",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Tr1hmt",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Tr1mt",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Traxisym1ht",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("TrPlaneStrRot",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("CCTplate",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("CCTplate3D",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Tria1PlateSubSoil",oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("QTrPlStr", [2,0,4,1,5,3], [[2,1,0],[0,5,4],[4,3,2]],[]))#checked
        oofem_elemProp.append(oofem_elementProperties("QTrPlStrSlip", [2,0,4,1,5,3], [[2,1,0],[0,5,4],[4,3,2]],[]))
        oofem_elemProp.append(oofem_elementProperties("Tria2PlateSubSoil", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("PlaneStress2D", [0,1,2,3], [[0,1],[1,2],[2,3],[3,0]],[]))#checked
        oofem_elemProp.append(oofem_elementProperties("PlaneStress2DXFEM", [0,1,2,3], [[0,1],[1,2],[2,3],[3,0]],[]))#checked
        oofem_elemProp.append(oofem_elementProperties("Quad1PlaneStrain", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Quad1ht", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Quad1hmt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Quad1mt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Quadaxisym1ht", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Quadaxisym1hmt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Quadaxisym1mt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("quad1platesubsoil", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Interface2dlin", [0,1,3,2], [[0,1],[],[3,2],[]], []))
        oofem_elemProp.append(oofem_elementProperties("IntElLine1", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("QPlaneStress2D", [2,4,6,0,3,5,7,1], [[2,3,4],[4,5,6],[6,7,0],[0,1,2]],[]))#checked
        oofem_elemProp.append(oofem_elementProperties("QPlaneStress2DSlip", [2,4,6,0,3,5,7,1], [[2,3,4],[4,5,6],[6,7,0],[0,1,2]],[]))
        oofem_elemProp.append(oofem_elementProperties("QQuad1ht", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("QQuad1mt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("QQuad1hmt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Quad2plateSubsoil", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LSpace", [4,7,6,5,0,3,2,1], [[4,7],[7,6],[6,5],[5,4],[4,0],[7,3],[6,2],[5,1],[0,3],[3,2],[2,0],[1,0]], [[4,7,6,5],[0,3,2,1],[4,0,3,7],[7,3,2,6],[6,2,1,5],[5,1,0,4]]))#checked
        oofem_elemProp.append(oofem_elementProperties("Brick1ht", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Brick1mt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Brick1hmt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LSpaceBB", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("QSpace", [12,18,16,14,0,6,4,2,19,17,15,13,7,5,3,1,8,11,10,9], [[12,19,18],[18,17,16],[16,15,14],[14,13,12],[12,8,0],[18,11,6],[16,10,4],[14,9,2],[0,7,6],[6,5,4],[4,3,2],[2,1,0]], [[12,19,18,17,16,15,14,13],[0,7,6,5,4,3,2,1],[12,8,0,7,6,11,18,19],[18,11,6,5,4,10,16,17],[16,10,4,3,2,9,14,15],[14,9,2,1,0,8,12,13]])) #checked [brick nodes], [edges nodes], [faces nodes]
        oofem_elemProp.append(oofem_elementProperties("QBrick1ht", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("QBrick1hmt", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("LTRSpace", [0,1,2,3], [[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]], [[0,1,2],[0,1,3],[1,2,3],[0,2,3]]))#checked
        oofem_elemProp.append(oofem_elementProperties("Tetrah1ht", oofem_elemProp[-1]))
        oofem_elemProp.append(oofem_elementProperties("Tet1supg", [0,1,2,3], [[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]], [[0,1,2],[0,1,3],[1,2,3],[0,2,3]]))
        oofem_elemProp.append(oofem_elementProperties("tet21stokes", [9, 2, 0, 4, 7, 1, 6, 8, 3, 5], [], [[2,7,9,6,0,1],[2,3,4,8,9,7],[4,3,2,1,0,5],[0,6,9,8,5,4]]))
        oofem_elemProp.append(oofem_elementProperties("qtrspace", [9, 2, 0, 4, 7, 1, 6, 8, 3, 5], [], [[2,7,9,6,0,1],[2,3,4,8,9,7],[4,3,2,1,0,5],[0,6,9,8,5,4]]))
        oofem_elemProp.append(oofem_elementProperties("tet21ghostsolid", [9, 2, 0, 4, 7, 1, 6, 8, 3, 5], [], [[2,7,9,6,0,1],[2,3,4,8,9,7],[4,3,2,1,0,5],[0,6,9,8,5,4]]))
        oofem_elemProp.append(oofem_elementProperties("Tet1BubbleStokes", [0,1,2,3], [[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]], [[0,1,2],[0,1,3],[1,2,3],[0,2,3]]))
        oofem_elemProp.append(oofem_elementProperties("Qwedge", [0,2,4,9,11,13,1,3,5,10,12,14,6,7,8], [[0,1,2],[2,3,4],[4,5,0],[9,10,11],[11,12,13],[13,14,9],[0,6,9],[2,7,11],[4,8,13]], [[0,1,2,3,4,5],[9,10,11,12,13,14],[0,1,2,7,11,10,9,6],[2,3,4,8,13,12,11,7],[0,6,9,14,13,8,4,5]]))
        return oofem_elemProp


    def scanfile(self):
        """ Read file & fill the section list"""
        loop=True
        while loop:
            line=self.file.readline()
            if len(line)>0:
                if line.startswith(self.startFlag):
                    # identify section & save offset
                    id=int(self.file.readline())
                    offset=self.file.tell()
                    self.sections.append([id,offset])
                    # ignore data until end of section
                    while not(self.file.readline().startswith(self.endFlag)):
                        pass
            else:
                loop=False
        # rewind file
        self.file.seek(0)
        
    def UNV2411Reader(self, file, FEM):
        """ reads an UNV2411 dataset (nodes) from file and store data in FEM structure
            return an updated FEM object"""
        endFlag='    -1'
        loop=True
        while loop:
            line1=file.readline()
            line2=file.readline()
            if len(line2)>0:
                if line1.startswith(endFlag):
                    loop=False
                    break
                else:
                    dataline=Line2Int(line1)
                    coords=Line2Float(line2)
                    #print (dataline)
                    FEM.nodes.append(Node(dataline[0],coords))
                    FEM.nnodes=FEM.nnodes+1
            else:
                loop=False
        return FEM
    
    def UNV2412Reader(self, file, FEM):
        """ reads an UNV2412 dataset (elements) from file and store data in FEM structure
            return an updated FEM object"""
        endFlag='    -1'
        loop=True
        while loop:
            line1=file.readline()
            line2=file.readline()
            if len(line2)>0:
                if line1.startswith(endFlag):
                    loop=False
                    break
                else:
                    dataline=Line2Int(line1)
                    eltype=dataline[1]
                    if eltype==11 or eltype==22:# types of elements which are defined on 3 lines
                        # 1D elements have an additionnal line in their definition
                        line3=file.readline()
                        cntvt=Line2Int(line3)
                    elif eltype==113:#Quadratic wedge have nodes on 2 lines
                        line3=file.readline()
                        cntvt = Line2Int(line2) + Line2Int(line3)
                        print(cntvt, type(cntvt))
                    elif eltype==116:#Quadratic brick element have nodes on 3 lines
                        line3=file.readline()
                        line4=file.readline()
                        cntvt = Line2Int(line2) + Line2Int(line3) + Line2Int(line4)
                        #print(cntvt, type(cntvt))
                    elif eltype==118: # Quadratic tetrahedron has nodes on two lines
                        line3=file.readline()
                        cntvt=Line2Int(line2) + Line2Int(line3)
                        FEM.nodes[cntvt[9]-1].quadratic=0
                        FEM.nodes[cntvt[2]-1].quadratic=0
                        FEM.nodes[cntvt[0]-1].quadratic=0
                        FEM.nodes[cntvt[4]-1].quadratic=0
                        # 9, 2, 0, 4
                    else:
                        # standard elements have their connectivities on the second line
                        cntvt=Line2Int(line2)
                    if(len(dataline)<6):
                        print ("I need at least 6 entries on dataline %s" % dataline)
                        exit(0)
                    FEM.elems.append(Element(dataline[0],dataline[1],0,0,dataline[5],cntvt))
                    FEM.nelems=FEM.nelems+1
            else:
                loop=False
        return FEM
    
    def UNV2467Reader(self, file, FEM):
        """ reads an UNV2467 dataset (groups) from file and store data in FEM structure a group may represent a nodeset, an elementset or an edgeset ...
        return an updated FEM object"""
        endFlag='    -1'
        loop=True
        while loop:
            line1=file.readline()
            if line1 == '\n':
                continue
            if line1.startswith(endFlag):
                loop=False
                break
            else:
                # read group
                line2=file.readline()
                dataline=Line2Int(line1)
                groupname=line2
                if(len(dataline)==0):
                    print ("Group %s is empty, did you remesh the object and lost the members?" % groupname)
                    exit(0)
                else:
                    id=0 # dataline[0]
                    nitems=dataline[7]
                nlines=(nitems+1)//2
                # read group items
                lst=[]
                for i in range(nlines):
                    dat=Line2Int(file.readline())
                    #print "dat = ", dat
                    lst.append(dat[0:3])
                    if len(dat)>4:
                        lst.append(dat[4:7])
                # split group in node (7) sets, element or edge sets (8)
                nset=Group(id,groupname)
                elset=Group(id,groupname)
                nset.type=7
                elset.type=8
                for item in lst:
                    if item[0]==7:
                        nset.items.append(item[1])
                    if item[0]==8:
                        elset.items.append(item[1])
                nset.nitems=len(nset.items)
                elset.nitems=len(elset.items)
                # store non empty groups
                if nset.nitems>0:
                    FEM.nodesets.append(nset)
                if elset.nitems>0:
                    FEM.elemsets.append(elset)
                #print "%u \n" % elset.id
                FEM.nnodesets=len(FEM.nodesets)
                FEM.nelemsets=len(FEM.elemsets)
        return FEM

    def parse(self):
        """ parse UNV file to fill the FEM data structure"""
        self.file=open(self.filename,'r')
        self.scanfile()
        for sectionId,offset in self.sections:
            if (sectionId in self.datasetsIds):
                self.file.seek(offset)
                func=self.datasetsHandlers[self.datasetsIds.index(sectionId)]
                self.FEM=func(self.file,self.FEM)
        self.file.close()
        return self.FEM


if __name__=='__main__':
    helpmsg=""" UNV2X: example of use, write data in separate text files
        usage: UNV2X unvfile prefix
        what it does: read unvfile, create an internal FEM object structure
                      in memory and writes the following files:
                      prefix.nodes, prefix.elems, prefix.groups
    """
    if len(sys.argv)==3:
        unvfile=sys.argv[1]
        prefix=sys.argv[2]
        # read UNV file in FEM object structure
        UNV=UNVParser(unvfile)
        FEM=UNV.parse()
        # write files
        ls=os.linesep
        # node file
        nf=open(prefix + '.nodes','w')
        for node in FEM.nodes:
            nf.write(('%5d, %8g, %8g, %8g'+ls) % (node.id,node.coords[0],node.coords[1],node.coords[2]))
        nf.close()
        # element file
        ef=open(prefix + '.elems','w')
        for elem in FEM.elems:
            dat=[elem.id,elem.type,elem.material,elem.color,elem.nnodes]
            dat.extend([x for x in elem.cntvt])
            format='%5d, ' * (4 + elem.nnodes) + ' %5d' + ls
            ef.write(format % tuple(dat))
        ef.close()
        # group file
        grouptypes={7:'NodeSet',8:'ElementSet'}
        gf=open(prefix + '.groups','w')
        for group in (FEM.nodesets + FEM.elemsets):
            gf.write(('%s(%d) Id: %d'+ls) % (grouptypes[group.type],group.type,group.id))
            gf.write('Name: %s' % (group.name))
            count=0
            lst=group.items
            for i in range(group.nitems):
                count=count+1
                if (count<8)&(i!=g<roup.nitems-1):
                    gf.write('%5d, ' % lst.pop(0))
                else:
                    gf.write(('%5d'+ls) % lst.pop(0))
                    count=0
        gf.close()

    else:
        print(helpmsg)

