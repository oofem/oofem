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
import os
import os.path
import sys

class FEM:
    """ a simple FEM object structure"""
    def __init__(self):
        self.nnodes=0
        self.nelems=0
        self.nnodesets=0
        self.nelemsets=0
        self.nodes=[]
        self.elems=[]
        self.nodesets=[]
        self.elemsets=[]
        
class Node:
    """ a single node object """
    def __init__(self,id,coords,quadratic=1):
        self.id=id
        self.coords=coords
        self.quadratic=quadratic
        
class Element:
    """ a single finite element object"""
    def __init__(self,id,type,material,color,nnodes,cntvt):
        self.id=id
        self.type=type
        self.material=material
        self.color=color
        self.nnodes=nnodes
        self.cntvt=cntvt

class Group:
    """ FE group object """
    def __init__(self,id,name):
        self.id=id
        self.type=0
        self.name=name
        self.nitems=0
        self.items=[]

def Line2Float(line):
    """Convert a string into a list of Float"""
    return map(float,line.split())

def Line2Int(line):
    """Convert a string into a list of Int"""
    return map(int,line.split())

def UNV2411Reader(file,FEM):
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
                FEM.nodes.append(Node(dataline[0],coords))
                FEM.nnodes=FEM.nnodes+1
        else:
            loop=False
    return FEM

def UNV2412Reader(file,FEM):
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
                elif eltype==118: # Quadratic tetrahedron has nodes on two lines
                    line3=file.readline()
                    cntvt=Line2Int(line2) + Line2Int(line3)
                    FEM.nodes[cntvt[9]-1].quadratic=0
                    FEM.nodes[cntvt[2]-1].quadratic=0
                    FEM.nodes[cntvt[0]-1].quadratic=0
                    FEM.nodes[cntvt[4]-1].quadratic=0
                    # 9, 2, 0, 4
                elif eltype==116:#Quadratic brick element has data on 4 lines
                    line3=file.readline()
                    line4=file.readline()
                    cntvt = Line2Int(line2) + Line2Int(line3) + Line2Int(line4)
                    #print cntvt, type(cntvt)
                else:
                    # standard elements have their connectivities on second line
                    cntvt=Line2Int(line2)
                if(len(dataline)<6):
                    print "I need at least 6 entries on dataline %s" % dataline
                    exit(0)
                FEM.elems.append(Element(dataline[0],dataline[1],0,0,dataline[5],cntvt))
                FEM.nelems=FEM.nelems+1
        else:
            loop=False
    return FEM

def UNV2467Reader(file,FEM):
    """ reads an UNV2467 dataset (groups) from file and store data in FEM structure a group may represent a nodeset, an elementset or an edgeset ...
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
                # read group
                dataline=Line2Int(line1)
                groupname=line2
                if(len(dataline)==0):
                    print "Group %s is empty, did you remesh the object and lost the members?" % groupname
                    exit(0)
                id=dataline[0]
                nitems=dataline[7]
                nlines=(nitems+1)/2
                # read group items
                lst=[]
                for i in range(nlines):
                    dat=Line2Int(file.readline())
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
                FEM.nnodesets=len(FEM.nodesets)
                FEM.nelemsets=len(FEM.elemsets)
        else:
            loop=False
    return FEM


class UNVParser:
    """ Universal file parser class"""
    def __init__(self,filename):
        self.file=None
        self.filename=filename
        self.FEM=FEM()
        self.startFlag='    -1'
        self.endFlag='    -1'
        # list of supported datasets and corresponding dataset handler functions
        self.datasetsIds=[2411,2412,2467]
        self.datasetsHandlers=[UNV2411Reader,UNV2412Reader,UNV2467Reader]
        self.sections=[]

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

    def parse(self):
        """ parse UNV file to fill the FEM data structure"""
        self.file=open(self.filename,'rb')
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
                if (count<8)&(i<>group.nitems-1):
                    gf.write('%5d, ' % lst.pop(0))
                else:
                    gf.write(('%5d'+ls) % lst.pop(0))
                    count=0
        gf.close()

    else:
        print(helpmsg)


