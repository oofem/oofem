# -*- coding: utf-8 -*-
# OOFEM CTRL FILE reader, B.Patzak, www.oofem.org, (c) 2009
# Licence: GPL
# Purpose: this is a set of objects & functions to read oofem ctrl File
#          into a simple FEM object structure in order to simplify the
#          conversion of Mesh definitions from UNV to OOFEM native format
# How it works:
#          This code is based on two main objects:
#           1) a simple FEM object structure storing nodes, elements and groups
#              as provided bu UNVParser from unv2x module 
#           2) CTRLParser reads and parses oofem ctrl file, extends the data structure in given FEM 
#              object to provide a simple & modular solution to read
#              some datasets from UNV file and store them in OOFEM data structure
#           3) It adds following records into FEM data structure:
#              FEM.Node.oofem_groups       - list of receiver groups
#              FEM.Element.oofem_groups    - list of receiver groups
#              FEM.Element.oofem_elemtype  - oofem element type, ID in the list of elements
#              FEM.Group.oofem_properties  - attributes related to components in that group
#                                            these are resolved to individual components
#                                            in CTRL.parse
#              FEM.Group.oofem_etypemap    - dictionary of element types
#                                            key is unv element id
#                                            value is corresponding oofem element type
#                                            if element mapping is not provided, element is skipped


import os
import os.path
import sys
import re

def Line2Float(line):
    """Convert a string into a list of Float"""
    return map(float,line.split())

def Line2Int(line):
    """Convert a string into a list of Int"""
    return map(int,line.split())

#Structure holding OOFEM's element characteristics
class oofem_elementProperties:
    def __init__(self,name,nodeMask,edgeMask,faceMask):
        self.name=name#string - name if OOFEM elements
        self.nodeMask=nodeMask#list of nodes in the sequence of OOFEM
        self.edgeMask=edgeMask#2D array expressing node masks for OOFEM's edge 1,2,..., original UNV node numbering
        self.faceMask=faceMask#2D array expressing node masks for OOFEM's face 1,2,..., original UNV node numbering

class CTRLParser:
    """ a simple CTRL object structure"""

#   Table of element characteristics
    oofem_elemProp = []
    oofem_elemProp.append(oofem_elementProperties("None", [0], [], []))#leave this line [0] as it is
    oofem_elemProp.append(oofem_elementProperties("BoundaryLoads", [0,1],[],[]))#Edge loads are treated as elements in UNV
    oofem_elemProp.append(oofem_elementProperties("Tr1ht", [0,2,1], [[0,2],[2,1],[1,0]],[]))
    oofem_elemProp.append(oofem_elementProperties("TrPlaneStress2d", [0,2,1], [[0,2],[2,1],[1,0]],[]))

    def __init__(self, filename):
        self.filename=filename
        self.file=None
        self.header=""
        self.footer=""
        self.ncrosssect=0
        self.nmat=0
        self.nbc=0
        self.nic=0
        self.nltf=0

    def getRecordLine(self):
        while True:
            line=self.file.readline()
            if not line.startswith('#'):
                break
        return line

    def getNodeGroup(self, FEM, gid):
        for i in FEM.nodesets:
            #print i.id, i.name, i.items
            if (i.name.strip() == gid):
                return i
        return None

    def getElementGroup(self, FEM, gid):
        for i in FEM.elemsets:
            #print i.id, i.name, i.items
            if (i.name.strip() == gid):
                return i
        return None

    def addGroupToComponent(self, comp, gid):
            comp.oofem_groups.append(gid)

    def parseGroup(self,FEM):
        #loop over lines until EOF
        while True:
            line=self.file.readline()
            if not line: break
            if line.startswith('#'):
                continue
            match=re.search('^group((\w| )+)', line)
            if match:
                groups = match.group(1).split()
                print "\tFound properties for group(s):", groups
                # parse group record
                while True:
                    line=self.file.readline()
                    if not line: 
                        break
                    match=re.search('^(nodeprop|elemprop|etype\[\d+\]|#)\s+((\w| |\t)+)',line)
                    if match:
                        if match.group(1) == 'nodeprop':
                            for igroup in groups:
                                __gr=self.getNodeGroup(FEM,igroup)
                                if __gr:
                                    __gr.oofem_properties=match.group(2)
                                    print "\t\tgroup \"",igroup, "\" nodeprops:", __gr.oofem_properties
                                else:
                                    print "\t\tWARNING: group \"",igroup, "\" no such node group found"
                        elif match.group(1) == 'elemprop':
                            for igroup in groups:
                                __gr=self.getElementGroup(FEM,igroup)
                                if __gr:
                                    __gr.oofem_properties=match.group(2)
                                    print "\t\tgroup \"",igroup, "\" elemprops:", __gr.oofem_properties
                                else:
                                    print "\t\tWARNING: group \"",igroup, "\" no such element group found"
                        elif match.group(1)[:5] == 'etype':
                            etmatch=re.search('etype\[(\d+)\]*', match.group(1))
                            unvetype = int(etmatch.group(1))
                            for igroup in groups:
                                __gr=self.getElementGroup(FEM,igroup)
                                if __gr:
                                    elemName = match.group(2).strip()
                                    #check that elemName exists in a list and assign
                                    for n in range(len(self.oofem_elemProp)):
                                        if(self.oofem_elemProp[n].name == elemName):
                                            __gr.oofem_etypemap[unvetype]= n
                                            break
                                    else:
                                        print "OOFEM element %s not found in OOFEM's list of eligible elements" % elemName
                                        sys.exit(0)

                                    print "\t\tgroup \"",igroup, "\" etype[", unvetype, "] =", __gr.oofem_etypemap[unvetype]
                                else:
                                    print "\t\tWARNING: group \"",igroup, "\" no such element group found"
                    else:
                        break


    def parse(self, FEM):

        self.file=open(self.filename,'rb')#'rb' mode for M$ compatibility - tell() and seek() functions
        # read header info ending with OutputManager
        while True:
            line=self.file.readline()
            if len(line)==1:
                print "I did not find keyword \"OutputManager\" or empty line is present in the header"
                sys.exit(0)
            self.header+=line
            if (line.split()[0].lower()=="outputmanager"):
                break

        # read component info (ncrossSect, nMat, mBc, nIc, nLtf) and extractor information to the end of the file
        data=self.getRecordLine()
        dataline = data.split()
        if(dataline[0].lower()!="ncrosssect" or dataline[2].lower()!="nmat" or dataline[4].lower()!="nbc" or dataline[6].lower()!="nic" or dataline[8].lower()!="nltf"):
            print "Error in entry: %s" % data
            print "Correct format: ncrosssect 1 nmat 8 nbc 2 nic 1 nltfs 1"
            sys.exit(0)
        self.ncrosssect=int(dataline[1])
        self.nmat=int(dataline[3])
        self.nbc=int(dataline[5])
        self.nic=int(dataline[7])
        self.nltf=int(dataline[9])

        #read crossSect, material, bc, ic, and lft records into footer
        count = 0
        while count < (self.ncrosssect+self.nmat+self.nbc+self.nic+self.nltf):
            line=self.file.readline()
            self.footer+= line
            if not line.startswith('#'):
                count+= 1

        #look, whether the next line contains extractor data
        pos = self.file.tell()
        line=self.file.readline()
        self.file.seek(pos,0)#return one line back
        if (len(line)>1 and line.split()[0] == '#%BEGIN_CHECK%'):
            while (line.split()[0] != '#%END_CHECK%'):
                line=self.file.readline()
                self.footer+= line
                if (len(line)==0 or len(line)==1):
                    print "EOF reached prematurely, missing #%END_CHECK%"
                    sys.exit(0)

        #init group properties
        for igroup in FEM.nodesets:
            igroup.oofem_properties=""
        for igroup in FEM.elemsets:
            igroup.oofem_properties=""
            igroup.oofem_elemtype=0
            igroup.oofem_etypemap={}
        # read and parse individual group records
        while True:
            if not self.parseGroup (FEM):
                break

        # build component maps and apply group info to components
        nodemap=dict([(x.id, x) for x in FEM.nodes])
        elemmap=dict([(x.id, x) for x in FEM.elems])
        #init component properties
        for i in FEM.nodes:
            i.oofem_groups=[]
        for i in FEM.elems:
            i.oofem_groups=[]
            i.oofem_elemtype=0
            i.oofem_outputData=[]#used for outputting OOFEM file

        # apply group membership to individual components
        for igroup in FEM.nodesets:
            for inode in igroup.items:
                self.addGroupToComponent(nodemap[inode], igroup)
        for igroup in FEM.elemsets:
            for ielem in igroup.items:
                self.addGroupToComponent(elemmap[ielem], igroup)
                #resolve elemtype
                if elemmap[ielem].type in igroup.oofem_etypemap:
                    elemmap[ielem].oofem_elemtype = igroup.oofem_etypemap[elemmap[ielem].type]


