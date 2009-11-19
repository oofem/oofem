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
#              FEM.Element.oofem_elemtype  - oofem element type 
#              FEM.Group.oofem_properties  - attributes related to components in that group
#                                            these are resolved to individual components
#                                            in CTRL.parse
#              FEM.Group.oofem_etypemap    - ddictionary of element types
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



class CTRLParser:
    """ a simple CTRL object structure"""

#   table of element names 
    elementNames=("None",
                  "TrPlaneStress2d");


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
                                    __gr.oofem_etypemap[unvetype]=int(match.group(2).strip())
                                    print "\t\tgroup \"",igroup, "\" etype[", unvetype, "] =", __gr.oofem_etypemap[unvetype]
                                else:
                                    print "\t\tWARNING: group \"",igroup, "\" no such element group found"
                    else:
                        break
                

                                
    def parse(self, FEM):

        self.file=open(self.filename,'r')
        count = 0
        # read header info
        while count < 5:
            line=self.file.readline()
            self.header+=line
            if not line.startswith('#'):
                count+= 1
        # read component info (ncrossSect, nMat, mBc, nIc, nLtf)
        
        dataline=Line2Int(self.getRecordLine())
        self.ncrosssect=dataline[0]
        self.nmat=dataline[1]
        self.nbc=dataline[2]
        self.nic=dataline[3]
        self.nltf=dataline[4]

        #read crossSect,material, bc, ic, and lft records into footer
        count = 0
        while count < (self.ncrosssect+self.nmat+self.nbc+self.nic+self.nltf):
            line=self.file.readline()
            self.footer+= line
            if not line.startswith('#'):
                count+= 1

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
