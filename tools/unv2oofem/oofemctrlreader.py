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
#                                            in CTRL parse
#                                            set INT can be specified as the last argument
#              FEM.Group.oofem_etypemap    - dictionary of element types
#                                            key is unv element id
#                                            value is corresponding oofem element type
#                                            if element mapping is not provided, element is skipped


import os
import os.path
import sys
import re

def remove_values_from_list(the_list, val):
    the_list = [element.lower() for element in the_list]
    val = val.lower()
    while val in the_list:
        the_list.remove(val)

def Line2Float(line):
    """Convert a string into a list of Float"""
    return map(float,line.split())

def Line2Int(line):
    """Convert a string into a list of Int"""
    return map(int,line.split())

#Structure holding OOFEM's element properties
class oofem_elementProperties:
    def __init__(self,name,nodeMaskOrPointer,edgeMask=None,faceMask=None):
        if (edgeMask!=None):
            self.name=name#string - name as OOFEM elements
            self.nodeMask=nodeMaskOrPointer#list of nodes in the sequence of OOFEM
            self.edgeMask=edgeMask#2D array expressing node masks for OOFEM's edge 1,2,..., original UNV node numbering
            self.faceMask=faceMask#2D array expressing node masks for OOFEM's face 1,2,..., original UNV node numbering
        else:
            self.name=name#string - name as OOFEM elements
            self.nodeMask=nodeMaskOrPointer.nodeMask#list of nodes in the sequence of OOFEM
            self.edgeMask=nodeMaskOrPointer.edgeMask#2D array expressing node masks for OOFEM's edge 1,2,..., original UNV node numbering
            self.faceMask=nodeMaskOrPointer.faceMask#2D array expressing node masks for OOFEM's face 1,2,..., original UNV node numbering

class CTRLParser:
    """ a simple CTRL object structure"""
    #oofem_elemProp = []
    #oofem_elemProp.append(oofem_elementProperties("None", [0], [], []))#leave this line [0] as it is
    # Node assembly
    
    def __init__(self, filename, mapping):
        self.filename=filename
        self.file=None
        self.header=""
        self.footer=""
        self.ncrosssect=0
        self.nmat=0
        self.nbc=0
        self.nic=0
        self.nltf=0
        self.nset=0
        self.wplist=[[], [], []] # wp, element, edge
        self.oofem_elemProp=mapping
        self.nxfemman=0

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
            match=re.search('^group(([\-\w]| )+)', line, re.IGNORECASE)
            if match:
                groups = match.group(1).split()
                print ("\tFound properties for group(s):", groups)
                # parse group record
                while True:
                    line=self.file.readline()
                    if not line:
                        break
                    lineSplit = line.split()
                    if len(lineSplit) == 0:
                        break;
                    #print("Line: ", line)
                    if lineSplit[0].lower() == 'nodeprop':
                        if (lineSplit[-2].lower()=='set'):
                            for igroup in groups:
                                __gr=self.getNodeGroup(FEM, igroup)
                                if __gr:
                                    __gr.oofem_sets.append(int(lineSplit[-1]));
                                    str = "\t\tGroup of nodes belongs to set %d" % int(lineSplit[-1])
                                    print (str)
                                    if (len(lineSplit)>3):#contains also other properties
                                        __gr=self.getNodeGroup(FEM, igroup)
                                        __gr.oofem_properties=' '.join(lineSplit[1:-2])
                                        str = "\t\tGroup of nodes \"%s\" has properties: %s" % (igroup, __gr.oofem_properties)
                                        if igroup.startswith('OOFEM-Hanging'):
                                            __gr.oofem_hangingNode = True
                                            str += ' and are HANGING nodes'
                                        print (str)
                                else:
                                    str = "WARNING: Group of nodes \"%s\" does not exist" % igroup
                                    print (str)
                        else:
                            for igroup in groups:
                                #print("igroup", igroup)
                                __gr=self.getNodeGroup(FEM,igroup)
                                if __gr:
                                    __gr.oofem_properties=' '.join(lineSplit[1:])
                                    str = "\t\tGroup of nodes \"%s\" has properties: %s" % (igroup, __gr.oofem_properties)
                                    if igroup.startswith('OOFEM-Hanging'):
                                        __gr.oofem_hangingNode = True
                                        str += ' and are HANGING nodes'
                                    print (str)
                                else:
                                    str = "WARNING: Group of nodes \"%s\" does not exist" % igroup
                                    print (str)
                    elif lineSplit[0].lower() == 'elemprop':
                        if (lineSplit[1].lower()=='bloadnum'):#check if the group represents a boundary load
                            for igroup in groups:
                                __gr=self.getElementGroup(FEM,igroup)
                                if __gr:
                                    __gr.oofem_boundaryLoadsNum=[int(s) for s in lineSplit[2:]]
                                    __gr.oofem_groupNameForLoads=igroup.lstrip()
                                    str = "\t\tGroup of elements \"%s\" has boundary loads with numbers: %s" % (igroup, __gr.oofem_boundaryLoadsNum)
                                    print (str)
                                else:
                                    str = "WARNING: Group of elements \"%s\" for boundary load does no exist" % (igroup)
                        elif (lineSplit[-2].lower()=='set'):#check if the group represents a set
                            for igroup in groups:
                                __gr=self.getElementGroup(FEM,igroup)
                                if __gr:
                                    # Save id of set
                                    __gr.oofem_sets.append(int(lineSplit[-1]))
                                    # Save orientation
                                    str = "\t\tGroup of elements \"%s\" belongs to set: %s" % (igroup, lineSplit[-1])
                                    print (str)
                                    if (len(lineSplit)>3):#contains also other properties
                                        __gr=self.getElementGroup(FEM,igroup)
                                        __gr.oofem_properties=' '.join(lineSplit[1:-2])
                                        str = "\t\tGroup of elements \"%s\" has properties: %s" % (igroup, __gr.oofem_properties)
                                        print (str)
                                else:
                                    str = "WARNING: Group of elements \"%s\" for set does no exist" % (igroup)
                                    print (str)
                        else:#not boundary loads neither set
                            for igroup in groups:
                                __gr=self.getElementGroup(FEM,igroup)
                                if __gr:
                                    __gr.oofem_properties=' '.join(lineSplit[1:])
                                    str = "\t\tGroup of elements \"%s\" has properties: %s" % (igroup, __gr.oofem_properties)
                                    print (str)
                                else:
                                    str = "WARNING: Group of elements \"%s\" does no exist" % (igroup)
                                    print (str)
                    elif lineSplit[0][:5].lower() == 'etype':
                        etmatch=re.search('^etype\[(\d+)\]', lineSplit[0], re.IGNORECASE)
                        unvetype = int(etmatch.group(1))
                        for igroup in groups:
                            __gr=self.getElementGroup(FEM,igroup)
                            if __gr:
                                if len(lineSplit) > 1:
                                    elemName = lineSplit[1].strip().rstrip('\n')
                                else:
                                    if (__gr.oofem_boundaryLoadsNum):
                                        elemName = 'RepresentsBoundaryLoad'#assign this name to an element so we know it represents a boundary load
                                    elif (__gr.oofem_sets):
                                        elemName = 'RepresentsBoundaryLoad'
                                    else:
                                        elemName = 'Unknown'
                                #check that elemName exists in a list and assign
                                for n in range(len(self.oofem_elemProp)):
                                    if(self.oofem_elemProp[n].name.lower() == elemName.lower()):
                                        __gr.oofem_etypemap[unvetype]= n
                                        break
                                else:
                                    print ("OOFEM element %s not found in OOFEM's list of eligible elements" % elemName)
                                    sys.exit(0)

                                str = "\t\tGroup of elements \"%s\" of unv_element_type[%d] = %s" % (igroup, unvetype, elemName)
                                print (str)
                            else:
                                str = "WARNING: Group of elements \"%s\" not found" % (igroup)
                                print (str)
                    elif lineSplit[0].lower() == '#':
                        continue;

                else:
                    break


    def parse(self, FEM):

        self.file=open(self.filename,'r')#'rb' mode for M$ compatibility
        # read header info ending with OutputManager
        while True:
            line=self.file.readline()
            if len(line)==1:
                print ("I did not find keyword \"OutputManager\" or empty line is present in the header")
                sys.exit(0)
            self.header+=line
            if (line.split()[0].lower()=="outputmanager"):
                break

        # read component info (ncrossSect, nMat, nBc, nIc, nLtf) and extractor information to the end of the file
        data=self.getRecordLine()
        dataline = data.split()
        if(dataline[0].lower()!="ncrosssect" or dataline[2].lower()!="nmat" or dataline[4].lower()!="nbc" or dataline[6].lower()!="nic" or dataline[8].lower()!="nltf"):
            print ("Error in entry: %s" % data)
            print ("Correct format: ncrosssect # nmat # nbc # nic # nltf #")
            print ("but keywords are %s # %s # %s # %s # %s"%(dataline[0], dataline[2], dataline[4], dataline[6], dataline[8]))
            sys.exit(0)

        if len(dataline)>10 and dataline[10].lower() == "nset":
            self.nset = int(dataline[11])
        self.ncrosssect=int(dataline[1])
        self.nmat=int(dataline[3])
        self.nbc=int(dataline[5])
        self.nic=int(dataline[7])
        self.nltf=int(dataline[9])
        
        if dataline.count("nxfemman"):
            xmanInd = dataline.index("nxfemman")
            if xmanInd != -1:
                self.nxfemman = int(dataline[xmanInd+1])

        #read crossSect, material, bc, ic, and lft records into footer
        count = 0
         
        while count < (self.ncrosssect+self.nmat+self.nbc+self.nic+self.nltf):
            line=self.file.readline()
            dataline=line.split()
            
            self.footer+= line
            if not line.startswith('#'):
                count+= 1

        #attach all following lines until empty line is hit. This is the end of OOFEM section
        while True:
            line=self.file.readline()
            dataline=line.split()
            
            if line.strip() == '':
                break
            
            #No need to increase nset
            #if dataline[0].lower() == 'set':
            #    self.nset=self.nset + 1;
            self.footer+= line

        #look, whether the next line contains extractor data
        #pos = self.file.tell()
        #line=self.file.readline()
        #self.file.seek(pos,0)#return one line back
        #if (len(line)>1 and line.split()[0] == '#%BEGIN_CHECK%'):
            #while (line.split()[0] != '#%END_CHECK%'):
                #line=self.file.readline()
                #self.footer+= line
                #if (len(line)==0 or len(line)==1):
                    #print "EOF reached prematurely, missing #%END_CHECK%"
                    #sys.exit(0)

        #init group properties on UNV data
        for igroup in FEM.nodesets:#defined in ctrl file
            igroup.oofem_properties=""
            igroup.oofem_hangingNode = False
            igroup.oofem_sets=[]
        for igroup in FEM.elemsets:#defined in ctrl file
            igroup.oofem_properties=""
            igroup.oofem_elemtype=0
            igroup.oofem_etypemap={}
            igroup.oofem_boundaryLoadsNum=[]#numbers of boundary loads
            igroup.oofem_sets=[]
            igroup.oofem_groupNameForLoads=""#CTRL-group name
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
            i.oofem_bLoads=[]#array storing BoundaryLoads IDs, they will be merged at the output routine

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
                    
