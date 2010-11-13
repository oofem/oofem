# -*- coding: utf-8 -*-
from unv2x import *
from oofemctrlreader import *
import time

if __name__=='__main__':
    helpmsg=""" 
Usage: unv2oofem.py unvfile ctrlfile oofemfile

What it does: read unvfile, create an internal FEM object structure
              in memory and writes the oofem native input file
              The ctrlfile specifies additional properties required by oofem

The format of ctrl file is following: (lines beginning with '#' are comments)

Output file record
Job description record
Analysis record
Domain record
Output manager record
# number of cross section records (ncrosssect), material records (nmat)
# number of initial conditions (nic), and number of load time functions (nltf)
ncrosssect nmat nbc nic nltf
# individual records
cross section records
material records
boundary condition records
initial condition records
load time function records

#The attributes to nodes and elements
#can be assigned on group by group basis
#for each group one can specify nodal and element
#properties, as well as specify the element types
#the first line specifies the group names for which the attributes are provided
#following lines specify individual attributes
#no empty line(s) between attribute records
#The syntax is following:
group name1 [name2] [name3] ...
[nodeprop "nodal_attributes_appended_to_nodal_records"]
[elemprop "element_attributes_appended_to_element_records"]
[etype[unv_etype] oofem_etype] #provides mapping between unv and oofem element types 

#The group sections are separated by empty line(s).


By default, all nodes will be exported, 
elements are exported only when associated to some group
with valid element mapping


Enjoy. 
"""
    print """
UNV2OOFEM: Converts UNV file from Salome to OOFEM native file format
                    (C) 2009 Borek Patzak
"""
    t1 = time.time()
    if len(sys.argv)==4:
        unvfile=sys.argv[1]
        ctrlfile=sys.argv[2]
        oofemfile=sys.argv[3]
        of=open(oofemfile,'w') 
        # read UNV file in FEM object structure
        UNV=UNVParser(unvfile)
        print 'Parsing unv file ....',
        FEM=UNV.parse()
        print "done"
        
        print "Detected node groups:",
        for i in FEM.nodesets:
            print i.name.strip(),
        print

        print "Detected element groups:",
        for i in FEM.elemsets:
            print i.name.strip(),
        print

        # read oofem ctrl file
        CTRL=CTRLParser(ctrlfile)
        print 'Parsing ctrl file ....'
        CTRL.parse(FEM)
        print "done"
        # write files in native oofem format
        
        print 'Writting oofem file ...',
        # write oofem header
        of.write(CTRL.header)

        #store elements in meshElements list. Reason: need to assign boundaryLoad to elements, which may be read after elements
        meshElements = []
        #create auxiliary array of elemtn numbers to be searched for boundaryLoads
        elemNotBoundary = []
        #counter = 0
        for elem in FEM.elems:
            #resolve element properties
            properties=""
            for igroup in elem.oofem_groups:
                properties+=igroup.oofem_properties
            #do output if oofem_elemtype resolved and not BoundaryLoads
            if elem.oofem_elemtype and CTRL.oofem_elemProp[elem.oofem_elemtype].name != 'BoundaryLoads':
                elemNotBoundary.append(elem)
                dat = elem.oofem_outputData
                dat.append(CTRL.oofem_elemProp[elem.oofem_elemtype].name)
                dat.append("%-5d" % elem.id)
                dat.append("nodes")
                dat.append("%-3d" % elem.nnodes)
                for n in range(elem.nnodes):
                    mask = CTRL.oofem_elemProp[elem.oofem_elemtype].nodeMask[n]
                    dat.append("%-3d" % elem.cntvt[mask])
                #dat.extend(["%-3d" % x for x in elem.cntvt])
                dat.append(properties)
                meshElements.append([])
        
        #assign BoundaryLoads to elements
        nboLoads = 0
        for belem in FEM.elems:
            #resolve element properties
            properties=""
            for igroup in belem.oofem_groups:
                properties+=igroup.oofem_properties
            if CTRL.oofem_elemProp[belem.oofem_elemtype].name == 'BoundaryLoads':
                nodesOnBoundary = belem.cntvt
                for elem in elemNotBoundary:
                    cnt=0
                    for n in range(len(nodesOnBoundary)):
                        if(elem.cntvt.count(int(nodesOnBoundary[n]))):
                            cnt = cnt+1
                    if (cnt==len(nodesOnBoundary)):#found eligible element
                        edgeMask = CTRL.oofem_elemProp[elem.oofem_elemtype].edgeMask
                        for i in range(len(edgeMask)):
                            nodesInMask = []#list of nodes which are extracted according to mask
                            for x in edgeMask[i]:
                                nodesInMask.append(elem.cntvt[x])
                            if (len(edgeMask[i]) == 2):#extract edges defined by two points
                                if(nodesOnBoundary.count(nodesInMask[0]) and nodesOnBoundary.count(nodesInMask[1])):
                                    elem.oofem_outputData.append(properties)
                                    elem.oofem_outputData.append("%d" % (i+1))
       
        #write component record
        of.write('ndofman %d nelem %d ncrosssect %d nmat %d nbc %d nic %d nltf %d\n' % (FEM.nnodes, len(elemNotBoundary), CTRL.ncrosssect, CTRL.nmat, CTRL.nbc, CTRL.nic, CTRL.nltf))
        #write nodes
        for node in FEM.nodes:
            #resolve nodal properties
            outputLine="node %-5d coords %-2d" % (node.id, len(node.coords))
            for coord in node.coords:
                outputLine+= "%-8g " % coord
            properties=""
            for igroup in node.oofem_groups:
                if(len(properties)>0 and properties[-1]!=" "):#insert white space if necessary
                    properties+=" "
                properties+=igroup.oofem_properties
            outputLine+=properties
            # write nodal record
            of.write(('%s\n') % (outputLine))

        for elem in elemNotBoundary:
            str = ' '.join(elem.oofem_outputData)
            of.write('%s\n' % str) 

        # write final sections
        of.write(CTRL.footer);
        of.close()
        #
        t2 = time.time()
        #
        print "done ( %d nodes %d elements)" % (FEM.nnodes, len(elemNotBoundary))
        print "Finished in %0.2f [s]" % ((t2-t1))
        
    else:
        print(helpmsg)

        


