#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function
from unv2x import *
from abaqus2x import *
from oofemctrlreader import *
import time
import sys
import json
from numpy.core.defchararray import splitlines


if __name__=='__main__':
    helpmsg=""" 
Usage: unv2oofem.py unvfile ctrlfile oofemfile auxiliaryFile

What it does: read unvfile, create an internal FEM object structure
              in memory and writes the oofem native input file
              The ctrlfile specifies additional properties required by oofem

See http://www.oofem.org/wiki/doku.php?id=unv2oofem:unv2oofem for more info.

The format of ctrl file is following: (lines beginning with '#' are comments)

Output file record
Job description record
Analysis record
Domain record
Output manager record
ncrosssect # nmat # nbc # nic # nltf # nset # nxfemman #
cross section records
material records
boundary condition records
initial condition records
load time function records
extractor records
set records

Assignment of properties to nodes and elements is based on association with some unv group. The same mechanism
is valid for assignment of boundary conditions (edge, surface) load. The syntax is following:
group name1 [name2] [name3] ...
nodeprop "nodal_attributes_appended_to_nodal_records" [set INT]
elemprop "element_attributes_appended_to_element_records" [set INT]
etype[unv_etype] oofem_etype #provides mapping between unv and oofem element types

By default, all nodes will be exported,
elements are exported only when associated to some group
with valid element mapping

Enjoy.
"""
    welcomeMsg = """
UNV2OOFEM: Converts UNV file from Salome to OOFEM native file format
                    (C) 2009 Borek Patzak
                    Running python version %s.%s
""" % (sys.version_info.major, sys.version_info.minor)
    print (welcomeMsg)
    t1 = time.time()
    if (len(sys.argv)==4):
        unvfile=sys.argv[1]
        ctrlfile=sys.argv[2]
        oofemfile=sys.argv[3]
        of=open(oofemfile,'w')
        
        # read file in FEM object structure
        fileExtension = unvfile.split('.')
        if (fileExtension[-1].lower()=='unv'): # Salome output file
            Parser=UNVParser(unvfile)
        elif (fileExtension[-1].lower()=='inp'): # Abaqus output file
            Parser=AbaqusParser(unvfile)
        else:
            print ("Unknown extension of input file %s" % fileExtension[-1].lower())
            exit(0)
        
        print ('Parsing mesh file %s' % sys.argv[1], end=' ')
        FEM=Parser.parse()
        print ("done")

        print ("Detected node groups:", end=' ')
        for i in FEM.nodesets:
            print (i.name.strip(), end=' ')
        print ()

        print ("Detected element groups:", end=' ')
        for i in FEM.elemsets:
            print (i.name.strip(), end=' ')
        print ()

        # read oofem ctrl file
        CTRL=CTRLParser(ctrlfile, Parser.mapping())
        print ('Parsing ctrl file %s' % sys.argv[2])
        CTRL.parse(FEM)
        print ("done")
        # write files in native oofem format

        for sets in FEM.nodesets:
            if sets.name.startswith('OOFEM-Hanging'):
                print("Found group of HANGING nodes, assigning them")


        print ('Writing oofem file %s' % sys.argv[3])
        # write oofem header
        of.write(CTRL.header)

        #store elements in meshElements list. Reason: need to assign boundaryLoad to elements, which may be read after elements
        meshElements = []
        #create auxiliary array of element numbers to be searched for boundaryLoads
        elemNotBoundary = []
        
        # List for sets containing boundaries
        boundarySets=[];

        for elem in FEM.elems:#loop through all unv elements
            #resolve element properties
            properties=""
            for igroup in elem.oofem_groups:
                #print igroup.name
                properties+=igroup.oofem_properties
                #print('Properties', properties)
            #Do output if oofem_elemtype resolved and not BoundaryLoads
            if ( elem.oofem_elemtype):
                if(CTRL.oofem_elemProp[elem.oofem_elemtype].name != 'RepresentsBoundaryLoad'):
                    #Check if unv element and OOFEM element have the same amount of nodes
                    if (elem.nnodes != len(CTRL.oofem_elemProp[elem.oofem_elemtype].nodeMask)):
                        print ("\nUnv element #%d has %d nodes, which should be mapped on OOFEM element \"%s\" with %d nodes" % \
                            (elem.id, elem.nnodes,CTRL.oofem_elemProp[elem.oofem_elemtype].name, len(CTRL.oofem_elemProp[elem.oofem_elemtype].nodeMask)))
                        exit(0)

                    elemNotBoundary.append(elem)
                    dat = elem.oofem_outputData
                    dat.append(CTRL.oofem_elemProp[elem.oofem_elemtype].name)
                    dat.append("%-5d" % elem.id)
                    dat.append("nodes")
                    dat.append("%-3d" % elem.nnodes)
                    for n in range(elem.nnodes):
                        mask = CTRL.oofem_elemProp[elem.oofem_elemtype].nodeMask[n]
                        try:
                            dat.append("%-3d" % elem.cntvt[mask])
                        except:
                            print ("Exception in mapping nodes in unv element number %d, nodes %s" % (elem.id, elem.cntvt))
                            exit(0)
                    #dat.extend(["%-3d" % x for x in elem.cntvt])
                    dat.append(properties)
                    meshElements.append([])

        #Assign BoundaryLoads to elements (corresponds to edge and face loads).
        #We need to loop over all elements and to check whether they have assigned loads. This is quite time consuming but robust algorithm.
        for belem in FEM.elems:#loop over all elements from unv file
            #resolve element properties
            #for igroup in elem.oofem_groups:#unv element with boundary load is assigned to some ctrl element group
                #print belem.id, belem.oofem_elemtype, CTRL.oofem_elemProp[belem.oofem_elemtype].name
                if CTRL.oofem_elemProp[belem.oofem_elemtype].name == 'RepresentsBoundaryLoad':#found element, which represents boundary load
                    nodesOnBoundary = belem.cntvt
                    nodesOnBoundary.sort()
                    for elem in elemNotBoundary: #loop over, e.g. triangular elements, in order to find which element belem is a boundary to
                        cnt=0
                        for n in range(len(nodesOnBoundary)):
                            if(elem.cntvt.count(int(nodesOnBoundary[n]))):
                                cnt = cnt+1
                        if (cnt==len(nodesOnBoundary)):#found eligible element to which assign b.c. Now find which edge/face it is.
                            success = 0
                            if(belem.type==11 or belem.type==22):#elements representing EDGE loads
                                mask = CTRL.oofem_elemProp[elem.oofem_elemtype].edgeMask
                            else:#face loads
                                mask = CTRL.oofem_elemProp[elem.oofem_elemtype].faceMask

                            for i in range(len(mask)):
                                nodesInMask = []#list of nodes which are extracted according to mask
                                for x in mask[i]:
                                    nodesInMask.append(elem.cntvt[x])
                                #We need to compare both arrays nodesInMask and nodesOnBoundary. If they contain the same node numbers, we found edge/face.
                                nodesInMask.sort()
                                if(nodesInMask==nodesOnBoundary):#both lists are sorted so they can be compared
                                    success = 1
                                    #since boundary element may be in more unv groups, we need to find corresponding ctrl group
                                    for bel in belem.oofem_groups:
                                        #print "%d '%s' '%s'" % (len(belem.oofem_groups), bel.name.rstrip(), bel.oofem_groupNameForLoads)
                                        if (bel.name.rstrip() == bel.oofem_groupNameForLoads):
                                            #continue
                                            #build a new int list, which reflects load numbers and edges/faces
                                            if (len(bel.oofem_boundaryLoadsNum) > 0):
                                                loadNum = bel.oofem_boundaryLoadsNum
                                                newList=[-1]*(2*len(loadNum))
                                                for j in range(len(loadNum)):
                                                    newList[2*j] = loadNum[j]
                                                    newList[2*j+1] = i+1
                                                #print newList
                                                elem.oofem_bLoads+=newList
                                                print ("Boundary load \"%s\" found for element %d " % (bel.name.rstrip('\n'), elem.id))
                                                #print bel.name, elem.id, elem.oofem_bLoads
                                                
                                        if (bel.oofem_sets):
                                            print ("Set \"%s\" found for boundary of element %d " % (bel.name.rstrip('\n'), elem.id))
                                            setNum = bel.oofem_sets;
                                            # setID, element id, element side
                                            for thisSet in setNum:
                                                boundarySets.append([thisSet, elem.id, i+1])
                                                
                            if(success==0):
                                print ("Can not assign edge/face load \"%s\" to unv element %d" % (bel.name, elem.id))

        #write component record
        of.write('ndofman %d nelem %d ncrosssect %d nmat %d nbc %d nic %d nltf %d nset %d nxfemman %d\n' % (FEM.nnodes, len(elemNotBoundary), CTRL.ncrosssect, CTRL.nmat, CTRL.nbc, CTRL.nic, CTRL.nltf, CTRL.nset, CTRL.nxfemman))
        
        #write nodes
        for node in FEM.nodes:
            #resolve nodal properties
            hangNodes = {}
            for sets in FEM.nodesets:
                if sets.name.startswith('OOFEM-Hanging'):
                    for key in sets.items:
                        hangNodes[key] = ''
            
            outputLine="%s %-5d coords %-2d" % ('hangingNode' if node.id in hangNodes else 'node', node.id, len(node.coords))
            
            for coord in node.coords:
                outputLine+= "% -8g " % coord
            
            if node.id in hangNodes:
                outputLine+= "dofType %d %s" % (len(node.coords), ' 2 2 ' if len(node.coords) == 2 else ' 2 2 2 ') 
            
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
            #Add the list of boundaryLoads if it exists
            if(elem.oofem_bLoads):
                str+=" BoundaryLoads %d " % len(elem.oofem_bLoads)
                str+= ' '.join(["%d" % el for el in elem.oofem_bLoads])
            of.write('%s\n' % str)

        # write final sections
        sl=CTRL.footer.splitlines()
        for s in sl:
            words=s.split()
            #if len(words)==0:#skip empty lines
                #continue
            if (words[0].lower()=='set'):
                setID=int(words[1])

                if (words[2].lower()=='nodes'):
                    nodelist=[];
                    for nodeset in FEM.nodesets:
                        for oofemset in nodeset.oofem_sets:
                            if (setID==oofemset):
                                nodelist.extend(nodeset.items)
                    setElements=list(set(nodelist))

                elif (words[2].lower()=='elements'):
                    ellist=[]
                    for elemset in FEM.elemsets:
                        #print elemset.id
                        if setID == elemset.id:
                            ellist.extend(elemset.items)

                        for oofemset in elemset.oofem_sets:
                            if (setID==oofemset):
                                ellist.extend(elemset.items)
                    setElements=list(set(ellist))

                elif (words[2].lower()=='elementboundaries' or words[2].lower()=='elementedges'):
                    setElements=[]
                    for thisSet in boundarySets:
                        if (thisSet[0]==int(words[1])):
                            setElements.extend([thisSet[1], thisSet[2]])

                of.write('%s %s %s %u ' % ( words[0], words[1], words[2], len(setElements)) )
                for setElement in setElements:
                    of.write('%u ' % setElement)
                of.write('\n')

            else:
                of.write('%s\n' % s)

        of.close()
        #
        t2 = time.time()
        #
        print ("done ( %d nodes %d elements)" % (FEM.nnodes, len(elemNotBoundary)))
        print ("Finished in %0.2f [s]" % ((t2-t1)))

    else:
        print(helpmsg)



