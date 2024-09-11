#!/usr/bin/python
#
# oofem2part.py       (c) 2013 Borek Patzak, www.oofem.org
#
# Description
#    
#     oofem2part reads serial oofem input file and performs 
#     the node-cut partitioning of the problem domain 
#     into a set of subdomains (using metis),
#     automatically identifying shared nodes.
#     On output, it produces a set of oofem input files for each 
#     subproblem, suitable for parallel oofem analysis.
#
# Dependencies
#
#     oofem2part requires "Metis to python" wrapper module
#     (see http://metis.readthedocs.org/en/latest/)
#     the METIS_DLL environment variable should point to 
#     Metis shared library location.
#
# Usage
#
#     oofem2part -f file -n #
#     where -f option allows to set path to oofem input input file, 
#           -n option sets the number of desired target partitions
#
# Note
#
#     Some adjustment of created parallel input files may be necessary,
#     typically the solver has to be changed into parallel one.
#

import sys
import re
#from sets import Set
import getopt

#debug flag
debug = 0

#some constants defining nodal status
UNDEF = 0
LOCAL = 1
SHARED= 2


#----do not change------
# number of partitions
nparts = 0
# input File Name 
inputFileName = ''

def print_usage():
    print ("\nUsage:\noofem2part -f file -n #")
    print ("where -f file sets path to input oofem (serial) file to be partitioned")
    print ("      -n #    allows to se the required number of target partitions")


# returns the value corresponding to given keyword and record
def getKeywordValue (record, kwd, optional = None):
    match = re.search (kwd+'\s+\"*([\\\.\+\-,:\w]+)\"*', record, re.IGNORECASE)
    if match:
        return match.group(1)
    else:
        #issue an error if argument compulsory
        if optional == None:
            print ("\nMissing keyword \"", kwd, "\" in\n", record)
            exit (1)
        else:
            return optional


def readRecord (file):
    line = file.readline()
    while (re.search (r"^#", line)):
        line = file.readline()

    #print line
    return line


# this pattern covers simple slave relations and rigid arm nodes
mastermask_re = re.compile(r'\s+mastermask\s+(\d+)\s+((\d+\s+)+)', re.IGNORECASE)
slavenodemasters_re = re.compile(r'^slavenode\s+.*\s+masterdofman\s+(\d+)\s+((\d+\s+)+)', re.IGNORECASE)
hangingnode_re = re.compile(r'^hangingnode', re.IGNORECASE)



# returns a list of node masters, empty list if node is not slave
def giveNodeMasters(inode):
    global nodes
    global nodemap
    global elemnodes
    global elemmap

    answer = []
    # this pattern covers simple slave relations and rigid arm nodes
    match = mastermask_re.search (nodes[inode])
    if (match):
        #print "master detected", inode
        nm = int(match.group(1))
        en = match.group(2).split()
        for i in range(nm):
            inode = int(en[i]) # global number
            answer.append(nodemap[inode]) # local numbering
        return answer
    # match slave nodes
    match = slavenodemasters_re.search (nodes[inode])
    if (match):
        #print "master detected", inode
        nm = int(match.group(1))
        en = match.group(2).split()
        for i in range(nm):
            inode = int(en[i]) # global number
            answer.append(nodemap[inode]) # local numbering
        return answer

    match = hangingnode_re.search (nodes[inode])
    if (match):
        pattern = r'masterElement\s+(\d+)\s+'
        match = re.search(pattern, nodes[inode], re.IGNORECASE)
        if (match):
            e = int(match.group(1))
            # add all elem nodes to master list
            answer.append(elemnodes[elemmap[e]])
            return answer
        else:
            # hanging node without masterElement
            print ("Warning: Support for hanging nodes wwithout masterElement not available")
            exit (1)

    return answer

# returns a list containing for each node follwing tuple:
# (statuses, partitions), where status is local or shared and
# partitions is a list containing node partitions
def classifyNodes (elem_part):
    global nodes, nelem
    global ndofman
    
    nodalstatuses = [0] * ndofman
    nodalpartitions=[]
    # loop over dofmans
    for i in range(ndofman):
        nodalpartitions.append(set())
        # loop over elements sharing node
        for j in nodalconnectivity[i]:
            nodalpartitions[i].add(elem_part[j])

        if (len(nodalpartitions[i]) == 1):
            nodalstatuses[i] = LOCAL
        else:
            nodalstatuses[i]=SHARED

    # account for master->slave conections (here only simple slave dofs handled)
    # loop over all nodes
    for inode in range(ndofman):
        masters = giveNodeMasters(inode) #masters empty if none 
        for m in masters:
            # add inode partition to master partition
            
            if (nodalstatuses[inode] == LOCAL):
                if (nodalpartitions[inode].issubset(nodalpartitions[m])):
                    pass # dependence already satisfied
                else:
                    # make sure that master is accesible from inode partition
                    nodalstatuses[m] = SHARED
                    nodalpartitions[m].update(nodalpartitions[inode])
            else:
                # inode is shared => master has to be shared as well on the same partitions
                nodalstatuses[m] = SHARED
                nodalpartitions[m].update(nodalpartitions[inode])

    return (nodalstatuses, nodalpartitions)



# returns (loc_nodes, shared_nodes) tuple containing
# loc_nodes dict of local nodes on partition i
# shared_nodes dict of shared nodes 
def getPartitionNodeList (i, elem_part):
    global nodes

    loc_nodes=set()
    shared_nodes={}
    for j in range(nelem):
        ep = elem_part[j]
        if (ep == i):
            
            # element on i-th partition found
            # loop over element nodal list
            for inode in elemnodes[j]:
                # now check if all elements shared by node are assigned to the same i-th partition => local node, shared otherwise
                local = True
                for k in nodalconnectivity[inode]:
                    if (elem_part[k] != i): 
                        local = False
                        
                        shared_nodes.setdefault(inode, set()).add(elem_part[k])
                if (local):
                    
                    loc_nodes.add(inode)

    # account for master->slave conections (here only simple slave dofs handled)
    # loop over all nodes
    for inode in nodes:
        masters = giveNodeMasters(inode)
        for m in masters:
            # mark masters as shared to make sure that master is accessible in i-th partition
            if (inode in loc_nodes) and not (m in loc_nodes): #master is on remote partition
                shared_nodes.setdefault(inode, set()).add(elem_part[k])

            if not ((inode in loc_nodes) and (m in loc_nodes)):
                # add master into a list of shared nodes
                shared_nodes.setdefault(m, set()).add(elem_part[k])
            
    return (loc_nodes, shared_nodes)
                    

def writePartition (i, part_vert, nodalstatuses, nodalpartitions):
    global header, bottom
    global ndofman, nodes, elements, nelem
    # open partition input file (for writting)
    pfile = open(inputFileName+'.'+str(i), "w")
    # write outpuf file name first
    pfile.write(header[0].rstrip('\r\n')+'.'+str(i)+'\n')
    # write rest of the header
    for j in range(len(header)-1):
        pfile.write(header[j+1])
    #
    #(ln, sn) = getPartitionNodeList(i, part_vert)
    #
    if debug:
        print ("Partition ", i)
        print ("  Local nodes:",)
        for j in range(ndofman):
            if ((nodalstatuses[j]==LOCAL) and (i in nodalpartitions[j])):
                print (j,)
        print ("  Shared nodes:",) 
        for j in range(ndofman):
            if ((nodalstatuses[j]==SHARED) and (i in nodalpartitions[j])):
                print (j,)
        print ("  Elements:",) 
        for k in range(len(part_vert)):
            if (part_vert[k]==i): print (k,)
        print ("\n")

    nloc = 0
    nshd = 0
    for j in range(ndofman):
        if ((nodalstatuses[j]==LOCAL) and (i in nodalpartitions[j])):
            nloc = nloc+1
        if ((nodalstatuses[j]==SHARED) and (i in nodalpartitions[j])):
            nshd = nshd+1

    # write component record
    pfile.write("ndofman "+str(nloc+nshd) + ' nelem '+str( part_vert.count(i)) + ' ncrosssect '+str(ncs)+' nmat '+str(nmat)+' nbc '+str(nbc)+' nic '+str(nic)+' nltf '+str(nltf)+'\n')
    # write local nodes first
    nloc = 0
    for j in range(ndofman):
        if ((nodalstatuses[j]==LOCAL) and (i in nodalpartitions[j])):
            pfile.write(nodes[j])
            nloc = nloc+1
    # write shared nodes
    nshared = 0
    for j in range(ndofman):
        if ((nodalstatuses[j]==SHARED) and (i in nodalpartitions[j])):
            nshared = nshared+1
            rec = nodes[j].rstrip('\r\n')+' shared partitions '+str(len(nodalpartitions[j]))
            for k in nodalpartitions[j]: rec=rec+' '+str(k)
            pfile.write(rec+'\n')
    # write element records
    for j in range(nelem):
        if (part_vert[j] == i): pfile.write(elements[j])
    # write the bottom part
    for j in bottom:
        pfile.write(j)
    pfile.close()
    
    #print some stats
    print ('{0:>9d} {1:>12d} {2:>12d} {3:>10d}'.format(i, nloc, nshared, part_vert.count(i)))


def parseInput(infile):
    global header, bottom
    global ndofman, nelem, ncs, nmat, nbc, nic, nltf
    global nodes, nodemap, elements, elemmap, elemnodes, nodalconnectivity
    global nodemap
    # loop over input file lines
    # first read first 5 lines (output file, description, emodel, domain, output manager) and store them in header list
    header = []
    for i in range (5):
        header.append(readRecord(infile))
    # check if there are an additional records in header section (export modules)
    match=re.search('\s+nmodules\s+(\d+)\s+', header[2], re.IGNORECASE)
    if match:
        nmodules = int (match.group(1))
    else:
        nmodules = 0

    # read additional nmodules records into header section
    for i in range (nmodules):
        header.append(readRecord(infile))

    # now component record should be there
    componentrec = readRecord(infile)
    ndofman = int (getKeywordValue (componentrec, 'ndofman'))
    nelem = int (getKeywordValue (componentrec, 'nelem'))
    ncs = int (getKeywordValue (componentrec, 'ncrosssect'))
    nmat = int (getKeywordValue (componentrec, 'nmat'))
    nbc = int (getKeywordValue (componentrec, 'nbc'))
    nic = int (getKeywordValue (componentrec, 'nic'))
    nltf = int (getKeywordValue (componentrec, 'nltf'))

    # read nodes
    nodes={}
    nodemap = {}
    for i in range(ndofman):
        rec = readRecord(infile)
        match=re.search('^\w+\s+(\d+)\s+', rec)
        if (match):
            num = int (match.group(1))
        else:
            exit(1);
        nodes[i]=rec
        nodemap[num]=i
    # read elements and create node connectivity array
    elements={}
    elemmap ={}
    elemnodes=[]
    nodalconnectivity=[[] for i in range(ndofman+1)]
    for ie in range(nelem):
        rec = readRecord(infile)
        # match component number
        match=re.search('^\w+\s+(\d+)\s+', rec)
        if (match):
            num = int (match.group(1))
        else:
            exit(1);
        elements[ie]=rec
        elemmap[num]=ie

        # read number of element nodes
        nnode = int (getKeywordValue (rec, 'nodes'))
        enodes=[]
        # match the nodes of element
        pattern = r'\s+nodes\s+\d+\s+((\d+\s+)+)'
        match = re.search (pattern, rec, re.IGNORECASE)
        if (match):
            en = match.group(1).split()
            for i in range(nnode):
                inode = int(en[i])
                nodalconnectivity[nodemap[inode]].append(ie)
                enodes.append(nodemap[inode]) # local numbering
        else:
            print ("Unable to parse element nodal rec: ", rec)
            exit (1);
        elemnodes.append(enodes)



    # store remaining records into bottom buffer (they will be replicated in every subdomain)
    bottom=[]
    for i in range (ncs+nmat+nbc+nic+nltf):
        bottom.append(readRecord(infile))


#
#
# ---------------------------begin here----------------------------
#
#
# parse command line 
options, remainder = getopt.getopt(sys.argv[1:], 'f:n:h', ['input=', 'np=', 'help'])
for opt, arg in options:
    if opt in ('-f', '--input'):
        inputFileName = arg
    elif opt in ('-n', '--np'):
        nparts = int(arg)
    elif opt in ('-h', '--help'):
        print_usage()
        exit(1)

if ((inputFileName == '') or (nparts == 0)):
    print_usage()
    exit(1)
    
# open input file to partition
infile = open(inputFileName, 'r')

# parse input file
parseInput (infile)

#now create adjacency list 
#nodecut assumed - elements graph verices, nodes represent edges
adj = [set() for i in range(nelem)] #emtpy adj list
for ie in range(nelem):
    #print "element", ie
    #loop over element nodes (local numbering)
    for i in elemnodes[ie]:
        #print "  node", i, " connectivity ",  nodalconnectivity[i]
        for k in nodalconnectivity[i]:
            if k != ie:
                adj[ie].add(k)
                adj[k].add(ie)

#print adj

#done; call metis now
try:
    from metis import part_graph
    cuts, part_vert = part_graph(adj, nparts)
except:
    print ("metis module not installed or internal metis error encountered")
    exit(0)

#print "Metis"
#print " part_vert :", part_vert
#print " cuts:", cuts

# write partitioned mesh on output
print ("Partition  local_nodes shared_nodes   elements")
print ("----------------------------------------------")
nodalstatuses, nodalpart = classifyNodes(part_vert)
for i in range (nparts):
    writePartition(i, part_vert, nodalstatuses, nodalpart)

#print nodalstatuses
#print nodalpart


print ("Done")
