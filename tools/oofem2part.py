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
from sets import Set
import getopt

#debug flag
debug = 0



#----do not change------
# number of partitions
nparts = 0
# input File Name 
inputFileName = ''

def print_usage():
    print "\nUsage:\noofem2part -f file -n #"
    print "where -f file sets path to input oofem (serial) file to be partitioned"
    print "      -n #    allows to se the required number of target partitions"


# returns the value corresponding to given keyword and record
def getKeywordValue (record, kwd, optional = None):
    match = re.search (kwd+'\s+\"*([\\\.\+\-,:\w]+)\"*', record)
    if match:
        return match.group(1)
    else:
        #issue an error if argument compulsory
        if optional == None:
            print "\nMissing keyword \"", kwd, "\" in\n", record
            exit (1)
        else:
            return optional


def readRecord (file):
    line = file.readline()
    while (re.search (r"^#", line)):
        line = file.readline()

    #print line
    return line


# returns (loc_nodes, shared_nodes) tuple containing
# loc_nodes dict of local nodes on partition i
# shared_nodes dict of shared nodes 
def getPartitionNodeList (i, elem_part):
    loc_nodes=Set()
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
                        
                        shared_nodes.setdefault(inode, Set()).add(elem_part[k])
                if (local):
                    
                    loc_nodes.add(inode)
    return (loc_nodes, shared_nodes)
                    

def writePartition (i, part_vert):
    global header, bottom
    global nodes, elements, nelem
    # open partition input file (for writting)
    pfile = open(inputFileName+'.'+str(i), "w")
    # write outpuf file name first
    pfile.write(header[0].rstrip('\r\n')+'.'+str(i)+'\n')
    # write rest of the header
    for j in range(len(header)-1):
        pfile.write(header[j+1])
    (ln, sn) = getPartitionNodeList(i, part_vert)

    if debug:
        print "Partition ", i
        print "  Local nodes:", ln
        print "  Shared nodes:", sn.keys()
        print "  Elements:", 
        for k in range(len(part_vert)):
            if (part_vert[k]==i): print k,
        print "\n"


    # write component record
    pfile.write("ndofman "+str( len(ln)+len(sn)) + ' nelem '+str( part_vert.count(i)) + ' ncrosssect '+str(ncs)+' nmat '+str(nmat)+' nbc '+str(nbc)+' nic '+str(nic)+' nltf '+str(nltf)+'\n')
    # write local nodes first
    for j in ln:
        pfile.write(nodes[j])
    # write shared nodes
    for j in sn:
        rec = nodes[j].rstrip('\r\n')+' shared partitions '+str(len(sn[j]))
        for k in sn[j]: rec=rec+' '+str(k)
        pfile.write(rec+'\n')
    # write element records
    for j in range(nelem):
        if (part_vert[j] == i): pfile.write(elements[j])
    # write the bottom part
    for j in bottom:
        pfile.write(j)
    pfile.close()
    
    #print some stats
    print '{0:>9d} {1:>12d} {2:>12d} {3:>10d}'.format(i, len(ln), len(sn), part_vert.count(i))


def parseInput(infile):
    global header, bottom
    global ndofman, nelem, ncs, nmat, nbc, nic, nltf
    global nodes, nodemap, elements, elemmap, elemnodes, nodalconnectivity
    # loop over input file lines
    # first read first 5 lines (output file, description, emodel, domain, output manager) and store them in header list
    header = []
    for i in range (5):
        header.append(readRecord(infile))
    # check if there are an additional records in header section (export modules)
    match=re.search('\s+nmodules\s+(\d+)\s+', header[2])
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
        match = re.search (pattern, rec)
        if (match):
            en = match.group(1).split()
            for i in range(nnode):
                inode = int(en[i])
                nodalconnectivity[nodemap[inode]].append(ie)
                enodes.append(nodemap[inode]) # local numbering
        else:
            print "Unable to parse element nodal rec: ", rec
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
adj = [Set() for i in range(nelem)] #emtpy adj list
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
    print "metis module not installed or internal metis error encountered"

#print "Metis"
#print " part_vert :", part_vert
#print " cuts:", cuts

# write partitioned mesh on output
print "Partition  local_nodes shared_nodes   elements"
print "----------------------------------------------"
for i in range (nparts):
    writePartition(i, part_vert)

print "Done"
