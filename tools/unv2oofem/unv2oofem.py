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

        #count elements
        nelems=0
        for elem in FEM.elems:
            if elem.oofem_elemtype:
                nelems+=1

        #write component record
        of.write('ndofman %d nelem %d ncrosssect %d nmat %d nbc %d nic %d nltf %d\n' % (FEM.nnodes, nelems, CTRL.ncrosssect, CTRL.nmat, CTRL.nbc, CTRL.nic, CTRL.nltf))
        #write nodes
        for node in FEM.nodes:
            #resolve nodal properties
            properties=""
            for igroup in node.oofem_groups:
                properties+=igroup.oofem_properties
            # write nodal record
            of.write(('node %5d coords 3 %8g %8g %8g %s\n') % (node.id,node.coords[0],node.coords[1],node.coords[2], properties))


        #write elements
        for elem in FEM.elems:
            #resolve element properties
            properties=""
            for igroup in elem.oofem_groups:
                properties+=igroup.oofem_properties
            #do output if elemtype resolved
            if elem.oofem_elemtype:
                dat=[CTRL.elementNames[elem.oofem_elemtype], elem.id,elem.nnodes]
                dat.extend([x for x in elem.cntvt])
                dat.append(properties)
                # format='%5d, ' * (4 + elem.nnodes) + ' %5d'            
                of.write('%s %5d nodes %d %5d %5d %5d %s\n' % tuple(dat))

        # write final sections
        of.write(CTRL.footer);
        of.close()
        #
        t2 = time.time()
        #
        print "done ( %d nodes %d elements)" % (FEM.nnodes, nelems)
        print "Finished in %0.2f [s]" % ((t2-t1))
        
    else:
        print(helpmsg)

        
