#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#  Extractor.py        (c) 2009 Borek Patzak, www.oofem.org
#
import re
import getopt, sys, os.path

rt_timestep = 9999
rt_dofman   = 9998
rt_dof      = 9997
rt_elem     = 9996
rt_reaction = 9995
rt_loadlevel= 9994

global rectime, recnumber, recdofnum, recvalue
global firstTimeStepFlag, debug
global tolerance

# TODO:
# - add reading of userrec from file
# - add support for extractor/checker mode (or split into separate scripts)
# - add support for additional records -


#record tuple format
#('nr', solution_step, node_id, dof_id, type, value) - dof manager record
#('er', solution_step, elem_id, ir_id, gp_id, 'keyword', keyword_indx, value) - element record
#('ber', solution_step, elem_id, 'keyword', keyword_indx, value) - beam element record
#('rr', solution_step, node_id, dof_id, value) - reaction
#('llr',solution_step, value) - load level record
#('time') - time, only extractor mode
#('include', result) - inclusion and processing of another file

#set default tolerance
tolerance = 1.0e-4

# debug flag, set to 1 for debugging info being printed
debug=0

#recursion level
recursion_level=0

#default mode 'e' - Extractor
#mode 'c' - Checker, can be turned on on command line with '-c' option
mode='e'


class Context:
    def __init__(self):
        self.userrec = []   #list of record tuples
        self.recVal = {}    #parsed values
        self.firstTimeStepFlag = 1 #flag indication that the first time step encountered



begin_re = re.compile(r"""
        ^\#%BEGIN_CHECK%
        """, re.X)

end_re   = re.compile(r"""
        ^           #beginning of the line
        \#%END_CHECK%
        """, re.X)

timeStep_re = re.compile(r"""
        ^           # begining of line
        Output\ for\ time\s*  # characteristic string
        ([-+ ]\d+\.\d+(e[+-]\d+)*) # time step
        """,re.X)

dofMan_re = re.compile(r"""
        ^           # beginning of line
        (?:Node|RigidArmNode|HangingNode)\s*         # char string
        (\d+).*         # node label
        """,re.X)

dof_re    = re.compile(r"""
        ^           # beginning of line
        \s*dof\s*       # char string
        (\d+)           # dof number
        """,re.X)

doftypes_re = re.compile(r"""
        (?:([dvatfp])\s+([-]?\d+\.\d+(?:e[+-]\d+)*))+  #dof types and values on one line
        """, re.X)

element_re = re.compile(r"""
        ^           # begining of line
        element\s*      # element string
        (\d+).*         # element number (label)
        """,re.X)

beamelement_re = re.compile(r"""
        ^           # begining of line
        (beam|spring)\ element\s*    # element string
        (\d+).*         # element number (label)
        """,re.X)


gp_re = re.compile(r"""
        ^           # begining of line
        \s*GP\s*        # gp string
        (\d+)           # irule number
        \.
        (\d+)           # gp number
        \s*:.*          # eat some chars
        """,re.X)

gpstress_re = re.compile (r"""
        stresses\s*
        ([\s+-e\d]+)
        """,re.X)

gpstrain_re = re.compile (r"""
        strains\s*
        ([\s+-e\d]+)
        """,re.X)

gpstrain_re = re.compile (r"""
        strains\s*
        ([\s+-e\d]+)
        """,re.X)

gpstatus_re = re.compile (r"""
        status\s.*
        """, re.X)

gpstate_re = re.compile (r"""
        state\s*
        ([\s+-e\d]+)
        """,re.X)

gpflow_re = re.compile (r"""
        flow\s*
        ([\s+-e\d]+)
        """,re.X)

gpDoH_re = re.compile (r"""
        DoH\s*
        ([\s+-e\d]+)
        """,re.X)

gpHeatPower_re = re.compile (r"""
        HeatPower\s*
        ([\s+-e\d]+)
        """,re.X)

beamrec_re  = re.compile (r"""
        displacements|forces
        """, re.X)

beamdispl_re = re.compile (r"""
        local\ displacements
        """, re.X)

beamforces_re = re.compile (r"""
        local\ end\ forces
        """, re.X)

springforces_re = re.compile (r"""
        spring\ force\ or\ moment
        """, re.X)


reaction_re = re.compile (r"""
        (\w+)\s+(\d+)\s+
        iDof
        \s+(\d+)\s+
        reaction\s+
        ([-]*\d+\.\d+(e[+-]\d+)?) # value
        """,re.X)

loadlevel_re = re.compile (r"""
        load\ level\s+:\s+
        ([-]*\d+\.\d+(e[+-]\d+)?) # value
        """,re.X)

include_re = re.compile (r"""
        ^\#(INCLUDE|include)\s+
        ([\w\.]+)
	\s*
        """, re.X)

# returns the value corresponding to given keyword and record
def getKeywordValue (infilename, record, kwd, optional = None):
    match = re.search (kwd+'\s+\"*([\\\.\+\-,:\w]+)\"*', record)
    if match:
        return match.group(1)
    else:
        #issue an error if argument compulsory
        if optional == None:
            print "\nChecker.py:%55s   "%os.path.basename(infilename)
            print "\nMissing keyword \"", kwd, "\" in\n", record
            exit (1)
        else:
            return optional

# parses the extractor/checker input record
def parse_input_rec (context, recline):
    if re.search('^#(DOFMAN|NODE)', recline):
        if (mode == 'c'): tstep = float(getKeywordValue(context.infilename, recline, 'tStep'))
        else: tstep = 0
        try:
            number= int(getKeywordValue(context.infilename, recline, 'number'))
            dof   = int(getKeywordValue(context.infilename, recline, 'dof'))
            type  = getKeywordValue(context.infilename, recline, '(?:type|unknown)')
            value = float(getKeywordValue(context.infilename, recline, 'value',0.0))
            return ('nr', tstep, number, dof, type, value)

        except ValueError:
            print "Input error on\n",recline
            return None
    elif re.search('^#ELEMENT', recline):
        if (mode == 'c'): tstep = float(getKeywordValue(context.infilename, recline, 'tStep'))
        else: tstep = 0
        try:
            number= int(getKeywordValue(context.infilename, recline, 'number'))
            irule = int(getKeywordValue(context.infilename, recline, 'irule', 1))
            gp    = int(getKeywordValue(context.infilename, recline, 'gp'))
            #gprec = int(getKeywordValue(context.infilename, recline, 'record'))
            kwd   = getKeywordValue(context.infilename, recline, 'keyword')
            cmpn  = int(getKeywordValue(context.infilename, recline, 'component'))
            value = float(getKeywordValue(context.infilename, recline, 'value', 0.0))
            return ('er', tstep, number, irule, gp, kwd, cmpn, value)

        except ValueError:
            print "Input error on\n",recline
            return None

    elif re.search('^#BEAM_ELEMENT', recline):
        if (mode == 'c'): tstep = float(getKeywordValue(context.infilename, recline, 'tStep'))
        else: tstep = 0
        try:
            number= int(getKeywordValue(context.infilename, recline, 'number'))
            #gprec = int(getKeywordValue(context.infilename, recline, 'record'))
            kwd   = getKeywordValue(context.infilename, recline, 'keyword')
            cmpn  = int(getKeywordValue(context.infilename, recline, 'component'))
            value = float(getKeywordValue(context.infilename, recline, 'value', 0.0))
            return ('ber', tstep, number, kwd, cmpn, value)
        except ValueError:
            print "Input error on\n",recline
            return None

    elif re.search('^#REACTION',recline):
        if (mode == 'c'): tstep = float(getKeywordValue(context.infilename, recline, 'tStep'))
        else: tstep = 0
        try:
            number= int(getKeywordValue(context.infilename, recline, 'number'))
            dof   = int(getKeywordValue(context.infilename, recline, 'dof'))
            value = float(getKeywordValue(context.infilename, recline, 'value', 0.0))
            return ('rr', tstep, number, dof, value)
        except ValueError:
            print "Input error on\n",recline
            return None

    elif re.search('^#LOADLEVEL',recline):
        if (mode == 'c'): tstep = float(getKeywordValue(context.infilename, recline, 'tStep'))
        else: tstep = 0
        try:
            value = float(getKeywordValue(context.infilename, recline, 'value', 0.0))
            return ('llr', tstep, value)
        except ValueError:
            print "Input error on\n",recline
            return None
    elif (mode == 'e') and re.search('^#TIME',recline):
        return ('time', 0.0, 0.0)

    else:
        return None


# extract dofman record if apply for actual record
def check_node_rec (context):
    global mode
    for irec, rec in enumerate(context.userrec):
        if (mode == 'e'): timeflag = 1
        else: timeflag = (rec[1] == context.rectime)
        for tv in context.rectypesvalues: # loop over list of dof values and their types
            (type,value) = tv
            if ((rec[0]=='nr') and timeflag and (context.recnumber == rec[2]) and (context.recdofnum == rec[3]) and (type == rec[4])):
                context.recVal[irec]=value;

#extract element record
def check_element_rec (context, kwd, line):
    for irec,rec in enumerate(context.userrec):
        if (mode == 'e'): timeflag = 1
        else: timeflag = (rec[1] == context.rectime)

        if ((rec[0]=='er') and timeflag and (context.recnumber == rec[2]) and
            (context.recirule==rec[3]) and (context.recgpnum == rec[4])):
#               print "Found er: looking for ",rec[5],"in ", line
            match=re.search(rec[5]+'\s*((([-]*\d+(\.\d+)?(e[+-]\d+)?)\s*)+)', line)
            if match:
                # print "match: ",match.group(1)
                context.recVal[irec]=re.split('\s+',match.group(1))[rec[6]-1]
                # print "found: ",recVal[irec]


#extract element record
def check_beam_rec (context, kwd, line):
    for irec,rec in enumerate(context.userrec):
        if (mode == 'e'): timeflag = 1
        else: timeflag = (rec[1] == context.rectime)

        if ((rec[0]=='ber') and timeflag and (context.recnumber == rec[2])):
#            print "Found ber: looking for ",rec[3],"in ", line
            match=re.search(rec[3]+'\s*(([-]*\d+(\.\d+)?(e[+-]\d+)?)\s*)+', line)
            if match:
                context.recVal[irec]=re.split('\s+',match.group(0))[rec[4]]
#                print "found\n"

#extract reaction record
def check_reaction_rec (context):
    for irec,rec in enumerate(context.userrec):
        if (mode == 'e'): timeflag = 1
        else: timeflag = (rec[1] == context.rectime)

        if ((rec[0] == 'rr') and timeflag and (rec[2]==context.recnumber) and (rec[3]==context.recdofnum)):
            context.recVal[irec]=context.recvalue

#extract load level record
def check_loadlevel_rec (context):
    for irec,rec in enumerate(context.userrec):

        if (mode == 'e'): timeflag = 1
        else: timeflag = (rec[1] == context.rectime)

        if ((rec[0] == 'llr') and timeflag):
            context.recVal[irec]=context.recvalue

#check time rec
def check_time_rec (context):
    for irec,rec in enumerate(context.userrec):

        if (mode == 'e') and (rec[0] == 'time'):
            context.recVal[irec]=context.rectime



def match_primary_rec (context, line):
    global mode, debug

    match=timeStep_re.search(line)
    if match:
        context.rectype = rt_timestep
        context.rectime=float(match.group(1))
        if debug: print "found time ", context.rectime
        if context.firstTimeStepFlag:
            context.firstTimeStepFlag = 0
        # print parsed record values from previous step
        elif mode == 'e': print_step_results(context)
        check_time_rec (context)
        return None
    match=loadlevel_re.search(line)
    if match:
        context.rectype = rt_loadlevel
        context.recvalue= float(match.group(1))
        if debug: print "found load level ", context.recvalue
        check_loadlevel_rec (context)
        return None


    match=dofMan_re.search(line)
    if match:
        context.rectype = rt_dofman
        context.recnumber= int(match.group(1))
        if debug: print "found  node", context.recnumber
        nline = match_dofrec(context)
        return nline

    match=element_re.search(line)
    if match:
        context.rectype = rt_elem
        context.recnumber = int(match.group(1))
        if debug: print "found element", context.recnumber
        nline = match_gprec(context)
        return nline

    match=beamelement_re.search(line)
    if match:
        context.rectype = rt_elem
        context.recnumber = int(match.group(2))
        if debug: print "found element", context.recnumber
        nline = match_beamrec(context)
        return nline

    match=reaction_re.search(line)
    if match:
        if debug: print "Rea: ", line
        context.rectype = rt_reaction
        context.recnumber = int(match.group(2))
        context.recdofnum = int(match.group(3))
        context.recvalue  = float(match.group(4))
        check_reaction_rec (context)
        return None

def match_dofrec (context):
    global debug

    # parse dof records
    for line in context.infile:
        match=dof_re.search(line)
        if match:
            context.rectype = rt_dof
            context.recdofnum= int(match.group(1))
            match=doftypes_re.findall(line)
            context.rectypesvalues = match
            #rectype  = match.group(2)
            #recvalue = float(match.group(3))
            check_node_rec (context);
            if debug: print "     dof", context.recdofnum
            continue
        else:
            return line


def match_gpsubrec (context, aline):
    global debug
    pmatch=gpstress_re.search(aline)
    if pmatch:
        check_element_rec (context, 'stresses', aline)
        if debug: print "     stress rec"
        return 1
    ppmatch = gpstrain_re.search(aline)
    if ppmatch:
        check_element_rec (context, 'strains', aline)
        if debug: print "     strain rec"
        return 1
    ppmatch = gpstatus_re.search(aline)
    if ppmatch:
        check_element_rec (context, 'status', aline)
        if debug: print "     status rec"
        return 1
    #state variables in transport problems
    ppmatch = gpstate_re.search(aline)
    if ppmatch:
        check_element_rec (context, 'state', aline)
        if debug: print "     state rec"
        return 1
    #flow vector in transport problems
    ppmatch = gpflow_re.search(aline)
    if ppmatch:
        check_element_rec (context, 'flow', aline)
        if debug: print "     flow rec"
        return 1
    #Degree of hydration in cement hydration models
    ppmatch = gpDoH_re.search(aline)
    if ppmatch:
        check_element_rec (context, 'DoH', aline)
        if debug: print "     DoH rec"
        return 1
    ppmatch = gpHeatPower_re.search(aline)
    if ppmatch:
        check_element_rec (context, 'HeatPower', aline)
        if debug: print "     HeatPower rec"
        return 1
    return 0

def match_singlegprec (context, line):
    global rectime, recnumber, recirule, recgpnum, debug
    match=gp_re.search(line)
    if match:
        context.recirule=int(match.group(1))
        context.recgpnum=int(match.group(2))
        if debug: print "  gp", context.recgpnum
        match_gpsubrec (context, line)
        for line in context.infile:
            match=gp_re.search(line)
            if match:
                return line;
            result = match_gpsubrec (context, line)
            if not result == 1 :
                return line
    return line

def match_gprec (context):
    # parse possible gp records
    for line in context.infile:
        match=gp_re.search(line)
        while (match):
            line = match_singlegprec (context, line)
            match=gp_re.search(line)

        else:
            return line

def match_beamrec (context):
    global rectime, recnumber, recgpnum, debug
    for line in context.infile:
        pmatch=beamdispl_re.search(line)
        if pmatch:
            check_beam_rec (context, 'displacements', line)
            if debug: print "     local displacements rec"
            continue
        ppmatch = beamforces_re.search(line)
        if ppmatch:
            check_beam_rec (context, 'forces', line)
            if debug: print "     local forces rec"
            continue
        ppmatch = springforces_re.search(line)
        if ppmatch:
            check_beam_rec (context, 'huhu', line)
            if debug: print "     spring rec"
            continue

        return line


def check_results (context, tolerance):
    success = 1
    if ( context.parentfilename):
        header= "Checker.py: "+(os.path.basename(context.parentfilename)+'->'+os.path.basename(context.infilename)).ljust(55)[:55]
    else :
        header= "Checker.py: "+(os.path.basename(context.infilename)).ljust(55)[:55]

    # loop over recVal items
    for irec, rec in enumerate(context.userrec):
        if (rec[0] == 'include'):
            if (rec[1] > 0): success = 0
            continue
        elif (context.recVal[irec] == '--'):
            if (success): print header # print header before reporting error for the first time
            print "\tError when checking rule ",irec,": not found"
            success = 0
            continue

        try:
            err = (float(rec[-1])-float(context.recVal[irec]))
        except ValueError:
            err = float(rec[-1])
        if (abs(err) > float(tolerance)):
            if (success): print header # print header before reporting error for the first time
            print "\tError when checking rule ",irec,": err = ",err, ", value is ",float(rec[-1])," and should be ",float(context.recVal[irec])
            # print rec, recVal
            success = 0


    if success:
        print header+"[  OK  ]"
        return 0
    else:
        print header+"[FAILED]"
        return 1


def print_step_results (context):
    # loop over recVal items
    for ival in sorted(context.recVal.keys()):
        try:
            value = float(context.recVal[ival])
            print "%15e"%value,
        except ValueError:
            print "%15s"%'---',
        context.recVal[ival]='--'
    print


def usage():
# prints usage
    print """
Checker (c) 2009 Borek Patzak

Usage: extractor.py -f input_file.in [-c]

Where: input_file.in is extractor input file
       -c turns on the checker mode

The syntax of extractor input file is following:
------BEGIN---------
oofem_output_file_name
....
#%BEGIN_CHECK% [tolerance #]
#DOFMAN    {tStep #} number # dof # type [dtfva] {value #}
#ELEMENT   {tStep #} number # [irule #] gp # keyword # component # {value #}
#REACTION  {tStep #} number # dof # {value #}
#LOADLEVEL {tStep #} {value #}
#INCLUDE slave_input_file.in
<#TIME>
#%END_CHECK%
------END-----------
Other lines not matching the syntax are ignored.
The records in {} are required only when in checker mode,
records in <> are available only in extractor mode.
the records in [] are optional.

The type value is a single character determining the type of dof
value, where 'd' stands for displavement, 'v' for velocity,
'a' for acceleration, 't' for temperature, and 'f' for flux.
The output is printed to stdout, one row per each solution step,
In particular columns, the extracted values are printed, preserving
their order in input file.
The include directive allows to process slave extractor input files
from a single master file (one level recursion allowed). If any
slave file fails, master fails as well.

Example (for checker mode):
patch100.out
#%BEGIN_CHECK% tolerance 1.e-4
#ELEMENT tStep 1 number 1 gp 1 component 1 keyword "stresses" value -8.3333e+00
#ELEMENT tStep 1 number 1 gp 1 component 1 keyword "strains"  value -5.2083e-01
#DOFMAN  tStep 1 number 4 dof 1 type d value -1.5625
#%END_CHECK%
"""

def process_file (infilename, parentfilename):
    #parentfilename name of master file name, none otherwise
    global userrec, mode, tolerance, recursion_level
    context = Context();

    result = 0


    # open oofem input file
    if debug: print "Opening oofem input file:", sys.argv[1]
    context.infile=open(infilename)
    context.infilename=infilename
    context.parentfilename=parentfilename

    #read output file record
    for line in context.infile:
        if line[0] != '#':
            oofemoutfilename = line
            break

    #parse rest of oofem input file to find extractor records
    if 1:
        begin=0
        for line in context.infile:
            match = begin_re.search(line)
            if match:
                begin=1
                match=re.search("tolerance\s+(\d+\.\d*(e[+-]\d+)*)", line)
                if match:
                    tolerance = match.group(1)
                    if debug: print "Tolerance = ", tolerance
                break

    if begin==0:
        print "No extrator records found"
	return 2

    end=0
    userrec = []
    for line in context.infile:
        match = include_re.search(line)
        if match:
            # include directive detected, parse given file recursively
            #print "Processing include ", match.group(2)
            recursion_level = recursion_level+1
            if (recursion_level < 2):

                result=result+process_file (match.group(2), context.infilename)
                context.userrec.append (('include',result)) #remember result
		recursion_level = recursion_level-1
            else:
                print "Allowed recursion level reached, file:",match.group(2)
		return 2

        match = parse_input_rec (context, line)
        if match:
            context.userrec.append (match)
        match = end_re.search(line)

        if match:
            end = 1
            break

    if not end:
        print "No end record found"
	return 2

    if debug: print userrec
    context.infile.close()

    #process oofem output file
    if debug: print "Opening oofem output file:", oofemoutfilename.rstrip('\r\n')
    context.infile = open(oofemoutfilename.rstrip('\r\n'))

    for i in range(len(context.userrec)):
        context.recVal[i]='--'
    #print context.recVal

    #parse output file
    for line in context.infile:
        while (line):
            line = match_primary_rec (context, line)


    #print records

    #print recVal
    if mode=='c':
        return result+check_results(context, tolerance)
    elif mode =='e':
        print_step_results(context)
	return 0



#########################################################
# main function
#########################################################
def main():
    global userrec, recval, mode, tolerance
    infilename = ""

    try:
        # options expecting value followed by ':', see also Unix getopt
        opts, args = getopt.getopt (sys.argv[1:],'cf:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        exit(1)

    # flag whether -f compulsory option provided
    infileflag = False
    #print opts
    for o, a in opts:
        if o == "-f":
            filename = a
            infileflag = True
        elif o == "-c":
            mode = 'c'
        else:
            usage()
            assert False, "unhandled option"

    if not infileflag:
        usage()
        assert False, "-f option required"

    return process_file (filename, None)



if __name__ == "__main__":
    sys.exit(main())
