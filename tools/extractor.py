#!/usr/bin/python

#
#  Extractor.py            (c) 2009 Borek Patzak, www.oofem.org
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
#('er', solution_step, elem_id, gp_id, 'keyword', keyword_indx, value) - element record
#('ber', solution_step, elem_id, 'keyword', keyword_indx, value) - beam element record
#('rr', solution_step, node_id, dof_id, value) - reaction 
#('llr',solution_step, value) - load level record

#set default tolerance
tolerance = 1.0e-4

# debug flag, set to 1 for debugging info beeing printed
debug=0

#default mode 'e' - Extractor
#mode 'c' - Checker, can be turned on on command line with '-c' option
mode='e'

#list of record tuples
userrec = []
#parsed values 
recVal = {}

#flag indication that the first time step encountered
firstTimeStepFlag = 1


begin_re = re.compile(r"""
                ^\#%BEGIN_CHECK%
                """, re.X)

end_re   = re.compile(r"""
                ^                   #beginning of the line
                \#%END_CHECK%
                """, re.X)

timeStep_re = re.compile(r"""
		^                   # begining of line
		Output\ for\ time\s*  # charcteristic string
		(\d+\.\d+(e[+-]\d+)*) # time step
		""",re.X)

dofMan_re = re.compile(r"""
                ^                   # beginning of line
                (?:Node|RigidArmNode)\s*             # char string
                (\d+).*             # node label
                """,re.X)

dof_re    = re.compile(r"""
                ^                   # beginning of line
                \s*dof\s*           # char string
                (\d+)               # dof number
                \s*([dtfpva])\s*    # dof_type [dtf]
                ([-]*\d+\.\d+(e[+-]\d+)*) # dof value
                """,re.X)                

element_re = re.compile(r"""
                ^                   # begining of line
                element\s*          # element string
                (\d+).*             # element number (label)
                """,re.X)

beamelement_re = re.compile(r"""
                ^                   # begining of line
                beam\ element\s*          # element string
                (\d+).*             # element number (label)
                """,re.X)


gp_re = re.compile(r"""
                ^                   # begining of line
                \s*GP\s*            # gp string
                (\d+)               # gp number
                \s*:.*              # eat some chars
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

beamrec_re  = re.compile (r"""
                displacements|forces
                """, re.X)

beamdispl_re = re.compile (r"""
                local\ displacements
                """, re.X)

beamforces_re = re.compile (r"""
                local\ end\ forces
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


# returns the value corresponding to given keyword and record
def getKeywordValue (record, kwd, optional = None):
	global infilename
	match = re.search (kwd+'\s+\"*([\.\+\-,:\w]+)\"*', record)
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
def parse_input_rec (recline):
	if re.search('^#(DOFMAN|NODE)', recline):
		if (mode == 'c'): tstep = float(getKeywordValue(recline, 'tStep')) 
		else: tstep = 0
		try:
			number= int(getKeywordValue(recline, 'number'))
			dof   = int(getKeywordValue(recline, 'dof'))
			type  = getKeywordValue(recline, '(?:type|unknown)')
			value = float(getKeywordValue(recline, 'value',0.0))
			return ('nr', tstep, number, dof, type, value)

		except ValueError:
			print "Input error on\n",recline
			return None
	elif re.search('^#ELEMENT', recline):
		if (mode == 'c'): tstep = float(getKeywordValue(recline, 'tStep')) 
		else: tstep = 0
		try:
			number= int(getKeywordValue(recline, 'number'))
			gp    = int(getKeywordValue(recline, 'gp'))
     		        #gprec = int(getKeywordValue(recline, 'record'))
			kwd   = getKeywordValue(recline, 'keyword')
			cmpn  = int(getKeywordValue(recline, 'component'))
			value = float(getKeywordValue(recline, 'value', 0.0))
			return ('er', tstep, number, gp, kwd, cmpn, value)
		except ValueError:
			print "Input error on\n",recline
			return None

	elif re.search('^#BEAM_ELEMENT', recline):
		if (mode == 'c'): tstep = float(getKeywordValue(recline, 'tStep')) 
		else: tstep = 0
		try:
			number= int(getKeywordValue(recline, 'number'))
     		        #gprec = int(getKeywordValue(recline, 'record'))
			kwd   = getKeywordValue(recline, 'keyword')
			cmpn  = int(getKeywordValue(recline, 'component'))
			value = float(getKeywordValue(recline, 'value', 0.0))
			return ('ber', tstep, number, kwd, cmpn, value)
		except ValueError:
			print "Input error on\n",recline
			return None
			
	elif re.search('^#REACTION',recline):
		if (mode == 'c'): tstep = float(getKeywordValue(recline, 'tStep')) 
		else: tstep = 0
		try:
			number= int(getKeywordValue(recline, 'number'))
			dof   = int(getKeywordValue(recline, 'dof'))
			value = float(getKeywordValue(recline, 'value', 0.0))
			return ('rr', tstep, number, dof, value)
		except ValueError:
			print "Input error on\n",recline
			return None
		
	elif re.search('^#LOADLEVEL',recline):
		if (mode == 'c'): tstep = float(getKeywordValue(recline, 'tStep')) 
		else: tstep = 0
		try:
			value = float(getKeywordValue(recline, 'value', 0.0))
			return ('llr', tstep, value)
		except ValueError:
			print "Input error on\n",recline
			return None

	else:
		return None

		
# extract dofman record if apply for actual record
def check_node_rec (time, node, dof, type, value):
	global recVal, mode
	for irec, rec in enumerate(userrec):
		if (mode == 'e'): timeflag = 1
		else: timeflag = (rec[1] == time)
		if ((rec[0]=='nr') and timeflag and (node == rec[2]) and (dof == rec[3]) and (type == rec[4])):
			recVal[irec]=value;

#extract element record
def check_element_rec (time, elem, gp, kwd, value, line):
	global recVal
	for irec,rec in enumerate(userrec):
		if (mode == 'e'): timeflag = 1
		else: timeflag = (rec[1] == time)
		
		if ((rec[0]=='er') and timeflag and (elem == rec[2]) and (gp == rec[3])):
#			print "Found er: looking for ",rec[5],"in ", line
			match=re.search(rec[4]+'\s*(([-]*\d+(\.\d+)?(e[+-]\d+)?)\s*)+', line)
			if match:
				recVal[irec]=re.split('\s+',match.group(0))[rec[5]]


#extract element record
def check_beam_rec (time, elem, kwd, value, line):
	global recVal
	for irec,rec in enumerate(userrec):
		if (mode == 'e'): timeflag = 1
		else: timeflag = (rec[1] == time)
		
		if ((rec[0]=='ber') and timeflag and (elem == rec[2])):
#			print "Found ber: looking for ",rec[3],"in ", line
			match=re.search(rec[3]+'\s*(([-]*\d+(\.\d+)?(e[+-]\d+)?)\s*)+', line)
			if match:
				recVal[irec]=re.split('\s+',match.group(0))[rec[4]]
#				print "found\n"
	
#extract reaction record
def check_reaction_rec (time, dofman, idof, value):
	global recVal
	for irec,rec in enumerate(userrec):
		if (mode == 'e'): timeflag = 1
		else: timeflag = (rec[1] == time)

		if ((rec[0] == 'rr') and timeflag and (rec[2]==dofman) and (rec[3]==idof)):
			recVal[irec]=value

#extract load level record
def check_loadlevel_rec (time, value):
	global recVal
	for irec,rec in enumerate(userrec):

		if (mode == 'e'): timeflag = 1
		else: timeflag = (rec[1] == time)

		if ((rec[0] == 'llr') and timeflag):
			recVal[irec]=value

def match_primary_rec (line):
	global rectime, recnumber, recdofnum, recvalue, mode, firstTimeStepFlag, debug

	match=timeStep_re.search(line)
	if match:
		rectype = rt_timestep
		rectime=float(match.group(1))
		if debug: print "found time ",rectime
		if firstTimeStepFlag:
			firstTimeStepFlag = 0
		# print parsed record values from previous step
		elif mode == 'e': print_step_results()
		return
	match=loadlevel_re.search(line)
	if match:
		rectype = rt_loadlevel
		recvalue= float(match.group(1))
		if debug: print "found load level ",recvalue
		check_loadlevel_rec (rectime, recvalue)
		return
	

	match=dofMan_re.search(line)
	if match:
		rectype = rt_dofman
		recnumber= int(match.group(1))
		if debug: print "found  node",recnumber
		nline = match_dofrec()
		match_primary_rec (nline)

	match=element_re.search(line)
	if match:
		rectype = rt_elem
		recnumber = int(match.group(1))
		if debug: print "found element", recnumber
		nline = match_gprec()
		match_primary_rec (nline)

	match=beamelement_re.search(line)
	if match:
		rectype = rt_elem
		recnumber = int(match.group(1))
		if debug: print "found element", recnumber
		nline = match_beamrec()
		match_primary_rec (nline)

	match=reaction_re.search(line)
	if match:
		if debug: print "Rea: ", line
		rectype = rt_reaction
		recnumber = int(match.group(2))
		recdofnum = int(match.group(3))
		recvalue  = float(match.group(4))
		check_reaction_rec (rectime, recnumber, recdofnum, recvalue)


def match_dofrec ():
	global rectime, recnumber, recdofnum, recvalue, debug

	# parse dof records
	for line in infile:
		match=dof_re.search(line)
		if match:
			rectype = rt_dof
			recdofnum= int(match.group(1))
			rectype  = match.group(2)
			recvalue = float(match.group(3))
			check_node_rec (rectime, recnumber, recdofnum, rectype, recvalue);
			if debug: print "         dof", recdofnum
			continue
		else:
			return line


def match_gpsubrec (aline):
	global rectime, recnumber, recgpnum, debug
	pmatch=gpstress_re.search(aline)
	if pmatch:
		check_element_rec (rectime, recnumber, recgpnum, 'stresses', 0.0, aline)
		if debug: print "     stress rec"
		return 1
	ppmatch = gpstrain_re.search(aline)
	if ppmatch:
		check_element_rec (rectime, recnumber, recgpnum, 'strains', 0.0, aline)
		if debug: print "     strain rec"
		return 1
	ppmatch = gpstatus_re.search(aline)
	if ppmatch:
		check_element_rec (rectime, recnumber, recgpnum, 'status', 0.0, aline)
		if debug: print "     status rec"
		return 1
	return 0

def match_singlegprec (line):
	global rectime, recnumber, recgpnum, debug
	match=gp_re.search(line)
	if match:
		recgpnum=int(match.group(1))
		if debug: print "  gp", recgpnum
		match_gpsubrec (line)
		for line in infile:
			match=gp_re.search(line)
			if match:
				return line;
			result = match_gpsubrec (line)
			if not result == 1 :
				return line
	return line

def match_gprec ():
	# parse possible gp records
	for line in infile:
		match=gp_re.search(line)
		while (match):
			line = match_singlegprec (line)
			match=gp_re.search(line)
                        
		else:
			return line

def match_beamrec ():
	global rectime, recnumber, recgpnum, debug
	for line in infile:
		pmatch=beamdispl_re.search(line)
		if pmatch:
			check_beam_rec (rectime, recnumber, 'displacements', 0.0, line)
			if debug: print "     local displacements rec"
			continue
		ppmatch = beamforces_re.search(line)
		if ppmatch:
			check_beam_rec (rectime, recnumber, 'forces', 0.0, line)
			if debug: print "     local forces rec"
			continue
			
		return line


def check_results ():
	global infilename, tolerance
	success = 1
	# loop over recVal items
	for irec, rec in enumerate(userrec):
		if (recVal[irec] == '--'):
			print "\tError when checking rule ",irec,": not found"
			success = 0
			continue
		try:
			err = (float(rec[-1])-float(recVal[irec]))
		except ValueError:
			err = float(rec[-1])
		if (abs(err) > tolerance):
			print "\tError when checking rule ",irec,": err = ",err
			success = 0
	print "Checker.py:%55s   "%os.path.basename(infilename),
	if success:
		print "[  OK  ]"
	else:
		print "[FAILED]"


def print_step_results ():
	# loop over recVal items
	for ival in sorted(recVal.keys()):
		try:
			value = float(recVal[ival])
			print "%15e"%value,
		except ValueError:
			print "%15s"%'---',
		recVal[ival]='--'
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
#ELEMENT   {tStep #} number # gp # keyword # component # {value #}
#REACTION  {tStep #} number # dof # {value #}
#LOADLEVEL {tStep #} {value #}
#%END_CHECK%
------END-----------
Other lines not matching the syntax are ignored.
The records in {} are required only when in checker mode.
The type value is a single character determining the type of dof 
value, where 'd' stands for displavement, 'v' for velocity,
'a' for acceleration, 't' for temperature, and 'f' for flux.
The output is printed to stdout, one row per each solution step,
In particular columns, the extracted values are printed, preserving
their order in input file

Example (for checker mode):
patch100.out
#%BEGIN_CHECK% tolerance 1.e-4
#ELEMENT tStep 1 number 1 gp 1 component 1 keyword "stresses" value -8.3333e+00
#ELEMENT tStep 1 number 1 gp 1 component 1 keyword "strains"  value -5.2083e-01
#DOFMAN  tStep 1 number 4 dof 1 type d value -1.5625
#%END_CHECK%
"""


#########################################################
# main function
#########################################################
def main():
	global userrec, recval, infile, infilename, mode, tolerance
	try:
		# options expecting value followed by ':', see also Unix getopt 
		opts, args = getopt.getopt (sys.argv[1:],'cf:')
	except getopt.GetoptError, err:
		print str(err)
		usage()
		exit(1)

	# flag whether -f compulsory option provided
	infileflag = 0
	#print opts
	for o, a in opts:
		if o == "-f":
			infilename = a
			infileflag = 1
		elif o == "-c":
			mode = 'c'
		else:
			usage()
			assert False, "unhandled option"

	if infileflag == 0:
		usage()
		assert False, "-f option required"


	# open oofem input file
	if debug: print "Opening oofem input file:", sys.argv[1]
	infile=open(infilename)

#read output file record 
	for line in infile:
		if line[0] != '#':
			oofemoutfilename = line
			break

#parse rest of oofem input file to find extractor records
	if 1:
		begin=0
		for line in infile:
			match = begin_re.search(line)
			if match:
				begin=1
				match=re.search("tolerance\s+(\d+\.\d*(e[+-]\d+)*)", line)
				if match:
					tolerance = match.group(1)
					if debug: print "Tolerance = ", tolerance
				break

	if begin==0:
		print "No beging record found"
		exit (1)
				
	end=0
	userrec = []
	for line in infile:
		match = parse_input_rec (line)
		if match:
			userrec.append (match)
		match = end_re.search(line)
		if match:
			end = 1
			break
	if not end:
		print "No end record found"
		exit(1)

	if debug: print userrec
	infile.close()

#process oofem output file
	if debug: print "Opening oofem output file:", oofemoutfilename[:-1]
	infile = open(oofemoutfilename[:-1])
		
	for i in range(len(userrec)):
		recVal[i]='--'
              
	#parse output file 
	for line in infile:
		match_primary_rec (line)


#print records

        #print recVal
	if mode=='c':
		check_results()
	elif mode =='e':
		print_step_results()



if __name__ == "__main__":
	main()
