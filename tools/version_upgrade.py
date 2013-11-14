#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This script upgrades the input files to the new format 
#
import re
import fileinput
import getopt, sys, os.path

# New naming convention: "problemtype_approximation"
keyword_renaming = {
"tr21stokes" : "stokes_tr21",
"tr1darcy" : "darcy_tr1",
"l4axisymm" : "structaxisym_quad1",
"quad_10_2d_supg" : "supg_quad10",
}
# ADD MORE!


# Converts from the domain type to the "nsd" parameter
nsd_domaintype = {
"2dplanestressrot" : 2,
"2dplanestress" : 2,
"planestrain" : 2,
"3d" : 3,
"3daxisymm" : 2,
"2dmindlinplate" : 2,
"3dshell" : 3,
"2dtruss" : 2,
"1dtruss" : 1,
"2dbeam" : 2,
"2dlattice" : 2,
"heattransfer" : 3,
"mass1transfer" : 3,
"hema1" : 3,
"2dincompflow" : 2,
"3dincompflow" : 3,
"3ddirshell" : 3,
"2dmasslatticetransport" : 2,
}

D_u = "1"
D_v = "2"
D_w = "3"
R_u = "4"
R_v = "5"
R_w = "6"
V_u = "7"
V_v = "8"
V_w = "9"
T_f = "10"
P_f = "11"
G_0 = "12"
G_1 = "13"
C_1 = "14"
W_u = "15"
W_v = "16"
W_w = "17"
Gamma = "18"
D_u_edge_const = "19"
D_u_edge_lin = "20"
D_v_edge_const = "21"
D_v_edge_lin = "22"

dofs_domaintype = {
"2dplanestressrot" : [D_u, D_v, R_w],
"2dplanestress" : [D_u, D_v],
"planestrain" : [D_u, D_v],
"3d" : [D_u, D_v, D_w],
"3daxisymm" : [D_u, D_v, R_w],
"2dmindlinplate" : [D_w, R_u, R_v],
"3dshell" : [D_u, D_v, D_w, R_u, R_v, R_w],
"2dtruss" : [D_u, D_v],
"1dtruss" : [D_u],
"2dbeam" : [D_u, D_w, R_v],
"2dlattice" : [D_u, D_v, R_w],
"heattransfer" : [T_f],
"mass1transfer" : [C_1],
"hema1" : [T_f, C_1],
"2dincompflow" : [V_u, V_v, P_f],
"3dincompflow" : [V_u, V_v, V_w, P_f],
"3ddirshell" : [D_u, D_v, D_w, W_u, W_v, W_w, Gamma],
"2dmasslatticetransport" : [P_f],
}

# 
def convertToLower(line):
    activeString = False
    for i in range(line):
        if line[i] == '"':
            activeString = not activeString
        if activeString:
            out[i] = line[i].lower()
        else:
            out[i] = line[i]
    return out


def necessarySets(bcs):
    


def main():
    inputfile = "test.in"
    outputfile = "mod_test.in"
    print("Opening oofem input file:", inputfile)

    fileinput = open(inputfile, 'r')
    fileoutput = open(outputile, 'w')

    reallinenumber = -1
    for line in fileinput.input():
        process(line)
        if line[0] = "#":
            outline = line
        else:
            reallinenumber += 1
            if reallinenumber < 2: # We're past filename and description now
                outline = line
            else:
                line = convertToLower(line)
                words = line.split()

                # RENAME RECORD KEYWORDS:
                if words[0] in keyword_renaming:
                    words[0] = keyword_renaming[words[0]]

                # FETCH THE INFO FROM THE DEPRECATED DOMAIN RECORD
                if words[0] == "domain":
                    nsd = nsd_domaintype[words[1]]
                    defaultdofs = dofs_domaintype[words[1]]
                    words[0] = "#domain"
                    axisymm = words[1] == "3daxisymm"

                # MODIFY THE DOMAIN RECORD LINE
                if words[0] == "ndofman":
                    words.insert(0, "domain")
                    words.extend(["nsd", "%d"%nsd])
                    if not "nset" in words:
                        words.extend(["nset", "0"])

                # Modify the node lines with the default dofidmasks and such
                if "coords" in words:
                    hasDofIdMask = "dofidmask" in words
                    hasMasterMask = "mastermask" in words
                    hasDofType = "doftype" in words
                    if not hasDofIdMask and ( hasMasterMask or hasDofType ):
                        words.extend(["dofidmask", len(defaultdofs)])
                        words.extend(defaultdofs)

                    # This is very messy to convert to the new format:
                    # Copy over the sets to this vector
                    if "bc" in hasBC:
                        index = words.index("bc")
                        size = int(words[index+1])
                        dofs = defaultdofs
                        if "dofidmask" in words:
                            index2 = words.index("dofidmask")
                            size2 = int(words[index2+1])
                            dofs = words[index2+2:index2+3+size2]
                        old_bc[words[1]] = zip((dofs, words[index+2:index+3+size])
                        words = words[:index] + words[index+3+size:]

                    if "load" in words:
                        index = words.index("load")
                        size = int(words[index+1])
                        dofs = defaultdofs
                        if "dofidmask" in words:
                            index2 = words.index("dofidmask")
                            size2 = int(words[index2+1])
                            dofs = words[index2+2:index2+3+size2]
                        old_nodaloads[words[1]] = zip((dofs, words[index+2:index+3+size])
                        words = words[:index] + words[index+3+size:]

                # Modify the element line
                if "nodes" in words:
                    if "mat" in words:
                        index = words.index("mat")
                        cs_index = words.index("crosssect")
                        cs_mat[words[cs_index+1]] = words[index+1]
                        words = words[:index] + words[index+2:]
                    if "bodyloads" in words:
                        index = words.index("bodyloads")
                        size = int(words[index+1])
                        vec = words[index+2:index+3+size]
                        old_bodyloads[words[1]] = zip(vec[:2:], vec[1:2:])
                        words = words[:index] + words[index+3+size:]
                    if "boundaryloads" in words:
                        index = words.index("boundaryloads")
                        size = int(words[index+1])
                        vec = words[index+2:index+3+size]
                        old_boundaryloads[words[1]] = zip(vec[:2:], vec[1:2:])
                        words = words[:index] + words[index+3+size:]
                    if "edgeloads" in words:
                        index = words.index("edgeloads")
                        size = int(words[index+1])
                        vec = words[index+2:index+3+size]
                        old_edgeloads[words[1]] = zip(vec[:2:], vec[1:2:])
                        words = words[:index] + words[index+3+size:]

                # Add the material to the cross-section if it doesn't already exist.
                if words[0] == "simplecs" or words[0] == "fluidcs" or words[0] == "transportcs":
                    if not "mat" in words:
                        words.extend(["mat",cs_mat[words[1]])
                    if not "mat" in words:
                        words.extend(["mat",cs_mat[words[1]])

                outline = ' '.join(words)
            fileoutput.write(outline + '\n')

if __name__ == "__main__":
    sys.exit(main())
