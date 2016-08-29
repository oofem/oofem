#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-

# Description: Creates a regular cubic grid of 27 node elements.
# Serves as an example of how to generate an input file for OOFEM from Python.
# Author: Mikael Öhman
# License: CC0

from __future__ import print_function, division
import os
import rveToolbox
import numpy as np


def printMesh(folder, name, rveSize, nelem):
    # elname, w = 'Hex21Stokes', 3
    # elname, w = 'Q27Space', 3
    # elname, w = 'LSpace', 2
    # elname, w = 'QBrick1HT', 3
    elname, w = 'Brick1HT', 2
    n = nelem*(w-1) + 1
    X, Y, Z = np.mgrid[-rveSize*0.5:rveSize*0.5:n*1j,
                       -rveSize*0.5:rveSize*0.5:n*1j,
                       -rveSize*0.5:rveSize*0.5:n*1j]

    # Convenient numbering of nodes (keeping track of all 27 nodes for the fixed grid is almost impossible otherwise).
    nX, nY, nZ = np.mgrid[0:n, 0:n, 0:n].astype(int)
    nC = nX + nY*n + nZ*n*n + 1  # type: np.ndarray

    # Create elements
    #print("Generating elements")
    eX, eY, eZ = np.mgrid[0:nelem, 0:nelem, 0:nelem].astype(int)
    e = eX + eY*nelem + eZ*nelem*nelem
    elem = np.zeros([nelem*nelem*nelem, w**3]).astype(int)
    if w == 3:
        elem[e, 0] = nC[0+eX*2, 0+eY*2, 2+eZ*2]
        elem[e, 1] = nC[0+eX*2, 2+eY*2, 2+eZ*2]
        elem[e, 2] = nC[2+eX*2, 2+eY*2, 2+eZ*2]
        elem[e, 3] = nC[2+eX*2, 0+eY*2, 2+eZ*2]

        elem[e, 4] = nC[0+eX*2, 0+eY*2, 0+eZ*2]
        elem[e, 5] = nC[0+eX*2, 2+eY*2, 0+eZ*2]
        elem[e, 6] = nC[2+eX*2, 2+eY*2, 0+eZ*2]
        elem[e, 7] = nC[2+eX*2, 0+eY*2, 0+eZ*2]

        elem[e,  8] = nC[0+eX*2, 1+eY*2, 2+eZ*2]
        elem[e,  9] = nC[1+eX*2, 2+eY*2, 2+eZ*2]
        elem[e, 10] = nC[2+eX*2, 1+eY*2, 2+eZ*2]
        elem[e, 11] = nC[1+eX*2, 0+eY*2, 2+eZ*2]

        elem[e, 12] = nC[0+eX*2, 1+eY*2, 0+eZ*2]
        elem[e, 13] = nC[1+eX*2, 2+eY*2, 0+eZ*2]
        elem[e, 14] = nC[2+eX*2, 1+eY*2, 0+eZ*2]
        elem[e, 15] = nC[1+eX*2, 0+eY*2, 0+eZ*2]

        elem[e, 16] = nC[0+eX*2, 0+eY*2, 1+eZ*2]
        elem[e, 17] = nC[0+eX*2, 2+eY*2, 1+eZ*2]
        elem[e, 18] = nC[2+eX*2, 2+eY*2, 1+eZ*2]
        elem[e, 19] = nC[2+eX*2, 0+eY*2, 1+eZ*2]

        elem[e, 20] = nC[1+eX*2, 1+eY*2, 2+eZ*2]
        elem[e, 21] = nC[1+eX*2, 1+eY*2, 0+eZ*2]
        elem[e, 22] = nC[0+eX*2, 1+eY*2, 1+eZ*2]
        elem[e, 23] = nC[1+eX*2, 2+eY*2, 1+eZ*2]
        elem[e, 24] = nC[2+eX*2, 1+eY*2, 1+eZ*2]
        elem[e, 25] = nC[1+eX*2, 0+eY*2, 1+eZ*2]

        elem[e, 26] = nC[1+eX*2, 1+eY*2, 1+eZ*2]
    else:
        elem[e, 0] = nC[0+eX, 0+eY, 1+eZ]
        elem[e, 1] = nC[0+eX, 1+eY, 1+eZ]
        elem[e, 2] = nC[1+eX, 1+eY, 1+eZ]
        elem[e, 3] = nC[1+eX, 0+eY, 1+eZ]
        elem[e, 4] = nC[0+eX, 0+eY, 0+eZ]
        elem[e, 5] = nC[0+eX, 1+eY, 0+eZ]
        elem[e, 6] = nC[1+eX, 1+eY, 0+eZ]
        elem[e, 7] = nC[1+eX, 0+eY, 0+eZ]
    
    f = open(folder + '/' + name, 'w')

    # Nodes:
    #X = X - rveSize*0.5
    #Y = Y - rveSize*0.5
    #Z = Z - rveSize*0.5
    for nz in range(n):
        for ny in range(n):
            for nx in range(n):
                node = nx + ny * n + nz * n * n + 1
                print('Node', node, 'coords 3', X[nx, ny, nz], Y[nx, ny, nz], Z[nx, ny, nz], file=f)

    # Elements:
    for ez in range(nelem):
        for ey in range(nelem):
            for ex in range(nelem):
                e = ex + ey*nelem + ez*nelem*nelem
                q = elem[e]
                print(elname, e+1, ' nodes', len(q), *q, file=f)

    xp = []  # x+
    for ez in range(0, nelem):
        for ey in range(0, nelem):
            ex = nelem-1
            e = ex + ey*nelem + ez*nelem*nelem + 1
            xp.append([e, 5])

    xm = []  # x-
    for ez in range(0, nelem):
        for ey in range(0, nelem):
            ex = 0
            e = ex + ey*nelem + ez*nelem*nelem + 1
            xm.append([e, 3])

    yp = []  # y+
    for ez in range(0, nelem):
        ey = nelem-1
        for ex in range(0, nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            yp.append([e, 4])

    ym = []  # y-
    for ez in range(0, nelem):
        ey = 0
        for ex in range(0, nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            ym.append([e, 6])

    zp = []  # z+
    ez = nelem-1
    for ey in range(0, nelem):
        for ex in range(0, nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            zp.append([e, 1])

    zm = []  # z-
    ez = 0
    for ey in range(0, nelem):
        for ex in range(0, nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            zm.append([e, 2])

    print('Set 1 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in xm]), file=f)
    print('Set 2 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in ym]), file=f)
    print('Set 3 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in zm]), file=f)
    print('Set 4 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in xp]), file=f)
    print('Set 5 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in yp]), file=f)
    print('Set 6 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in zp]), file=f)

    print('Set 7 nodes 1', (n*n*n + n*n + n) // 2, file=f)

    tot = np.concatenate([xm, ym, zm, xp, yp, zp])
    print('Set 8 elementboundaries', 6*2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in tot]), file=f)
    
    tot = np.concatenate([xp, yp, zp])
    print('Set 9 elementboundaries', 3*2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in tot]), file=f)
    tot = np.concatenate([xm, ym, zm])
    print('Set 10 elementboundaries', 3*2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in tot]), file=f)

def printRVE(folder, name, rveSampleNumber, rveSize, rvePosition, nelem, bctype, k):
    # elname, w = 'Hex21Stokes', 3
    # elname, w = 'Q27Space', 3
    # elname, w = 'LSpace', 2
    # elname, w = 'QBrick1HT', 3
    elname, w = 'Brick1HT', 2

    rveInclusions = rveToolbox.getInclusionsInBox(rvePosition, rveSize, inclusions)
    print('Generating RVE {:3d}, {:.1f}, {}, {:.0f}, {} inclusions in RVE'.format(rveSampleNumber, rveSize, bctype, k, len(rveInclusions)))

    n = nelem*(w-1) + 1
    #X, Y, Z = np.mgrid[-rveSize*0.5:rveSize*0.5:n*1j,
                       #-rveSize*0.5:rveSize*0.5:n*1j,
                       #-rveSize*0.5:rveSize*0.5:n*1j]

    # Convenient numbering of nodes (keeping track of all 27 nodes for the fixed grid is almost impossible otherwise).
    #nX, nY, nZ = np.mgrid[0:n, 0:n, 0:n].astype(int)
    #nC = nX + nY*n + nZ*n*n + 1  # type: np.ndarray

    # Create elements
    #eX, eY, eZ = np.mgrid[0:nelem, 0:nelem, 0:nelem].astype(int)
    #e = eX + eY*nelem + eZ*nelem*nelem
    #elem = np.zeros([nelem*nelem*nelem, w**3]).astype(int)
    #if w == 3:
        #elem[e, 0] = nC[0+eX*2, 0+eY*2, 2+eZ*2]
        #elem[e, 1] = nC[0+eX*2, 2+eY*2, 2+eZ*2]
        #elem[e, 2] = nC[2+eX*2, 2+eY*2, 2+eZ*2]
        #elem[e, 3] = nC[2+eX*2, 0+eY*2, 2+eZ*2]

        #elem[e, 4] = nC[0+eX*2, 0+eY*2, 0+eZ*2]
        #elem[e, 5] = nC[0+eX*2, 2+eY*2, 0+eZ*2]
        #elem[e, 6] = nC[2+eX*2, 2+eY*2, 0+eZ*2]
        #elem[e, 7] = nC[2+eX*2, 0+eY*2, 0+eZ*2]

        #elem[e,  8] = nC[0+eX*2, 1+eY*2, 2+eZ*2]
        #elem[e,  9] = nC[1+eX*2, 2+eY*2, 2+eZ*2]
        #elem[e, 10] = nC[2+eX*2, 1+eY*2, 2+eZ*2]
        #elem[e, 11] = nC[1+eX*2, 0+eY*2, 2+eZ*2]

        #elem[e, 12] = nC[0+eX*2, 1+eY*2, 0+eZ*2]
        #elem[e, 13] = nC[1+eX*2, 2+eY*2, 0+eZ*2]
        #elem[e, 14] = nC[2+eX*2, 1+eY*2, 0+eZ*2]
        #elem[e, 15] = nC[1+eX*2, 0+eY*2, 0+eZ*2]

        #elem[e, 16] = nC[0+eX*2, 0+eY*2, 1+eZ*2]
        #elem[e, 17] = nC[0+eX*2, 2+eY*2, 1+eZ*2]
        #elem[e, 18] = nC[2+eX*2, 2+eY*2, 1+eZ*2]
        #elem[e, 19] = nC[2+eX*2, 0+eY*2, 1+eZ*2]

        #elem[e, 20] = nC[1+eX*2, 1+eY*2, 2+eZ*2]
        #elem[e, 21] = nC[1+eX*2, 1+eY*2, 0+eZ*2]
        #elem[e, 22] = nC[0+eX*2, 1+eY*2, 1+eZ*2]
        #elem[e, 23] = nC[1+eX*2, 2+eY*2, 1+eZ*2]
        #elem[e, 24] = nC[2+eX*2, 1+eY*2, 1+eZ*2]
        #elem[e, 25] = nC[1+eX*2, 0+eY*2, 1+eZ*2]

        #elem[e, 26] = nC[1+eX*2, 1+eY*2, 1+eZ*2]
    #else:
        #elem[e, 0] = nC[0+eX, 0+eY, 1+eZ]
        #elem[e, 1] = nC[0+eX, 1+eY, 1+eZ]
        #elem[e, 2] = nC[1+eX, 1+eY, 1+eZ]
        #elem[e, 3] = nC[1+eX, 0+eY, 1+eZ]
        #elem[e, 4] = nC[0+eX, 0+eY, 0+eZ]
        #elem[e, 5] = nC[0+eX, 1+eY, 0+eZ]
        #elem[e, 6] = nC[1+eX, 1+eY, 0+eZ]
        #elem[e, 7] = nC[1+eX, 0+eY, 0+eZ]

    # Generate sets:
    #xp = []  # x+
    #for ez in range(0, nelem):
        #for ey in range(0, nelem):
            #ex = nelem-1
            #e = ex + ey*nelem + ez*nelem*nelem + 1
            #xp.append([e, 5])

    #xm = []  # x-
    #for ez in range(0, nelem):
        #for ey in range(0, nelem):
            #ex = 0
            #e = ex + ey*nelem + ez*nelem*nelem + 1
            #xm.append([e, 3])

    #yp = []  # y+
    #for ez in range(0, nelem):
        #ey = nelem-1
        #for ex in range(0, nelem):
            #e = ex + ey*nelem + ez*nelem*nelem + 1
            #yp.append([e, 4])

    #ym = []  # y-
    #for ez in range(0, nelem):
        #ey = 0
        #for ex in range(0, nelem):
            #e = ex + ey*nelem + ez*nelem*nelem + 1
            #ym.append([e, 6])

    #zp = []  # z+
    #ez = nelem-1
    #for ey in range(0, nelem):
        #for ex in range(0, nelem):
            #e = ex + ey*nelem + ez*nelem*nelem + 1
            #zp.append([e, 1])

    #zm = []  # z-
    #ez = 0
    #for ey in range(0, nelem):
        #for ex in range(0, nelem):
            #e = ex + ey*nelem + ez*nelem*nelem + 1
            #zm.append([e, 2])

    ccoord = np.mgrid[(0.5*rveSize/nelem + rvePosition[0]):(rvePosition[0]+rveSize-0.5*rveSize/nelem):nelem*1j,
                      (0.5*rveSize/nelem + rvePosition[1]):(rvePosition[1]+rveSize-0.5*rveSize/nelem):nelem*1j,
                      (0.5*rveSize/nelem + rvePosition[2]):(rvePosition[2]+rveSize-0.5*rveSize/nelem):nelem*1j]
    ccoord = np.transpose([ccoord[0].flatten('F'), ccoord[1].flatten('F'), ccoord[2].flatten('F')])

    emat = np.ones([nelem*nelem*nelem]).astype(int)
    for inc in rveInclusions:
        emat[np.linalg.norm(ccoord - inc[1], axis=1) <= inc[0]] = 2
    
    #print("Starting printing to file")
    # Write the input file
    f = open(folder + '/' + name + '.in', 'w')
    print(name+'.out', file=f)
    print('Hex grid RVE', file=f)
    if elname == 'Hexa21Stokes':
        print('StokesFlow nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nonlinform 1 nmodules 0', file=f)
        print('VTKXML tstep_all domain_all primvars 2 4 5 cellvars 1 103', file=f)
    elif elname == 'LSpace' or elname == 'Q27Space':
        print('StaticStructural nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nmodules 1', file=f)
        print('VTKXML tstep_all domain_all primvars 1 1 cellvars 3 103 1 4', file=f)
    else:
        print('StationaryProblem nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nmodules 1', file=f)
        print('VTKXML tstep_all domain_all primvars 1 6 cellvars 3 103 56 41', file=f)
    print('Domain 3d', file=f)
    print('OutputManager', file=f)
    print('ndofman', n*n*n, 'nelem', nelem*nelem*nelem, 'ncrosssect 2 nmat 2 nbc 2 nic 0 nltf 1 nset 12 nsd 3', file=f)

    # Nodes:
    #for nz in range(n):
        #for ny in range(n):
            #for nx in range(n):
                #node = nx + ny * n + nz * n * n + 1
                #print('Node', node, 'coords 3', X[nx, ny, nz], Y[nx, ny, nz], Z[nx, ny, nz], file=f)

    ## Elements:
    #for ez in range(nelem):
        #for ey in range(nelem):
            #for ex in range(nelem):
                #e = ex + ey*nelem + ez*nelem*nelem
                #q = elem[e]
                #print(elname, e+1, ' nodes', len(q), *q, file=f)

    print('@include "hex.grid"', file=f)

    #print('Set 1 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in xm]), file=f)
    #print('Set 2 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in ym]), file=f)
    #print('Set 3 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in zm]), file=f)
    #print('Set 4 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in xp]), file=f)
    #print('Set 5 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in yp]), file=f)
    #print('Set 6 elementboundaries', 2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in zp]), file=f)

    #print('Set 7 nodes 1', (n*n*n + n*n + n) // 2, file=f)

    #tot = np.concatenate([xm, ym, zm, xp, yp, zp])
    #print('Set 8 elementboundaries', 6*2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in tot]), file=f)
    
    #tot = np.concatenate([xp, yp, zp])
    #print('Set 9 elementboundaries', 3*2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in tot]), file=f)
    #tot = np.concatenate([xm, ym, zm])
    #print('Set 10 elementboundaries', 3*2*nelem*nelem, ' '.join([str(a) + ' ' + str(b) for a, b in tot]), file=f)
    
    print('Set 11 elements', np.count_nonzero(emat == 1), ' '.join([str(a+1) for a in np.nonzero(emat == 1)[0]]), file=f)
    print('Set 12 elements', np.count_nonzero(emat == 2), ' '.join([str(a+1) for a in np.nonzero(emat == 2)[0]]), file=f)

    if elname == 'Hexa21Stokes':
        print('FluidCS 1 mat 1 set 11', file=f)
        print('FluidCS 1 mat 2 set 12', file=f)
        print('NewtonianFluid 1 d 1 mu 1', file=f)
        print('NewtonianFluid 2 d {} mu 1'.format(k), file=f)
    elif elname == 'LSpace' or elname == 'Q27Space':
        print('SimpleCS 1 material 1 set 11', file=f)
        print('SimpleCS 1 material 2 set 12', file=f)
        print('IsoLE 1 d 1 E 1 n 0.3 talpha 0', file=f)
        print('IsoLE 2 d 1 E {} n 0.3 talpha 0'.format(k), file=f)
    else:
        print('SimpleTransportCS 1 mat 1 set 11', file=f)
        print('SimpleTransportCS 2 mat 2 set 12', file=f)
        print('IsoHeat 1 d 1 k 1 c 0', file=f)
        print('IsoHeat 2 d 1 k {} c 0'.format(k), file=f)

    if bctype == 'd' or bctype == 'md':
        # print('PrescribedGradient 1 loadTimeFunction 1 set 1 dofs 6 1 2 3 7 8 9 gradient 3 3 {0 0 1; 0 0 0; 1 0 0}', file=f)
        # print('MixedGradientPressureDirichlet 1 loadTimeFunction 1 set 1 dofs 6 1 2 3 7 8 9 devgradient 6 0 0 0 1 1 1 pressure 0', file=f)
        if bctype == 'md':
            print('TMGradDirichlet 1 loadTimeFunction 1 centercoords 3 0 0 0 gradient 3 0 0 1 dofs 1 10 set 8 surfsets 6 1 2 3 4 5 6 usexi', file=f)
        else:
            print('TMGradDirichlet 1 loadTimeFunction 1 centercoords 3 0 0 0 gradient 3 0 0 1 dofs 1 10 set 8', file=f)
        print('BoundaryCondition 2 loadTimeFunction 1 values 1 0 dofs 1 10 set 0', file=f)
    elif bctype == 'n' or bctype == 'mn':
        # print('MixedGradientPressureNeumann 1 loadTimeFunction 1 set 1 devgradient 6 0 0 0 1 1 1 pressure 0', file=f)
        if bctype == 'mn':
            print('TMGradNeumann 1 loadTimeFunction 1 centercoords 3 0 0 0 gradient 3 0 0 1 dofs 1 10 surfsets 6 1 2 3 4 5 6 useeta', file=f)
        else:
            print('TMGradNeumann 1 loadTimeFunction 1 centercoords 3 0 0 0 gradient 3 0 0 1 dofs 1 10 surfsets 6 1 2 3 4 5 6', file=f)
        print('BoundaryCondition 2 loadTimeFunction 1 values 1 0 dofs 1 10 set 0', file=f)
    elif bctype == 'p':
        print('TMGradPeriodic 1 loadTimeFunction 1 centercoords 3 0 0 0 gradient 3 0 0 1 jump 3 {0} {0} {0} dofs 1 10 set 9 masterset 10'.format(rveSize), file=f)
        print('BoundaryCondition 2 loadTimeFunction 1 values 1 0 dofs 1 10 set 0', file=f)
    else:
        print('BoundaryCondition 1 loadTimeFunction 1 values 1 0 dofs 1 10 set 1', file=f)
        print('BoundaryCondition 2 loadTimeFunction 1 values 1 1 dofs 1 10 set 4', file=f)

    print('ConstantFunction 1 f(t) 1.0', file=f)

    f.close()


seed = 0  # Controlled randomness
# Generate inclusions:
density = 0.20  # Minimum density
averageRadius = 1  # Average inclusion radius
boxSize = 100  # Max size of the domain with randomized particles. Computationally expensive to set it high.
inclusions = rveToolbox.generateSphericalInclusions(density, boxSize, averageRadius, averageRadius,
                                                    averageRadius*0.2, 3, seed)

import struct
outfile = open('inclusions.data', 'wb')
outfile.write( struct.pack('d', boxSize) )
outfile.write( struct.pack('i', len(inclusions) ) )
for inc in inclusions:
    outfile.write( struct.pack('d'*4, inc[1][0], inc[1][1], inc[1][2], inc[0])  );
outfile.close()

bcs = ['d', 'n', 'p', 'md', 'mn']
ks = [1e-3, 1e-1, 0.5, 2, 1e1, 1e3]
rveSizes = [5, 7.5, 10, 12.5, 15]
#rveSamples = [500, 250, 100, 50, 25]
rveSamples = [100, 75, 50, 25, 10]

# Generating examples:
#rveSizes = [10]
#rveSamples = [1]
#rvePosition = np.random.random(3) * (boxSize - 10)

for rveSize, samples in zip(rveSizes, rveSamples):
    nelem = int(8*rveSize)
    for k in ks:
        for bc in bcs:
            folder = 'SVEs_{:.0e}_{:.1f}_{}'.format(k, rveSize, bc)
            if not os.path.exists(folder):
                os.makedirs(folder)
            printMesh(folder, 'hex.grid', rveSize, nelem)
            for sample in range(samples):
                rvePosition = np.random.random(3) * (boxSize - rveSize)
                printRVE(folder, 'hex_sample{}'.format(sample), sample, rveSize, rvePosition, nelem, bc, k)
