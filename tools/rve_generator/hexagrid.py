#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-

# Description: Creates a regular cubic grid of 27 node elements. Serves as an example of how to generate an input file for OOFEM from Python.
# Author: Mikael Öhman
# License: CC0

from __future__ import print_function
import rveToolbox
from numpy import *

def printRVE(rveSampleNumber, rveSize, rvePosition):
    # 3 nodes along each element axis
    #elname, w = 'Hex21Stokes', 3
    #elname, w = 'Q27Space', 3
    #elname, w = 'LSpace', 2
    #elname, w = 'QBrick1HT', 3
    elname, w = 'Brick1HT', 2

    nelem = 5*rveSize; # 5 elements per unit cell

    rveInclusions = rveToolbox.getInclusionsInBox(rvePosition, rveSize, inclusions)
    print('%d inclusions in RVE'%len(rveInclusions))

    n = nelem*(w-1) + 1;
    X, Y, Z = mgrid[rvePosition[0]:(rvePosition[0]+rveSize):n*1j,
                    rvePosition[1]:(rvePosition[1]+rveSize):n*1j,
                    rvePosition[2]:(rvePosition[2]+rveSize):n*1j ]

    # Convenient numbering of nodes (keeping track of all 27 nodes for the fixed grid is almost impossible otherwise).
    nX, nY, nZ = mgrid[0:n, 0:n, 0:n ].astype(int)
    nC = nX + nY*n + nZ*n*n + 1; # + 1 because of oofems numbering

    # Create elements
    elem = zeros([nelem*nelem*nelem, w**2]).astype(int)
    emat = zeros([nelem*nelem*nelem]).astype(int)
    for ez in range(0,nelem):
        for ey in range(0,nelem):
            for ex in range(0,nelem):
                e = ex + ey*nelem + ez*nelem*nelem

                if w == 3:
                    elem[e, 0] = nC[0+ex*2,0+ey*2,2+ez*2]
                    elem[e, 1] = nC[0+ex*2,2+ey*2,2+ez*2]
                    elem[e, 2] = nC[2+ex*2,2+ey*2,2+ez*2]
                    elem[e, 3] = nC[2+ex*2,0+ey*2,2+ez*2]
            
                    elem[e, 4] = nC[0+ex*2,0+ey*2,0+ez*2]
                    elem[e, 5] = nC[0+ex*2,2+ey*2,0+ez*2]
                    elem[e, 6] = nC[2+ex*2,2+ey*2,0+ez*2]
                    elem[e, 7] = nC[2+ex*2,0+ey*2,0+ez*2]
            
                    elem[e, 8] = nC[0+ex*2,1+ey*2,2+ez*2]
                    elem[e, 9] = nC[1+ex*2,2+ey*2,2+ez*2]
                    elem[e,10] = nC[2+ex*2,1+ey*2,2+ez*2]
                    elem[e,11] = nC[1+ex*2,0+ey*2,2+ez*2]
            
                    elem[e,12] = nC[0+ex*2,1+ey*2,0+ez*2]
                    elem[e,13] = nC[1+ex*2,2+ey*2,0+ez*2]
                    elem[e,14] = nC[2+ex*2,1+ey*2,0+ez*2]
                    elem[e,15] = nC[1+ex*2,0+ey*2,0+ez*2]
            
                    elem[e,16] = nC[0+ex*2,0+ey*2,1+ez*2]
                    elem[e,17] = nC[0+ex*2,2+ey*2,1+ez*2]
                    elem[e,18] = nC[2+ex*2,2+ey*2,1+ez*2]
                    elem[e,19] = nC[2+ex*2,0+ey*2,1+ez*2]
            
                    elem[e,20] = nC[1+ex*2,1+ey*2,2+ez*2]
                    elem[e,21] = nC[1+ex*2,1+ey*2,0+ez*2]
                    elem[e,22] = nC[0+ex*2,1+ey*2,1+ez*2]
                    elem[e,23] = nC[1+ex*2,2+ey*2,1+ez*2]
                    elem[e,24] = nC[2+ex*2,1+ey*2,1+ez*2]
                    elem[e,25] = nC[1+ex*2,0+ey*2,1+ez*2]
            
                    elem[e,26] = nC[1+ex*2,1+ey*2,1+ez*2]
                else:
                    elem[e, 0] = nC[0+ex*1,0+ey*1,1+ez*1]
                    elem[e, 1] = nC[0+ex*1,1+ey*1,1+ez*1]
                    elem[e, 2] = nC[1+ex*1,1+ey*1,1+ez*1]
                    elem[e, 3] = nC[1+ex*1,0+ey*1,1+ez*1]
            
                    elem[e, 4] = nC[0+ex*1,0+ey*1,0+ez*1]
                    elem[e, 5] = nC[0+ex*1,1+ey*1,0+ez*1]
                    elem[e, 6] = nC[1+ex*1,1+ey*1,0+ez*1]
                    elem[e, 7] = nC[1+ex*1,0+ey*1,0+ez*1]


                # Check if center coordinate is within some sphere:
                emat[e] = 1
                ccoord = array([
                        X[1+ex*2,1+ey*2,1+ez*2],
                        Y[1+ex*2,1+ey*2,1+ez*2],
                        Z[1+ex*2,1+ey*2,1+ez*2]
                        ])
                for inc in rveInclusions:
                    if linalg.norm(ccoord - inc[1]) <= inc[0]:
                        emat[e] = 2
                        break


    # Write the input file
    fname = 'hex_size%d_sample%d'%(rveSize,rveSampleNumber);
    f = open(fname+'.in','w')
    print(fname+'.out', file=f)
    print('Hex grid RVE', file=f)
    if elname == 'Hexa21Stokes':
        print('StokesFlow nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nonlinform 1 nmodules 0', file=f)
        print('VTKXML tstep_all domain_all primvars 2 4 5 cellvars 1 46', file=f)
    elif elname == 'LSpace' or elname == 'Q27Space':
        print('StaticStructural nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nmodules 1', file=f)
        print('VTKXML tstep_all domain_all primvars 1 1 cellvars 2 46 4', file=f)
    else:
        print('StationaryProblem nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nmodules 1', file=f)
        print('VTKXML tstep_all domain_all primvars 1 1 cellvars 2 46 4', file=f)
    print('StaticStructural nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nmodules 1', file=f)
    print('VTKXML tstep_all domain_all primvars 1 1 cellvars 2 46 4', file=f)
    print('Domain 3d', file=f)
    print('OutputManager', file=f)
    print('ndofman', n*n*n, 'nelem', nelem*nelem*nelem, 'ncrosssect 1 nmat 2 nbc 4 nic 0 nltf 1 nset 4', file=f)

    # Nodes:
    for nz in range(0,n):
        for ny in range(0,n):
            for nx in range(0,n):
                node = nx + ny * n + nz * n * n + 1
                print('Node', node, 'coords 3', X[nx,ny,nz], Y[nx,ny,nz], Z[nx,ny,nz], file=f)


    # Elements:
    for ez in range(0,nelem):
        for ey in range(0,nelem):
            for ex in range(0,nelem):
                e = ex + ey*nelem + ez*nelem*nelem
                q = elem[e];
                print(elname, e+1, 'crosssect', emat[e], ' nodes',size(q), *q, file=f)

    if elname == 'Hexa21Stokes':
        print('FluidCS 1 mat 1', file=f)
        print('FluidCS 1 mat 2', file=f)
    elif elname == 'LSpace' or elname == 'Q27Space':
        print('SimpleCS 1 material 1', file=f)
        print('SimpleCS 1 material 2', file=f)
    else:
        print('SimpleTransportCS 1 mat 1', file=f)
        print('SimpleTransportCS 2 mat 2', file=f)

    print('IsoLE 1 d 1 E 1 n 0.3 talpha 0', file=f)
    print('IsoLE 2 d 1 E 5 n 0.3 talpha 0', file=f)
    #print('MixedGradientPressureNeumann 1 loadTimeFunction 1 set 1 devgradient 6 0 0 0 1 1 1 pressure 0', file=f)
    #print('MixedGradientPressureDirichlet 1 loadTimeFunction 1 set 1 dofs 6 1 2 3 7 8 9 devgradient 6 0 0 0 1 1 1 pressure 0', file=f)
    #print('PrescribedGradient 1 loadTimeFunction 1 set 1 dofs 6 1 2 3 7 8 9 gradient 3 3 {0 0 1; 0 0 0; 1 0 0}', file=f)
    print('BoundaryCondition 2 loadTimeFunction 1 values 3 0 0 0 dofs 3 1 2 3 set 0', file=f)
    print('BoundaryCondition 3 loadTimeFunction 1 values 3 0 0 0 dofs 3 1 2 3 set 0', file=f)
    print('ConstantFunction 1 f(t) 1.0', file=f)
    print('Set 1 elementboundaries', 2*6*nelem*nelem, end=' ', file=f)

    ez = 0
    for ey in range(0,nelem):
        for ex in range(0,nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            print(e, 2, end=' ', file=f)

    ez = nelem-1
    for ey in range(0,nelem):
        for ex in range(0,nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            print(e, 1, end=' ', file=f)

    for ez in range(0,nelem):
        ey = 0
        for ex in range(0,nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            print(e, 6, end=' ', file=f)

    for ez in range(0,nelem):
        ey = nelem-1
        for ex in range(0,nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            print(e, 4, end=' ', file=f)

    for ez in range(0,nelem):
        for ey in range(0,nelem):
            ex = 0
            e = ex + ey*nelem + ez*nelem*nelem + 1
            print(e, 3, end=' ', file=f)

    for ez in range(0,nelem):
        for ey in range(0,nelem):
            ex = nelem-1
            e = ex + ey*nelem + ez*nelem*nelem + 1
            print(e, 5, end=' ', file=f)

    print('', file=f)

    print('Set 2 elementboundaries', 2*nelem*nelem, end=' ', file=f)
    ez = nelem-1
    for ey in range(0,nelem):
        for ex in range(0,nelem):
            e = ex + ey*nelem + ez*nelem*nelem + 1
            print(e, 1, end=' ', file=f)
    print('', file=f)

    print('Set 3 nodes 1', nC[nelem][nelem][nelem], file=f)

    print('Set 4 elementRanges {(', 1, nelem*nelem*nelem,')}', file=f)

    f.close()




# Controlled randomness, specify the seed value
seed = 0
# Generate inclusions:
density = 0.1 # Minimum density
averageRadius = 0.5 # Average inclusion radius
boxSize = 50 # Max size of the domain with randomized particles. Computationally expensive to set it high.
inclusions = rveToolbox.generateSphericalInclusions(density, boxSize, averageRadius*2/3, averageRadius*4/3, averageRadius*0.05, 3, seed)

# Debugging code (testing to see that chosen sizes work well):
#rveToolbox.plotSphericalInclusions(inclusions, boxSize, 3, True)
#rveSize = 10;
#rvePosition = array([0,0,0])
#rveInclusions = rveToolbox.getInclusionsInBox(rvePosition, rveSize, inclusions)
#rveToolbox.plotSphericalInclusions(rveInclusions, boxSize, 3, True)


totalSamples = 200
rveSize = 3

for rveSampleNumber in range(totalSamples):
    print('Generating RVE %d'%rveSampleNumber)
    rvePosition = random.random(3) * (boxSize - rveSize)
    printRVE(rveSampleNumber, rveSize, rvePosition)



