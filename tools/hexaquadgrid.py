#!/usr/bin/python3

# Description: Creates a regular cubic grid of 27 node elements. Serves as an example of how to generate an input file for OOFEM from Python.
# Author: Mikael Öhman
# License: CC0

from numpy import *

nelem = 13;
n = nelem*2 + 1;

X, Y, Z = mgrid[0:nelem:n*1j, 0:nelem:n*1j, 0:nelem:n*1j ]

# Convenient numbering of nodes (keeping track of all 27 nodes for the fixed grid is almost impossible otherwise).
nX, nY, nZ = mgrid[0:n, 0:n, 0:n ].astype(int)
nC = nX + nY*n + nZ*n*n + 1; # + 1 because of oofems numbering

# Create elements
elem = zeros([nelem*nelem*nelem, 27]).astype(int)
emat = zeros([nelem*nelem*nelem]).astype(int)
for ez in range(0,nelem):
    for ey in range(0,nelem):
        for ex in range(0,nelem):
            e = ex + ey*nelem + ez*nelem*nelem

            emat[e] = 1;

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


# Write the input file
fname = 'hexstokes';
f = open(fname+'.in','w')
print(fname+'.out', file=f)
print('Hex grid for incompressible stokes flow', file=f)
print('StokesFlow nsteps 1 rtolf 1e-6 lstype 3 smtype 7 nonlinform 1 nmodules 1', file=f)
print('VTKXML tstep_all domain_all primvars 2 4 5', file=f)
print('Domain 2dIncompFlow', file=f)
print('OutputManager tstep_all dofman_all element_all', file=f)
print('ndofman', n*n*n, 'nelem', nelem*nelem*nelem, 'ncrosssect 1 nmat 2 nbc 2 nic 0 nltf 1', file=f)

# Nodes:
for nz in range(0,n):
    for ny in range(0,n):
        for nx in range(0,n):
            node = nx + ny * n + nz * n * n + 1
            no_p = mod(nx,2) or mod(ny,2) or mod(nz,2)
            bc = zeros([4 - no_p]).astype(int)
            if nz == 0:
                #if nx == 0 and ny == 0:
                bc[0] = 1
                bc[1] = 1
                bc[2] = 1
            elif nz == n-1:
                bc[0] = 1
                bc[1] = 1
                bc[2] = 2

            if size(bc) == 3:
                print('Node', node, 'coords 3', X[nx,ny,nz], Y[nx,ny,nz], Z[nx,ny,nz], 
                    'ndofs 3 dofidmask 3 7 8 9    bc 3', bc[0],bc[1],bc[2],file=f)
            else:
                print('Node', node, 'coords 3', X[nx,ny,nz], Y[nx,ny,nz], Z[nx,ny,nz], 
                    'ndofs 4 dofidmask 4 7 8 9 11 bc 4', bc[0],bc[1],bc[2],bc[3], file=f)


# Elements:
for ez in range(0,nelem):
    for ey in range(0,nelem):
        for ex in range(0,nelem):
            e = ex + ey*nelem + ez*nelem*nelem
            q = elem[e];
            print('Hexa21Stokes', e+1, 'mat', emat[e], 'crossSect 1 nodes',size(q),
                q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15],q[16],
                q[17],q[18],q[19],q[20],q[21],q[22],q[23],q[24],q[25],q[26], file=f)

print('EmptyCS 1', file=f)
print('NewtonianFluid 1 d 1 mu 1', file=f)
print('NewtonianFluid 2 d 1 mu 5', file=f)
print('BoundaryCondition 1 d 0.0 loadTimeFunction 1 defaultDofs 3 7 8 9', file=f)
print('BoundaryCondition 2 d 1.0 loadTimeFunction 1 defaultDofs 3 7 8 9', file=f)
print('ConstantFunction 1 f(t) 1.0', file=f)

f.close()
