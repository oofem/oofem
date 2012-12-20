/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "tr_dirshell.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "fei2dtrquad.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"


namespace oofem {
    FEI2dTrQuad TrDirShell :: interpolation_quad(1, 2);
    FEI3dTrQuad TrDirShell :: interpolation;
	IntArray TrDirShell :: ordering_x(18);
	IntArray TrDirShell :: ordering_m(18);
	IntArray TrDirShell :: ordering_gam(6);
    IntArray TrDirShell :: ordering_all(42);
	
	bool TrDirShell :: __initialized = TrDirShell :: initOrdering();

    FloatMatrix testmatrix(42,42);


    IntegrationRule **layerIntegrationRulesArray;

TrDirShell :: TrDirShell(int n, Domain *aDomain) : NLStructuralElement(n, aDomain),  LayeredCrossSectionInterface()
{
	this->numberOfIntegrationRules = 1;
	this->numberOfDofMans = 6;
    this->numberOfGaussPoints = 1;
    //this->computeGaussPoints();
}

IRResultType TrDirShell :: initializeFrom(InputRecord *ir)
{
    this->NLStructuralElement :: initializeFrom(ir);
	this->setupInitialNodeDirectors();
    return IRRT_OK;
}


int 
TrDirShell :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
	
	//FloatMatrix G;
	//interpolation.local2global(answer, lcoords, FEIElementGeometryWrapper(this));


	//this->evalContravarBaseVectors(G, lcoords, cellgeo, zeta)
	

    return 1;
    
}


void 
TrDirShell :: evalInitialCovarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3)
{
	int i;
	double x, y, z, Mx, My, Mz, zeta;
    FloatArray lcoords = *gp->giveCoordinates();
    zeta = giveLocalZetaCoord(gp);
	FloatArray N, M;
	FloatMatrix dNdxi, Mmat;
	
	// In plane base vectors
    this->interpolation.evaldNdxi(dNdxi, lcoords);

	G1.resize(3); G2.resize(3);
	G1.zero(); G2.zero();
	for ( i = 1; i <= 6; i++ ) {
		FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
		x = nodeI->at(1); y = nodeI->at(2); z = nodeI->at(3);

		M=this->giveInitialNodeDirector(i);
		Mx = M.at(1); My = M.at(2); Mz = M.at(3);

        G1.at(1) += dNdxi.at(i,1)*( x + zeta*Mx ); 
		G1.at(2) += dNdxi.at(i,1)*( y + zeta*My ); 
		G1.at(3) += dNdxi.at(i,1)*( z + zeta*Mz ); 
        G2.at(1) += dNdxi.at(i,2)*( x + zeta*Mx ); 
		G2.at(2) += dNdxi.at(i,2)*( y + zeta*My ); 
		G2.at(3) += dNdxi.at(i,2)*( z + zeta*Mz ); 
    }
	// Out of plane base vector = director
	this->evalInitialDirectorAt(gp, G3); // G3=M
}

void 
TrDirShell :: edgeEvalInitialCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &G1, FloatArray &G3)
{
	double x, y, z, Mx, My, Mz, zeta;
    FloatArray lcoords = *gp->giveCoordinates();
    zeta = 0.0; // no variation i z (yet)
	FloatArray N, M, dNdxi;
	FloatMatrix Mmat;
	IntArray edgeNodes;
    
    this->interpolation.computeLocalEdgeMapping(edgeNodes, iedge);
    this->interpolation.edgeEvaldNdxi(dNdxi, lcoords);

    // Base vector along edge
	G1.resize(3); G1.zero(); 
	for (int i = 1; i <= 3; i++ ) {
        FloatArray *nodeI = this->giveNode( edgeNodes.at(i) )->giveCoordinates();
		x = nodeI->at(1); y = nodeI->at(2); z = nodeI->at(3);
		M = this->giveInitialNodeDirector( edgeNodes.at(i) );
		Mx = M.at(1); My = M.at(2); Mz = M.at(3);
        G1.at(1) += dNdxi.at(i)*( x + zeta*Mx ); 
		G1.at(2) += dNdxi.at(i)*( y + zeta*My ); 
		G1.at(3) += dNdxi.at(i)*( z + zeta*Mz ); 
    }
	// Director will be the second base vector
	//this->evalInitialDirectorAt(gp, G3); // G3=M
    this->edgeEvalInitialDirectorAt(gp, G3, iedge);
}

void
TrDirShell :: evalInitialContravarBaseVectorsAt(GaussPoint * gp, FloatArray &Gcon1, FloatArray &Gcon2, FloatArray &Gcon3)
{	
	FloatArray Gcov1, Gcov2, Gcov3;
	this->evalInitialCovarBaseVectorsAt(gp, Gcov1, Gcov2, Gcov3);
	this->giveDualBase(Gcov1, Gcov2, Gcov3, Gcon1, Gcon2, Gcon3);
}

void
TrDirShell :: evalContravarBaseVectorsAt(GaussPoint *gp, FloatArray &gcon1, FloatArray &gcon2, FloatArray &gcon3, TimeStep *tStep, FloatArray &solVec)
{	
	FloatArray gcov1, gcov2, gcov3;
	this->evalCovarBaseVectorsAt(gp, gcov1, gcov2, gcov3, tStep, solVec);
	this->giveDualBase(gcov1, gcov2, gcov3, gcon1, gcon2, gcon3);
}

void
TrDirShell :: giveDualBase(const FloatArray &G1, const FloatArray &G2, const FloatArray &G3, FloatArray &g1, FloatArray &g2, FloatArray &g3 )
{	
	FloatMatrix gmat, ginv, test;
	
	gmat.resize(3,3);
	gmat.at(1,1) = G1.dotProduct(G1); gmat.at(1,2) = G1.dotProduct(G2); gmat.at(1,3) = G1.dotProduct(G3);
	gmat.at(2,2) = G2.dotProduct(G2); gmat.at(2,3) = G2.dotProduct(G3); gmat.at(3,3) = G3.dotProduct(G3);
	gmat.symmetrized();
	
	ginv.beInverseOf(gmat);
	g1.resize(3); g1.zero(); g2.resize(3); g2.zero(); g3.resize(3); g3.zero();

	g1.add(ginv.at(1,1),G1); g1.add(ginv.at(1,2),G2); g1.add(ginv.at(1,3),G3); 
	g2.add(ginv.at(2,1),G1); g2.add(ginv.at(2,2),G2); g2.add(ginv.at(2,3),G3);
	g3.add(ginv.at(3,1),G1); g3.add(ginv.at(3,2),G2); g3.add(ginv.at(3,3),G3);
}

void 
TrDirShell :: evalInitialDirectorAt(GaussPoint *gp, FloatArray &answer)
{	// Interpolates between the node directors
	//FloatArray lcoords = *((*gp).giveCoordinates()); // old
    
    FloatArray &lcoords = *gp->giveCoordinates();
	
    FloatArray N;
	this->interpolation_quad.evalN(N, lcoords, FEIElementGeometryWrapper(this));

	answer.resize(3); answer.zero();
	for (int i = 1; i <= 6; i++ ) {
		answer.add( N.at(i), this->giveInitialNodeDirector(i) ); 
	}
    
}

void 
TrDirShell :: edgeEvalInitialDirectorAt(GaussPoint *gp, FloatArray &answer, const int iEdge)
{	// Interpolates between the node directors
	//FloatArray lcoords = *((*gp).giveCoordinates()); // old
    
    FloatArray &lcoords = *gp->giveCoordinates();
	
    FloatArray N;
    IntArray edgeNodes;
    this->interpolation.computeLocalEdgeMapping(edgeNodes, iEdge);
	this->interpolation_quad.edgeEvalN(N, lcoords, FEIElementGeometryWrapper(this));

	answer.resize(3); answer.zero();
	for (int i = 1; i <= 3; i++ ) {
        answer.add( N.at(i), this->giveInitialNodeDirector( edgeNodes.at(i) ) ); 
	}
}

void 
TrDirShell :: setupInitialNodeDirectors()
{	// If the directors are not present in the input file, then they should be approximated as the normal to the initial surface.
	FloatMatrix dNdxi;
	FloatArray M, G1, G2, lcoords, nodeLocalXiCoords, nodeLocalEtaCoords;

	// Compute directors as normals to the surface

	// Set the local coordinates for the element nodes
	nodeLocalXiCoords.setValues( 6, 1., 0., 0., .5, 0., .5); // corner nodes then midnodes, uncertain of node numbering
	nodeLocalEtaCoords.setValues(6, 0., 1., 0., .5, .5, 0.);
	lcoords.resize(2); G1.resize(3); G2.resize(3); M.resize(3);		

    //double thickness = this->giveCrossSection()->give(CS_Thickness);
    double thickness = 1.0;
	for (int node = 1; node <= 6; node++ ) {

		this->initialNodeDirectors[node-1].resize(3);
		this->initialNodeDirectors[node-1].zero();

		lcoords.at(1) = nodeLocalXiCoords.at(node); 
		lcoords.at(2) = nodeLocalEtaCoords.at(node);
        this->interpolation.evaldNdxi(dNdxi, lcoords);
		
		G1.zero(); G2.zero(); M.zero();
		for ( int i = 1; i <= 6; i++){	// base vectors of the initial surface
			FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
			G1.add( dNdxi.at(i,1), *nodeI );
			G2.add( dNdxi.at(i,2), *nodeI );
		}
		M.beVectorProductOf(G1,G2);
		M.normalize();
        M.times(thickness); // Initialize M with constant thickness
		this->initialNodeDirectors[node-1].add(M); 	
	}
	
}

void
TrDirShell :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, TimeStep *tStep, FloatArray &genEps)
{
	FloatArray lcoords = *gp->giveCoordinates();
	double zeta = giveLocalZetaCoord(gp);
	
	FloatArray dxdxi, dxdxi1, dxdxi2, m, dmdxi, dmdxi1, dmdxi2, dgamdxi,  test;
    double dgamdxi1, dgamdxi2, gam;
    this->giveGeneralizedStrainComponents(genEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m, dgamdxi1, dgamdxi2, gam);

	double fac1 = ( zeta + 0.5*gam*zeta*zeta );
	double fac2 = ( 0.5*zeta*zeta );
	double fac3 = ( 1.0 + zeta*gam );

    g1.resize(3), g2.resize(3), g3.resize(3);
    g1.zero(); g2.zero(); g3.zero();
    g1.add(dxdxi1); g1.add(fac1,dmdxi1); g1.add(fac2*dgamdxi1,m);
    g2.add(dxdxi2); g2.add(fac1,dmdxi2); g2.add(fac2*dgamdxi2,m);
    g3.add(fac3,m);

}

void
TrDirShell :: edgeEvalCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &g1, FloatArray &g3, TimeStep *tStep)
{
	double zeta=0.0, gam;
	FloatArray lcoords = *gp->giveCoordinates();
	FloatArray a;
	FloatMatrix B;
    IntArray edgeNodes;

    this->interpolation.computeLocalEdgeMapping(edgeNodes, iedge);
	this->edgeComputeBmatrixAt(gp, B, 1, ALL_STRAINS);
	this->edgeGiveUpdatedSolutionVector(a, iedge, tStep);
	FloatArray eps;	      // generalized strain
	eps.beProductOf(B,a); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

	
	FloatArray dxdxi, m, dmdxi, test;
    double dgamdxi;
	dxdxi.setValues(3, eps.at(1), eps.at(2), eps.at(3) );
    dmdxi.setValues(3, eps.at(4), eps.at(5), eps.at(6) );
	m.setValues(3, eps.at(7), eps.at(8), eps.at(9) );
    dgamdxi = eps.at(10);
    gam = eps.at(11);

	g1.resize(3),  g3.resize(3);
	double fac1 = ( zeta + 0.5*gam*zeta*zeta );
	double fac2 = ( 0.5*zeta*zeta );
	double fac3 = ( 1.0 + zeta*gam );
	
	g1.at(1) = dxdxi.at(1) + fac1*dmdxi.at(1) + fac2*m.at(1)*dgamdxi;
	g1.at(2) = dxdxi.at(2) + fac1*dmdxi.at(2) + fac2*m.at(2)*dgamdxi;
	g1.at(3) = dxdxi.at(3) + fac1*dmdxi.at(3) + fac2*m.at(3)*dgamdxi;

	g3.at(1) = fac3*m.at(1); 
	g3.at(2) = fac3*m.at(2); 
	g3.at(3) = fac3*m.at(3);

}

void
TrDirShell :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer.setValues(7, D_u, D_v, D_w, w_u, w_v, w_w, gam);
}

void
TrDirShell :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at stepN.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    
	/*
	double dens, dV, load;
    GaussPoint *gp = NULL;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, stepN, mode);

    if ( force.giveSize() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        dens = this->giveMaterial()->give('d', gp);
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness);

        answer.resize(18);
        answer.zero();

        load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(7) = load;
        answer.at(13) = load;

        load = force.at(2) * dens * dV / 3.0;
        answer.at(2) = load;
        answer.at(8) = load;
        answer.at(14) = load;

        load = force.at(3) * dens * dV / 3.0;
        answer.at(3) = load;
        answer.at(9) = load;
        answer.at(15) = load;

        // transform result from global cs to local element cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }
    } else {
        answer.resize(0);          // nil resultant
    }
	*/
	answer.resize(0);
}


void
TrDirShell :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    /*
	int i, j;
    GaussPoint *gp;
    FloatArray v;

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

	
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        for ( j = 0; j < integrationRulesArray [ i ]->getNumberOfIntegrationPoints(); j++ ) {
            gp = integrationRulesArray [ i ]->getIntegrationPoint(j);

            // gp   -> printOutputAt(file,stepN) ;
			fprintf( file, "  GP %2d.%-2d :", i + 1, gp->giveNumber() );

            this->giveIPValue(v, gp, IST_ShellStrainCurvatureTensor, tStep);
            fprintf(file, "  strains ");
            fprintf( file,
                    " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                    v.at(1), v.at(2), v.at(3),  2. * v.at(4), 2. * v.at(5), 2. * v.at(6),
                    v.at(7), v.at(8), v.at(9),  2. * v.at(10), 2. * v.at(11), 2. * v.at(12) );

            this->giveIPValue(v, gp, IST_ShellForceMomentumTensor, tStep);
            fprintf(file, "\n              stresses");
            fprintf( file,
                    " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                    v.at(1), v.at(2), v.at(3),  v.at(4), v.at(5), v.at(6),
                    v.at(7), v.at(8), v.at(9),  v.at(10), v.at(11), v.at(12) );

            fprintf(file, "\n");
        }
    }
	*/
}


void
TrDirShell :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
	answer.resize(7, 42); answer.zero();
    FloatArray &lcoords = *gp->giveCoordinates();
	FloatArray N;
    this->interpolation.evalN(N, lcoords, FEIElementGeometryWrapper(this));

	/*    18   18    6
  	   3 [N_x   0    0
	   3   0   N_m   0
	   1   0    0  N_gmm ]
    */
    int i, j;
	for( i = 1, j = 0; i<=6; i++, j+=7  ){ 
		answer.at(1,1+j) = N.at(i);
		answer.at(2,2+j) = N.at(i);
		answer.at(3,3+j) = N.at(i);
		answer.at(4,4+j) = N.at(i);
		answer.at(5,5+j) = N.at(i);
		answer.at(6,6+j) = N.at(i);
        answer.at(7,7+j) = N.at(i);
	}
}

void
TrDirShell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li , int ui)
/* Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
   B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
*/
{
	answer.resize(18, 42); answer.zero();
    FloatArray lcoords = *gp->giveCoordinates();
	FloatArray N;
	FloatMatrix dNdxi;

    this->interpolation.evalN(N, lcoords, FEIElementGeometryWrapper(this));
    this->interpolation.evaldNdxi(dNdxi, lcoords);

	/*    18   18   6
  	   6 [B_u   0   0
	   6   0   B_w  0
	   3   0   N_w  0
	   2   0    0  B_gam 
	   1   0    0  N_gam] 
	*/
	  
    int i, j, pos;
	// First row
	for( i = 1, j = 0; i<=6; i++, j+=7  ){ 
		answer.at(1,1+j) = dNdxi.at(i,1);
		answer.at(2,2+j) = dNdxi.at(i,1);
		answer.at(3,3+j) = dNdxi.at(i,1);
		answer.at(4,1+j) = dNdxi.at(i,2);
		answer.at(5,2+j) = dNdxi.at(i,2);
		answer.at(6,3+j) = dNdxi.at(i,2);
	}

	// Second row
	pos =18;
	for( i = 1, j = 0; i<=6; i++, j+=7  ){ 
		answer.at(7,4+j) = dNdxi.at(i,1);
		answer.at(8,5+j) = dNdxi.at(i,1);
		answer.at(9,6+j) = dNdxi.at(i,1);
		answer.at(10,4+j) = dNdxi.at(i,2);
		answer.at(11,5+j) = dNdxi.at(i,2);
		answer.at(12,6+j) = dNdxi.at(i,2);
		answer.at(13,4+j) = N.at(i);
		answer.at(14,5+j) = N.at(i);
		answer.at(15,6+j) = N.at(i);
	}

	// Third row
	pos =36;
	for( i = 1, j = 0; i<=6; i++, j+=7  ){ 
		answer.at(16,7+j) = dNdxi.at(i,1);
		answer.at(17,7+j) = dNdxi.at(i,2);
		answer.at(18,7+j) = N.at(i);
	}

}

void
TrDirShell :: edgeComputeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint along one edge
{
	answer.resize(7, 21);
    answer.zero();

    FloatArray lcoords = *gp->giveCoordinates();
	FloatArray N;
	FloatMatrix dNdxi;
	this->interpolation.edgeEvalN(N, lcoords, FEIElementGeometryWrapper(this));


	/*    9   9    3
  	   3 [N_x   0    0
	   3   0   N_m   0
	   1   0    0  N_gmm ]
    */
    int i, j;
	for( i = 1, j = 0; i<=3; i++, j+=7  ){ 
		answer.at(1,1+j) = N.at(i);
		answer.at(2,2+j) = N.at(i);
		answer.at(3,3+j) = N.at(i);
		answer.at(4,4+j) = N.at(i);
		answer.at(5,5+j) = N.at(i);
		answer.at(6,6+j) = N.at(i);
        answer.at(7,7+j) = N.at(i);
	}

	
}

void
TrDirShell :: edgeComputeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li , int ui)
/* Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
   B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
*/
{
	answer.resize(11, 21); answer.zero();
    FloatArray lcoords = *gp->giveCoordinates();
	FloatArray N, dNdxi;
	//FloatMatrix dNdxi;

    this->interpolation.edgeEvalN(N, lcoords, FEIElementGeometryWrapper(this));
    this->interpolation.edgeEvaldNdxi(dNdxi, lcoords);

	/*     9    9   3
  	   3 [B_u   0   0
	   3   0   B_w  0
	   3   0   N_w  0
	   1   0    0  B_gam 
	   1   0    0  N_gam] 
	*/
	  
    int i, j;
	// First row
	for( i = 1, j = 0; i<=3; i++, j+=7  ){ 
		answer.at(1,1+j) = dNdxi.at(i);
		answer.at(2,2+j) = dNdxi.at(i);
		answer.at(3,3+j) = dNdxi.at(i);

	}

	// Second row
	for( i = 1, j = 0; i<=3; i++, j+=7  ){ 
		answer.at(4,4+j) = dNdxi.at(i);
		answer.at(5,5+j) = dNdxi.at(i);
		answer.at(6,6+j) = dNdxi.at(i);
		answer.at(7,4+j) = N.at(i);
		answer.at(8,5+j) = N.at(i);
		answer.at(9,6+j) = N.at(i);
	}

	// Third row
	for( i = 1, j = 0; i<=3; i++, j+=7  ){ 
		answer.at(10,7+j) = dNdxi.at(i);
		answer.at(11,7+j) = N.at(i);
	}

}




void 
TrDirShell :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        
        int nPointsTri = 6;
        int nPointsThickness = 2;
        int nPointsEdge = 2;
        integrationRulesArray = new IntegrationRule * [ 3 ];
        // Midplane and thickness
		integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
        integrationRulesArray[0]->SetUpPointsOnWedge2(nPointsTri, nPointsThickness, _3dMat);

        // Midplane only (Mass matrix integrated analytically through the thickness)
        integrationRulesArray [ 1 ] = new GaussIntegrationRule(1, this);
        integrationRulesArray[1]->SetUpPointsOnWedge2(nPointsTri, 1, _3dMat); 

        // Edge
        integrationRulesArray [ 2 ] = new GaussIntegrationRule(1, this);
        integrationRulesArray[2]->SetUpPointsOnLine2(nPointsEdge, _3dMat); 


        // Layered cross section
        LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >(TrDirShell::giveCrossSection());
        if( layeredCS == NULL ){
            OOFEM_ERROR("TrDirShell only supports layered cross section");
        }
        int numberOfLayers = layeredCS->give(CS_NumLayers);
        layerIntegrationRulesArray = new IntegrationRule * [ numberOfLayers ];
        for( int i = 1; i <= numberOfLayers; i++ ){
            layerIntegrationRulesArray[ i-1 ]= new GaussIntegrationRule(1, this);
            layerIntegrationRulesArray[ i-1 ]->SetUpPointsOnWedge2(nPointsTri, layeredCS->giveNumIntegrationPointsInLayer(), _3dMat); 
        }

    }

}


double 
TrDirShell :: giveLocalZetaCoord(GaussPoint *gp){

    return (*gp->giveCoordinates()).at(3) * this->giveCrossSection()->give(CS_Thickness)*0.5 ;
    //return (*gp->giveCoordinates()).at(3)  ;

}

double 
TrDirShell :: giveLayerZetaCoord(GaussPoint *gp, int layer){

    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >(TrDirShell::giveCrossSection());
    
    double layerMidCoord  = layeredCS->giveLayerMidZ(layer);
    double layerThickness = layeredCS->giveLayerThickness(layer);

    return layerMidCoord + layerThickness*(*gp->giveCoordinates()).at(3)*0.5;
   

}

FloatArray 
TrDirShell :: giveInitialNodeDirector(int i){
    FloatArray temp;
    temp = this->initialNodeDirectors[i-1];
    //temp.times(this->giveCrossSection()->give(CS_Thickness)*0.5);
    return temp;

}



void 
TrDirShell :: giveGeneralizedStrainComponents(FloatArray genEps, FloatArray &dphidxi1, FloatArray &dphidxi2, FloatArray &dmdxi1, 
         FloatArray &dmdxi2, FloatArray &m, double &dgamdxi1, double &dgamdxi2, double &gam){
        
        // generealized strain vector  [dxdxi, dmdxi, m, dgamdxi, gam]^T
        dphidxi1.setValues( 3, genEps.at(1) , genEps.at(2) , genEps.at(3)  );
        dphidxi2.setValues( 3, genEps.at(4) , genEps.at(5) , genEps.at(6)  );
		dmdxi1.setValues(   3, genEps.at(7) , genEps.at(8) , genEps.at(9)  );
		dmdxi2.setValues(   3, genEps.at(10), genEps.at(11), genEps.at(12) );
		     m.setValues(   3, genEps.at(13), genEps.at(14), genEps.at(15) );
		dgamdxi1 = genEps.at(16);
		dgamdxi2 = genEps.at(17);
		     gam = genEps.at(18);


}


// Tangent matrices

void
TrDirShell :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
   
    this->computeBulkTangentMatrix(answer, rMode, tStep);
        
    
    // Add contribution due to pressure load
    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    
    for ( int i = 1; i <= nLoads; i++ ) {  // For each pressre load that is applied
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int iSurf = this->boundaryLoadArray.at(2 * i); // load_id
        Load *load;
        load = this->domain->giveLoad(load_number);

        if ( load->giveClassID() == ConstantPressureLoadClass ){
            FloatMatrix K_pressure;
            this->computePressureTangentMatrix(K_pressure, load, iSurf, tStep);
            //answer.subtract(K_pressure);
            answer.add(K_pressure);
        }
    }
    
}

void
TrDirShell :: computeBulkTangentMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    GaussPoint *gp;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [0];

    double dV, a, b, c;
    FloatMatrix B, Bt, BtL, BtLB, L(18,18);
    FloatArray lcoords, BF;

    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec,tStep);
    

#if 0
    answer.resize(42,42); answer.zero();
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = *gp->giveCoordinates();
        double zeta = giveLocalZetaCoord(gp);

        FloatArray g1, g2, g3, S1g(3), S2g(3), S3g(3), m(3), dm1(3), dm2(3), temp1, temp2;
        double gam, dg1, dg2;

		this->computeBmatrixAt(gp, B);
        FloatArray genEps;	      
		genEps.beProductOf(B,solVec); // [dxdxi, dmdxi, m, dgamdxi, gam]^T
        this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);

        // Material stiffness
        FloatMatrix A[3][3];
        TrDirShell :: computeLinearizedStiffness(gp, tStep, S1g, S2g, S3g, A, genEps);
        
        // Tangent stiffness
        FloatArray   f4(3), f5(3), f3[2],temp(3), t1(3), t2(3), t3(3);
        FloatMatrix F1[3], F1T[2], F2, F2T, F2f5, temp3(3,3);

        // thickness coefficients
        a = zeta + 0.5*gam*zeta*zeta;  b = 0.5*zeta*zeta;  c = 1. + gam*zeta;

        // f1(alpha) = b*A(alpha,beta)*dg(beta) + c*A(alpha,3); 
        F1[0].add(dg1,A[0][0]); F1[0].add(dg2,A[0][1]); F1[0].times(b); F1[0].add(c,A[0][2]);
        F1[1].add(dg1,A[1][0]); F1[1].add(dg2,A[1][1]); F1[1].times(b); F1[1].add(c,A[1][2]);

        F2.add(dg1,A[2][0]);    F2.add(dg2,A[2][1]);    F2.times(b);    F2.add(c,A[2][2]);

        F1T[0].beTranspositionOf(F1[0]);
        F1T[1].beTranspositionOf(F1[1]);


        //f2 = c*A(3,3) + b*dg(alpha)*A(alpha,3)        
        F2f5.add(dg1,A[0][2]);    F2f5.add(dg2,A[1][2]);    F2f5.times(b);    F2f5.add(c,A[2][2]);
        F2T.beTranspositionOf(F2);

        //f3(alpha) = b*A(alpha,beta)*dm(beta) + zeta*A(alpha,3)*m; 
        t1.beProductOf(A[0][0],dm1); t2.beProductOf(A[0][1],dm2); t3.beProductOf(A[0][2],m);
        f3[0].add(b,t1); f3[0].add(b,t2); f3[0].add(zeta,t3);

        //f32 = b*( A.at(2,1)*dm1 + A.at(2,2)*dm2 )  + zeta*A.at(2,3)*m;
        t1.beProductOf(A[1][0],dm1); t2.beProductOf(A[1][1],dm2); t3.beProductOf(A[1][2],m);
        f3[1].add(b,t1); f3[1].add(b,t2); f3[1].add(zeta,t3);


        //f4 = b*dm(alpha)*A(alpha,3) + zeta*A(3,3)*m; 
        t1.beTProductOf(A[0][2],dm1); t2.beTProductOf(A[1][2],dm2); t3.beTProductOf(A[2][2],m);
        f4.add(b,t1); f4.add(b,t2);  f4.add(zeta,t3);

        // f5 = b*F1(alpha)*dm(alpha) + zeta*F2*m + zeta*N3
        t1.beProductOf(F1T[0],dm1); t2.beProductOf(F1T[1],dm2); t3.beProductOf(F2T,m);
        f5.add(b,t1); f5.add(b,t2); f5.add(zeta,t3); f5.add(zeta,S3g);


 /*
        
        A     a*A            F1           b*A*m             f3
       a*A    a^2*A          a*F1         a*b*A*m         a*f3+b*N
       F1     a*F1     b*F1*dgam+c*F2     b*(F1*m+N)        f5
       b*A*m  a*b*A*m    b*(F1*m+N)       b^2*m*A*m     b*m*f3
       f3     a*f3+b*N       f5           b*m*f3    b*dm*f3+zeta*m*f4
       */

        // L(1,1) = A(alpha,beta)
        FloatMatrix L11(6,6); L11.zero();
        L11.addSubMatrix(A[0][0],1,1); L11.addSubMatrix(A[0][1],1,4);
        L11.addSubMatrix(A[1][0],4,1); L11.addSubMatrix(A[1][1],4,4);

        // L(1,2) = a*A(alpha,beta)
        FloatMatrix L12(6,6); L12.zero(); L12.add(a,L11);
        
        // L(1,3) = F1(alpha)
        FloatMatrix L13(6,3); L13.zero();
        L13.addSubMatrix(F1[0],1,1); L13.addSubMatrix(F1[1],4,1);
        
        // L(1,4) = b*A*m
        FloatMatrix L14(6,2); L14.zero();
        t1.beProductOf(A[0][0],m); t2.beProductOf(A[0][1],m);
        L14.addSubVectorCol(t1,1,1); L14.addSubVectorCol(t2,1,2);
        t1.beProductOf(A[1][0],m); t2.beProductOf(A[1][1],m);
        L14.addSubVectorCol(t1,4,1); L14.addSubVectorCol(t2,4,2);
        L14.times(b);

        // L(1,5) = f3(alpha)
        FloatMatrix L15(6,1); L15.zero();
        L15.addSubVectorCol(f3[0],1,1); L15.addSubVectorCol(f3[1],4,1);


        // L(2,2) = a^2*A(alpha,beta)= a*L12
        FloatMatrix L22(6,6); L22.zero(); L22.add(a,L12);
        
        // L(2,3) = a*F1(alpha) = a*L13
        FloatMatrix L23(6,3); L23.zero(); L23.add(a,L13);
        
        // L(2,4) = a*F1(alpha)=a*L14
        FloatMatrix L24(6,2); L24.zero(); L24.add(a,L14);

        // L(2,5) = a*f3(alpha) + b*N(alpha)
        FloatMatrix L25(6,1); L25.zero(); 
        L25.addSubVectorCol(S1g,1,1); L25.addSubVectorCol(S2g,4,1); L25.times(b); L25.add(a,L15);
        

        // L(3,3) = b*F1(beta)*dgam(beta) + c*F2
        FloatMatrix L33(3,3); L33.zero();
        L33.add(c,F2T); temp3.add(dg1,F1T[0]); temp3.add(dg2,F1T[1]); temp3.times(b); L33.add(temp3);

                
        // L(3,4) = b*( F1(beta)*m + N(beta) )
        FloatMatrix L34(3,2); L34.zero(); 
        t1.beTProductOf(F1[0],m); t1.add(S1g); t1.times(b); 
        t2.beTProductOf(F1[1],m); t2.add(S2g); t2.times(b);
        L34.addSubVectorCol(t1,1,1); L34.addSubVectorCol(t2,1,2);

        // L(3,5) = f5
        FloatMatrix L35(3,1); L35.zero(); L35.addSubVectorCol(f5,1,1);


        // L(4,4) = b^2*m*A*m (2x2)
        FloatMatrix L44(2,2); L44.zero(); 
        t1.beProductOf(A[0][0],m); L44.at(1,1) = t1.dotProduct(m);
        t1.beProductOf(A[0][1],m); L44.at(1,2) = t1.dotProduct(m);
        t1.beProductOf(A[1][0],m); L44.at(2,1) = t1.dotProduct(m);
        t1.beProductOf(A[1][1],m); L44.at(2,2) = t1.dotProduct(m);
        L44.times(b*b);

        // L(4,5) = b*m*f3(alpha)
        FloatMatrix L45(2,1); 
        L45.at(1,1) = b*m.dotProduct(f3[0]);
        L45.at(2,1) = b*m.dotProduct(f3[1]);

        // L(5,5) = b*m*f3(alpha)
        FloatMatrix L55(1,1); 
        L55.at(1,1) = b*dm1.dotProduct(f3[0]) + b*dm2.dotProduct(f3[1]) + zeta*m.dotProduct(f4);
        
        FloatMatrix  L21, L31, L41, L43, L51, L32, L42, L52, L54, L53;
        L21.beTranspositionOf(L12); L31.beTranspositionOf(L13); L51.beTranspositionOf(L15); L32.beTranspositionOf(L23);
        L42.beTranspositionOf(L24); L43.beTranspositionOf(L34); L52.beTranspositionOf(L25); L54.beTranspositionOf(L45);
        L53.beTranspositionOf(L35); L41.beTranspositionOf(L14);

        // Assemble into L
        L.zero();
        L.addSubMatrix(L11,1,1);  
        L.addSubMatrix(L22,7,7);  
        L.addSubMatrix(L33,13,13); 
        L.addSubMatrix(L44,16,16);
        L.addSubMatrix(L55,18,18); 
        
        L.addSubMatrix(L12,1,7);   L.addSubMatrix(L21,7,1);       
        L.addSubMatrix(L13,1,13);  L.addSubMatrix(L31,13,1); 
        L.addSubMatrix(L14,1,16);  L.addSubMatrix(L41,16,1);  
        L.addSubMatrix(L15,1,18);  L.addSubMatrix(L51,18,1); 
        L.addSubMatrix(L23,7,13);  L.addSubMatrix(L32,13,7); 
        L.addSubMatrix(L24,7,16);  L.addSubMatrix(L42,16,7);
        L.addSubMatrix(L25,7,18);  L.addSubMatrix(L52,18,7);
        L.addSubMatrix(L34,13,16); L.addSubMatrix(L43,16,13);
        L.addSubMatrix(L35,13,18); L.addSubMatrix(L53,18,13);
        L.addSubMatrix(L45,16,18); L.addSubMatrix(L54,18,16); 
       
        //L.symmetrized();


        // Tangent matrix (K = B^T*L*B*dV)
        BtL.beTProductOf(B,L);
        BtLB.beProductOf(BtL,B);
        dV = this->computeVolumeAround(gp);
        BtLB.times(dV);
        answer.add(BtLB);

        
    }

    FloatMatrix test;
    test.beSubMatrixOf(L,1,6,1,6);
    printf("\n old \n");
    test.printYourself();
#endif



    answer.resize(42,42); answer.zero();
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >(this->giveCrossSection());

    int numberOfLayers = layeredCS->give(CS_NumLayers);

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [layer-1]; 
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

	    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);
            lcoords = *gp->giveCoordinates();

            double zeta = giveLayerZetaCoord(gp, layer);

            FloatArray g1, g2, g3, S1g(3), S2g(3), S3g(3), m(3), dm1(3), dm2(3), temp1, temp2;
            double gam, dg1, dg2;

		    this->computeBmatrixAt(gp, B);
            FloatArray genEps;	      
		    genEps.beProductOf(B,solVec); // [dxdxi, dmdxi, m, dgamdxi, gam]^T
            this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);

            // Material stiffness
            FloatMatrix A[3][3];
            TrDirShell :: computeLinearizedStiffness(gp, mat, tStep, S1g, S2g, S3g, A, genEps);
        
            // Tangent stiffness
            FloatArray   f4(3), f5(3), f3[2],temp(3), t1(3), t2(3), t3(3);
            FloatMatrix F1[3], F1T[2], F2, F2T, F2f5, temp3(3,3);

            // thickness coefficients
            a = zeta + 0.5*gam*zeta*zeta;  b = 0.5*zeta*zeta;  c = 1. + gam*zeta;

            // f1(alpha) = b*A(alpha,beta)*dg(beta) + c*A(alpha,3); 
            F1[0].add(dg1,A[0][0]); F1[0].add(dg2,A[0][1]); F1[0].times(b); F1[0].add(c,A[0][2]);
            F1[1].add(dg1,A[1][0]); F1[1].add(dg2,A[1][1]); F1[1].times(b); F1[1].add(c,A[1][2]);

            F2.add(dg1,A[2][0]);    F2.add(dg2,A[2][1]);    F2.times(b);    F2.add(c,A[2][2]);

            F1T[0].beTranspositionOf(F1[0]);
            F1T[1].beTranspositionOf(F1[1]);


            //f2 = c*A(3,3) + b*dg(alpha)*A(alpha,3)        
            F2f5.add(dg1,A[0][2]);    F2f5.add(dg2,A[1][2]);    F2f5.times(b);    F2f5.add(c,A[2][2]);
            F2T.beTranspositionOf(F2);

            //f3(alpha) = b*A(alpha,beta)*dm(beta) + zeta*A(alpha,3)*m; 
            t1.beProductOf(A[0][0],dm1); t2.beProductOf(A[0][1],dm2); t3.beProductOf(A[0][2],m);
            f3[0].add(b,t1); f3[0].add(b,t2); f3[0].add(zeta,t3);

            //f32 = b*( A.at(2,1)*dm1 + A.at(2,2)*dm2 )  + zeta*A.at(2,3)*m;
            t1.beProductOf(A[1][0],dm1); t2.beProductOf(A[1][1],dm2); t3.beProductOf(A[1][2],m);
            f3[1].add(b,t1); f3[1].add(b,t2); f3[1].add(zeta,t3);


            //f4 = b*dm(alpha)*A(alpha,3) + zeta*A(3,3)*m; 
            t1.beTProductOf(A[0][2],dm1); t2.beTProductOf(A[1][2],dm2); t3.beTProductOf(A[2][2],m);
            f4.add(b,t1); f4.add(b,t2);  f4.add(zeta,t3);

            // f5 = b*F1(alpha)*dm(alpha) + zeta*F2*m + zeta*N3
            t1.beProductOf(F1T[0],dm1); t2.beProductOf(F1T[1],dm2); t3.beProductOf(F2T,m);
            f5.add(b,t1); f5.add(b,t2); f5.add(zeta,t3); f5.add(zeta,S3g);


     /*
        
            A     a*A            F1           b*A*m             f3
           a*A    a^2*A          a*F1         a*b*A*m         a*f3+b*N
           F1     a*F1     b*F1*dgam+c*F2     b*(F1*m+N)        f5
           b*A*m  a*b*A*m    b*(F1*m+N)       b^2*m*A*m     b*m*f3
           f3     a*f3+b*N       f5           b*m*f3    b*dm*f3+zeta*m*f4
           */

            // L(1,1) = A(alpha,beta)
            FloatMatrix L11(6,6); L11.zero();
            L11.addSubMatrix(A[0][0],1,1); L11.addSubMatrix(A[0][1],1,4);
            L11.addSubMatrix(A[1][0],4,1); L11.addSubMatrix(A[1][1],4,4);

            // L(1,2) = a*A(alpha,beta)
            FloatMatrix L12(6,6); L12.zero(); L12.add(a,L11);
        
            // L(1,3) = F1(alpha)
            FloatMatrix L13(6,3); L13.zero();
            L13.addSubMatrix(F1[0],1,1); L13.addSubMatrix(F1[1],4,1);
        
            // L(1,4) = b*A*m
            FloatMatrix L14(6,2); L14.zero();
            t1.beProductOf(A[0][0],m); t2.beProductOf(A[0][1],m);
            L14.addSubVectorCol(t1,1,1); L14.addSubVectorCol(t2,1,2);
            t1.beProductOf(A[1][0],m); t2.beProductOf(A[1][1],m);
            L14.addSubVectorCol(t1,4,1); L14.addSubVectorCol(t2,4,2);
            L14.times(b);

            // L(1,5) = f3(alpha)
            FloatMatrix L15(6,1); L15.zero();
            L15.addSubVectorCol(f3[0],1,1); L15.addSubVectorCol(f3[1],4,1);


            // L(2,2) = a^2*A(alpha,beta)= a*L12
            FloatMatrix L22(6,6); L22.zero(); L22.add(a,L12);
        
            // L(2,3) = a*F1(alpha) = a*L13
            FloatMatrix L23(6,3); L23.zero(); L23.add(a,L13);
        
            // L(2,4) = a*F1(alpha)=a*L14
            FloatMatrix L24(6,2); L24.zero(); L24.add(a,L14);

            // L(2,5) = a*f3(alpha) + b*N(alpha)
            FloatMatrix L25(6,1); L25.zero(); 
            L25.addSubVectorCol(S1g,1,1); L25.addSubVectorCol(S2g,4,1); L25.times(b); L25.add(a,L15);
        

            // L(3,3) = b*F1(beta)*dgam(beta) + c*F2
            FloatMatrix L33(3,3); L33.zero();
            L33.add(c,F2T); temp3.add(dg1,F1T[0]); temp3.add(dg2,F1T[1]); temp3.times(b); L33.add(temp3);

                
            // L(3,4) = b*( F1(beta)*m + N(beta) )
            FloatMatrix L34(3,2); L34.zero(); 
            t1.beTProductOf(F1[0],m); t1.add(S1g); t1.times(b); 
            t2.beTProductOf(F1[1],m); t2.add(S2g); t2.times(b);
            L34.addSubVectorCol(t1,1,1); L34.addSubVectorCol(t2,1,2);

            // L(3,5) = f5
            FloatMatrix L35(3,1); L35.zero(); L35.addSubVectorCol(f5,1,1);


            // L(4,4) = b^2*m*A*m (2x2)
            FloatMatrix L44(2,2); L44.zero(); 
            t1.beProductOf(A[0][0],m); L44.at(1,1) = t1.dotProduct(m);
            t1.beProductOf(A[0][1],m); L44.at(1,2) = t1.dotProduct(m);
            t1.beProductOf(A[1][0],m); L44.at(2,1) = t1.dotProduct(m);
            t1.beProductOf(A[1][1],m); L44.at(2,2) = t1.dotProduct(m);
            L44.times(b*b);

            // L(4,5) = b*m*f3(alpha)
            FloatMatrix L45(2,1); 
            L45.at(1,1) = b*m.dotProduct(f3[0]);
            L45.at(2,1) = b*m.dotProduct(f3[1]);

            // L(5,5) = b*m*f3(alpha)
            FloatMatrix L55(1,1); 
            L55.at(1,1) = b*dm1.dotProduct(f3[0]) + b*dm2.dotProduct(f3[1]) + zeta*m.dotProduct(f4);
        
            FloatMatrix  L21, L31, L41, L43, L51, L32, L42, L52, L54, L53;
            L21.beTranspositionOf(L12); L31.beTranspositionOf(L13); L51.beTranspositionOf(L15); L32.beTranspositionOf(L23);
            L42.beTranspositionOf(L24); L43.beTranspositionOf(L34); L52.beTranspositionOf(L25); L54.beTranspositionOf(L45);
            L53.beTranspositionOf(L35); L41.beTranspositionOf(L14);

            // Assemble into L
            L.zero();
            L.addSubMatrix(L11,1,1);  
            L.addSubMatrix(L22,7,7);  
            L.addSubMatrix(L33,13,13); 
            L.addSubMatrix(L44,16,16);
            L.addSubMatrix(L55,18,18); 
        
            L.addSubMatrix(L12,1,7);   L.addSubMatrix(L21,7,1);       
            L.addSubMatrix(L13,1,13);  L.addSubMatrix(L31,13,1); 
            L.addSubMatrix(L14,1,16);  L.addSubMatrix(L41,16,1);  
            L.addSubMatrix(L15,1,18);  L.addSubMatrix(L51,18,1); 
            L.addSubMatrix(L23,7,13);  L.addSubMatrix(L32,13,7); 
            L.addSubMatrix(L24,7,16);  L.addSubMatrix(L42,16,7);
            L.addSubMatrix(L25,7,18);  L.addSubMatrix(L52,18,7);
            L.addSubMatrix(L34,13,16); L.addSubMatrix(L43,16,13);
            L.addSubMatrix(L35,13,18); L.addSubMatrix(L53,18,13);
            L.addSubMatrix(L45,16,18); L.addSubMatrix(L54,18,16); 
       
            //L.symmetrized();


            // Tangent matrix (K = B^T*L*B*dV)
            BtL.beTProductOf(B,L);
            BtLB.beProductOf(BtL,B);
            //dV = this->computeVolumeAround(gp);
            dV = this->computeVolumeAroundLayer(gp, layer);
            BtLB.times(dV);
            answer.add(BtLB);

        
        }
    }

    //answer.printYourself();
        //FloatMatrix test;
    //test.beSubMatrixOf(L,1,3,1,3);
    //test.beSubMatrixOf(L,4,6,4,6);
    //test.beSubMatrixOf(L,1,3,4,6);
    //test.beSubMatrixOf(L,1,6,7,12);
    //test.beSubMatrixOf(L,1,6,13,18);
    //test.beSubMatrixOf(L,1,6,16,18);
    //test.beSubMatrixOf(L,7,12,13,15);
    
    //test.beSubMatrixOf(L,7,12,16,18);
    
    //test.beSubMatrixOf(L,13,15,13,15);
    //test.beSubMatrixOf(L,13,15,16,18);
    //test.beSubMatrixOf(L,16,18,16,18);

    //test.beSubMatrixOf(L,1,6,1,6);
    //printf("\n new \n");
    //test.printYourself();
    //printf("\n Tangent \n");
    //test.printYourself();
   
    //testmatrix = L;

    }

void
TrDirShell :: computePressureTangentMatrix(FloatMatrix &answer, Load *load, const int iSurf, TimeStep *tStep)
{
    // Computes tangent matrix associated with the linearization of pressure loading. Assumes constant pressure.
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [1]; // rule #2 for mid-plane integration only

    double zeta, dA, a, b, c;
    FloatMatrix N, B, Nt, NtL, NtLB, L(7,18);
    FloatArray lcoords, BF;
    if( iSurf == 1 ){
        zeta = this->giveCrossSection()->give(CS_Thickness)*(-0.5); // bottom surface -> iSurf = 1
    }else if( iSurf == 2 ){         
        zeta = 0.0;                                                 // midplane surface -> iSurf = 2
    }else if( iSurf == 3 ){
        zeta = this->giveCrossSection()->give(CS_Thickness)*0.5;    // top surface -> iSurf = 3
    }else{
        _error("computePressureForce: incompatible load surface must be 1, 2 or 3");
    }   

    FloatArray solVec, pressure;
    this->giveUpdatedSolutionVector(solVec,tStep);
       


    answer.resize(42,42); answer.zero();
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        
        gp = iRule->getIntegrationPoint(i);
        lcoords = *gp->giveCoordinates();

        FloatArray g1, g2, g3, m(3), dm1(3), dm2(3), t1, t2;
        double gam, dg1, dg2;

        this->computeNmatrixAt(gp, N);
		this->computeBmatrixAt(gp, B);
        
        FloatArray genEps, temp1;	  
		genEps.beProductOf(B,solVec);  // [dxdxi, dmdxi, m, dgamdxi, gam]^T
        this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);
        this->evalCovarBaseVectorsAt(gp, g1, g2, g3, tStep, genEps); 
        FloatArray n; n.beVectorProductOf(g1,g2); // Compute normal (should not be normalized)
        
        // thickness coefficients
        a = zeta + 0.5*gam*zeta*zeta;  b = 0.5*zeta*zeta;  c = 1.0 + gam*zeta;       
        
        FloatMatrix W1(3,3), W2(3,3), temp3(3,3);

        // W = skew-symmetric matrix 
        W2.resize(3,3); W2.zero(); 
        W2.at(2,3) = -g2.at(1); W2.at(3,2) =  g2.at(1);
        W2.at(1,3) =  g2.at(2); W2.at(3,1) = -g2.at(2);
        W2.at(1,2) = -g2.at(3); W2.at(2,1) =  g2.at(3);
        W2.negated();

        W1.resize(3,3); W1.zero();
        W1.at(2,3) = -g1.at(1); W1.at(3,2) =  g1.at(1);
        W1.at(1,3) =  g1.at(2); W1.at(3,1) = -g1.at(2);
        W1.at(1,2) = -g1.at(3); W1.at(2,1) =  g1.at(3);
        W1.negated();



        // Tangent stiffness

       /*
       [    W2     - W1      a*W2      -a*W1       b(W2*dg1-W1*dg2)          b*W2*m1      -b*W1*m2,      b(W2*dm1 - W1*dm2)       (3x18)
          a*W2    -a*W1    a^2*W2      -a*W1      ab(W2*dg1-W1*dg2)        a*b*W2*m1    -a*b*W1*m2,     ab(W2*dm1 - W1*dm2)+bn    (3x18)
        m^t*W2  -m^t*W1  a*m^t*W2  -a*m^t*W1   m^t*b(W2*dg1-W1*dg2)+n^t  b*m^t*W2*m1  -b*m^t*W1*m2  b*m^t*(W2*dm1 - W1*dm2)    ]  (1x18)
       */
        // L(1,1) = [W2, -W1] 
        FloatMatrix L11(3,6); L11.zero();
        temp3.add(-1,W1); 
        L11.addSubMatrix(W2,1,1); L11.addSubMatrix(temp3,1,4); 

        // L(1,2) = a*[W2, -W1] = a*L11
        FloatMatrix L12(3,6); L12.zero(); L12.add(a,L11);

        // L(1,3) =  b(W2*dg1-W1*dg2)
        FloatMatrix L13(3,3); L13.zero();
        L13.add(dg1,W2); L13.add(-dg2,W1); L13.times(b);

        // L(1,4) = [b*W2*m, -b*W1*m]
        FloatMatrix L14(3,2); L14.zero();
        t1.beProductOf(W2,m); t2.beProductOf(W1,m); t2.negated();
        L14.addSubVectorCol(t1,1,1); L14.addSubVectorCol(t2,1,2); 
        L14.times(b);

        // L(1,5) = b(W2*dm1 - W1*dm2)
        FloatArray L15(3); L15.zero();
        t1.beProductOf(W2,dm1); t2.beProductOf(W1,dm2);
        L15.add(b,t1); L15.add(-b,t2); //L15.times(b);
        
        FloatMatrix L1(3,18); L1.zero();
        L1.addSubMatrix(L11,1,1); L1.addSubMatrix(L12,1,7);
        L1.addSubMatrix(L13,1,13); L1.addSubMatrix(L14,1,16);
        L1.addSubVectorCol(L15,1,18);


        FloatMatrix L2(3,18); L2.zero();
        L2.add(a,L1); t1.zero(); t1.add(b,n);
        L2.addSubVectorCol(t1,1,18);


        // L(3,:) = m^t*L(:,1) + n^t*deltam 
        FloatArray L3; L3.beTProductOf(L1,m); L3.times(b);
        L3.at(13) += b*n.at(1); L3.at(14) += b*n.at(2); L3.at(15) += b*n.at(3);

        L.zero(); L.addSubMatrix(L1,1,1); L.addSubMatrix(L2,4,1); L.addSubVectorRow(L3,7,1);

        load->computeValueAt(pressure, tStep, * ( gp->giveCoordinates() ), VM_Total); // pressure component 
        L.times(-pressure.at(1));


        // Tangent matrix (K = N^T*L*B*dA)
        NtL.beTProductOf(N,L);
        NtLB.beProductOf(NtL,B);
        dA = this->computeAreaAround(gp);
        NtLB.times(dA);
        answer.add(NtLB);

        
    }
    
    FloatMatrix test, compmat;


    //printf("\n K(7,7,1:18) \n");
    //test.beSubMatrixOf(L,4,6,13,18);
    //test.printYourself();
    
       /*        
       printf("\n K(1:6,1:6) \n");
       test.beSubMatrixOf(answer,1,6,1,6);
       test.printYourself();

       printf("\n K(7:12,1:6) \n");
       test.beSubMatrixOf(answer,7,12,1,6);
       test.printYourself();

       printf("\n K(1:6,7:12) \n");
       test.beSubMatrixOf(answer,1,6,7,12);
       test.printYourself();
       */


    /*
       //compareMatrices(answer, testmatrix, compmat);
       compmat.resize(42,42); compmat.zero();
       compmat.add(testmatrix);
       compmat.subtract(answer);

       printf("\n compmat(1:6,1:6) \n");
       test.beSubMatrixOf(compmat,37,42,1,6);
       test.printYourself();

       printf("\n compmat(1:6,:12) \n");
       test.beSubMatrixOf(compmat,37,42,7,12);
       test.printYourself();

       printf("\n compmat(1:6,13:18) \n");
       test.beSubMatrixOf(compmat,37,42,13,18);
       test.printYourself();

       printf("\n compmat(1:6,19:24) \n");
       test.beSubMatrixOf(compmat,37,42,19,24);
       test.printYourself();

       printf("\n compmat(1:6,25:30) \n");
       test.beSubMatrixOf(compmat,37,42,25,30);
       test.printYourself();

       printf("\n compmat(1:6,31:36) \n");
       test.beSubMatrixOf(compmat,37,42,31,36);
       test.printYourself();

       printf("\n compmat(1:6,37:42) \n");
       test.beSubMatrixOf(compmat,37,42,37,42);
       test.printYourself();
       */
}



// Strain and stress

void 
TrDirShell :: computeFAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *stepN, FloatArray &genEps){
	// Compute deformation gradient as open product(g_i, G_i)
	FloatArray gcov1, gcov2, gcov3, Gcon1, Gcon2, Gcon3;
	this->evalCovarBaseVectorsAt(gp, gcov1, gcov2, gcov3, stepN, genEps);
	this->evalInitialContravarBaseVectorsAt(gp, Gcon1, Gcon2, Gcon3);

	FloatMatrix F1, F2, F3;
	F1.beDyadicProductOf(gcov1,Gcon1);
	F2.beDyadicProductOf(gcov2,Gcon2);
	F3.beDyadicProductOf(gcov3,Gcon3);

    answer.resize(3,3); answer.zero();
	answer.add(F1);	answer.add(F2);	answer.add(F3);
   // Gcon3.printYourself();
}

void
TrDirShell :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &genEps){
    // Computes the Green-Lagrange strain tensor: E=0.5(C-I)
	FloatMatrix F, E;
	this->computeFAt(gp, F, stepN, genEps); // Deformation gradient

    E.beTProductOf(F,F); // C-Right Caucy-Green deformation tensor
    E.at(1,1) += -1;
	E.at(2,2) += -1;
	E.at(3,3) += -1;
	E.times(0.5);


    FloatArray temp(6);
    temp.beReducedVectorForm(E); // Convert to Voight form Todo: add enum strain/stress
    answer=temp;
    answer.at(4) = temp.at(4)*2.0;
    answer.at(5) = temp.at(5)*2.0;
    answer.at(6) = temp.at(6)*2.0;
    //answer.printYourself();
}

void
TrDirShell :: computeStressVector(FloatArray &answer, FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *stepN )
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step stepN. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{

    FloatArray Epsilon;
    StructuralCrossSection *cs = ( StructuralCrossSection * ) this->giveCrossSection();

    this->computeStrainVector(Epsilon, gp, stepN, genEps);

    //( ( StructuralMaterial * ) gp->giveElement()->giveMaterial() )
    //->giveRealStressVector(answer, ReducedForm, gp, Epsilon, stepN);

    ( ( StructuralMaterial * ) mat )->giveRealStressVector(answer, ReducedForm, gp, Epsilon, stepN);
    

}

void 
TrDirShell :: computeStressResultantsAt(GaussPoint *gp, FloatArray &Svec, FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, TimeStep *tStep, FloatArray &genEps){

     FloatArray g1, g2, g3;
     this->evalCovarBaseVectorsAt(gp, g1, g2, g3, tStep, genEps);

     FloatMatrix S; 
	 S.beMatrixForm(Svec);
     
     // Sig =S(i,j)*g_j, - stress vectors on the surfaces given by g_j?
     for( int j=1; j<=3; j++){
        S1g.at(j) = S.at(1,1)*g1.at(j) + S.at(1,2)*g2.at(j) + S.at(1,3)*g3.at(j);
        S2g.at(j) = S.at(2,1)*g1.at(j) + S.at(2,2)*g2.at(j) + S.at(2,3)*g3.at(j);
        S3g.at(j) = S.at(3,1)*g1.at(j) + S.at(3,2)*g2.at(j) + S.at(3,3)*g3.at(j);
     }

 }



// Base vectors and change of bases

void
TrDirShell :: transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatArray &VoightMatrix, FloatArray &answer){
    // Transform stress in cart system to curvilinear sytem
    // New: uses Bond transformation matrix. No need to go from matrix to Voigt form and back.
    FloatArray Gcon1, Gcon2, Gcon3;
    this->evalInitialContravarBaseVectorsAt(gp, Gcon1, Gcon2, Gcon3);
	FloatMatrix GE, M;
    this->giveCoordTransMatrix(GE, Gcon1, Gcon2, Gcon3);  
    this->giveBondTransMatrix(M, GE);
    answer.beProductOf(M,VoightMatrix);

    /*
    FloatMatrix EG(3,3);
    EG.at(1,1) = Gcon1.at(1); EG.at(1,2) = Gcon2.at(1); EG.at(1,3) = Gcon3.at(1);
    EG.at(2,1) = Gcon1.at(2); EG.at(2,2) = Gcon2.at(2); EG.at(2,3) = Gcon3.at(2);
	EG.at(3,1) = Gcon1.at(3); EG.at(3,2) = Gcon2.at(3); EG.at(3,3) = Gcon3.at(3);
	

	
	// transform according to: EG^T*mat*EG
    FloatMatrix temp(3,3), temp2(3,3), mat(3,3);
	mat.beMatrixForm(VoightMatrix);

    temp.beTProductOf(EG,mat);
    temp2.beProductOf(temp,EG);
	answer.beReducedVectorForm(temp2);
    */

}

void
TrDirShell :: transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatMatrix &stiffness, FloatMatrix &answer){
    // Transform stifness tensor in cart system to curvilinear sytem
    // New: uses Bond transformation matrix. No need to go from matrix to Voigt form and back.
    FloatArray Gcon1, Gcon2, Gcon3;
    this->evalInitialContravarBaseVectorsAt(gp, Gcon1, Gcon2, Gcon3);
	FloatMatrix GE, M, temp;
    this->giveCoordTransMatrix(GE, Gcon1, Gcon2, Gcon3);  
    this->giveBondTransMatrix(M, GE);

    temp.beProductTOf(stiffness,M);
    answer.beProductOf(M,temp);


}





void
TrDirShell :: giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep){
	
	FloatArray *Xi, Mi;
	this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, answer);

	// Add reference position and directors to computed update
	for( int i = 1, j = 0; i <= 6; i++, j += 7 ){
		Xi = this->giveNode(i)->giveCoordinates();
		Mi = this->giveInitialNodeDirector(i);
		answer.at(1+j) += Xi->at(1);
		answer.at(2+j) += Xi->at(2);
		answer.at(3+j) += Xi->at(3);
		answer.at(4+j) += Mi.at(1);
		answer.at(5+j) += Mi.at(2);
		answer.at(6+j) += Mi.at(3);
		//answer.at(7+j) =... gam(t=0)=0 so no update necessary. Well this assumes gam=0 at t=0
        //Mi.printYourself();
	}

}


void
TrDirShell :: edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep)
{	
	FloatArray *Xi, Mi;
	// should extract solution vector from three nodes only
    //this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, answer);
    this->computeBoundaryVectorOf(iedge, EID_MomentumBalance, VM_Total, tStep, answer);
    IntArray edgeNodes;
    this->interpolation.computeLocalEdgeMapping(edgeNodes, iedge);


	// Add reference position and directors to computed update
	for( int i = 1, j = 0; i <= 3; i++, j += 7 ){
		Xi = this->giveNode( edgeNodes.at(i) )->giveCoordinates();
		Mi = this->giveInitialNodeDirector( edgeNodes.at(i) );
		answer.at(1+j) += Xi->at(1);
		answer.at(2+j) += Xi->at(2);
		answer.at(3+j) += Xi->at(3);
		answer.at(4+j) += Mi.at(1);
		answer.at(5+j) += Mi.at(2);
		answer.at(6+j) += Mi.at(3);
		//answer.at(7+j) =... gam(t=0)=0 so no update necessary
	}

}



// Nodal forces
void
TrDirShell :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//
// Computes internal forces as a summation of: sectional forces + convective mass force + pressure load
{
    FloatArray fM, fP, fT;
    
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec,tStep);
    //solVec.printYourself();
    this->computeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord);
    //answer.printYourself();

    //this->computeConvectiveMassForce(fM, tStep);
    //answer.subtract(fM);

//    this->computePressureForce(fP, tStep);
//    answer.add(fP);


    // -----------------------------------
    // Numerical tangent

    /* 
    FloatMatrix Knum(42,42); Knum.zero(); 
    FloatArray sectForce_h(42);
    double h=1.0e-8;
    for( int i=1; i<=42; i++ ){
        //solVec.zero();
        this->giveUpdatedSolutionVector(solVec,tStep);
        solVec.at(i) += h;
        this->computeSectionalForces(sectForce_h, tStep, solVec, useUpdatedGpRecord);
       
        Knum.addSubVectorCol(sectForce_h,1,i);


    }
    Knum.times(1.0/h);
    
    FloatMatrix test(3,3);
    test.beSubMatrixOf(Knum,1,3,1,3);
    //test.beSubMatrixOf(Knum,4,6,4,6);
    //test.beSubMatrixOf(Knum,1,3,4,6);
    printf("\n Num tangent \n");
    test.printYourself();
    */
    
    
}

void
TrDirShell :: computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord)
//
{
    GaussPoint *gp;
    //Material *mat = this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [0];


    double dV;
    FloatMatrix B, Bt;
    FloatArray lcoords, cartStressVector, contravarStressVector, sectionalForces, BF, f;
    FloatArray genEps;	      // generalized strain

#if 0
    answer.resize(42); answer.zero();

	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
                
		this->computeBmatrixAt(gp, B);
		genEps.beProductOf(B,solVec); // [dxdxi, dmdxi, m, dgamdxi, gam]^T
        double zeta = giveLocalZetaCoord(gp);
        this->computeSectionalForcesAt(f, gp, tStep, genEps, zeta);

        BF.beTProductOf(B,f);
        dV = this->computeVolumeAround(gp);
        BF.times(dV);
        answer.add(BF);

    }
    //answer.printYourself();
#endif       

    answer.resize(42); answer.zero();
    
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >(this->giveCrossSection());
    int numberOfLayers = layeredCS->give(CS_NumLayers);

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = layerIntegrationRulesArray [layer-1]; 
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

	    for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {    
            gp = iRuleL->getIntegrationPoint(j-1);      
            dV = this->computeVolumeAroundLayer(gp, layer);

            this->computeBmatrixAt(gp, B);
		    genEps.beProductOf(B,solVec); 
            double zeta = giveLayerZetaCoord(gp,layer);
            this->computeSectionalForcesAt(f, gp, mat, tStep, genEps, zeta);
            BF.beTProductOf(B,f);

            BF.times(dV);
            answer.add(BF);

        }
    }

    //answer.printYourself();

   /*
    answer.resize(42); answer.zero();
    IntegrationRule *iRuleL = integrationRulesArray [1]; 
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >(this->giveCrossSection());
    GaussPoint *mastergp, *layergp, *sublayergp;
    int numberOfLayers = layeredCS->give(CS_NumLayers);

	for ( int i = 0; i < iRuleL->getNumberOfIntegrationPoints(); i++ ) {    // for each master integration point
        mastergp = iRuleL->getIntegrationPoint(i);
        
        dA = this->computeAreaAround(mastergp);
        for ( int j = 0; j < numberOfLayers; j++ ) {

            // integration
            layergp = layeredCS->giveSlaveGaussPointNew(mastergp, j);
            for( int k = 0; k < layeredCS->giveNumIntegrationPointsInLayer(); k++ ){
                sublayergp = layeredCS->giveSlaveGaussPointNew(layergp, k);
            }
            
            //dz = this->computeThicknessAroundLayer(layergp, j);
            dz = layeredCS->giveLayerThickness(j) * layergp->giveWeight();


            this->computeBmatrixAt(layergp, B);
		    genEps.beProductOf(B,solVec); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

            this->computeSectionalForcesAt(f, layergp, tStep, genEps);
            BF.beTProductOf(B,f);

            BF.times(dz*dA);
            answer.add(BF);

        }
    }

    //answer.printYourself();

*/





    // -----------------------------------
    // Numerical tangent

    /*
    FloatMatrix Knum(18,18); Knum.zero(); 
    FloatArray f_hp(18), f_hm(18), genEps_hp, genEps_hm, f0;
    f0 = f; f0.negated();
    double h=1.0e-8;
    for( int i=1; i<=18; i++ ){
        genEps_hp = genEps; genEps_hm = genEps;
        genEps_hp.at(i) += h; genEps_hm.at(i) -= h;
        this->computeSectionalForcesAt(f_hp, gp, tStep, genEps_hp);
        this->computeSectionalForcesAt(f_hm, gp, tStep, genEps_hm);
        f_hm.negated();
        Knum.addSubVectorCol(f_hm,1,i);
        Knum.addSubVectorCol(f_hp,1,i);


    }
    Knum.times(1.0/(2.0*h));
    
    FloatMatrix test, diff(18,18), temp(18,18), KnumT(18,18);
    //test.beSubMatrixOf(Knum,1,3,1,3);
    //test.beSubMatrixOf(Knum,4,6,4,6);
    //test.beSubMatrixOf(Knum,1,3,4,6);
    //test.beSubMatrixOf(Knum,1,6,7,12);
    test.beSubMatrixOf(Knum,1,6,13,18);
    //test.beSubMatrixOf(Knum,1,6,16,18);
    //test.beSubMatrixOf(Knum,7,12,13,15);
    //test.beSubMatrixOf(Knum,13,15,13,15);
    //test.beSubMatrixOf(Knum,13,15,16,18);
    //test.beSubMatrixOf(Knum,16,18,16,18);
    
    //test.beSubMatrixOf(Knum,7,12,16,18);
    
    test.beSubMatrixOf(Knum,13,18,13,18);
    
    
    //printf("\n Num tangent \n");
    //test.printYourself();

    
    //diff=testmatrix;
    //diff.subtract(Knum);

    this->compareMatrices( testmatrix, Knum, diff);
    
    test.beSubMatrixOf(diff,1,6,1,6);
    printf("\n diff 1-6,1-6 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,1,6,7,12);
    printf("\n diff 1-6,7-12 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,1,6,13,18);
    printf("\n diff 1-6,13-18 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,7,12,1,6);
    printf("\n diff 7-12,1-6 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,7,12,7,12);
    printf("\n diff 7-12,7-12 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,7,12,13,18);
    printf("\n diff 7-12,13-18 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,13,18,1,6);
    printf("\n diff 13-18,1-6 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,13,18,7,12);
    printf("\n diff 13-18,7-12 \n");
    test.printYourself();

    test.beSubMatrixOf(diff,13,18,13,18);
    printf("\n diff 13-18,13-18 \n");
    test.printYourself();
   */
}

void
TrDirShell :: computeSectionalForcesAt(FloatArray &answer, GaussPoint *gp, Material *mat, TimeStep *tStep, FloatArray &genEps, double zeta){

        FloatArray f(18), g1, g2, g3, S1g(3), S2g(3), S3g(3), m(3), dm1(3), dm2(3), temp1, temp2;
        FloatArray lcoords, cartStressVector, contravarStressVector, sectionalForces, BF, a;
        double fac1, fac2, fac3, gam, dg1, dg2;

        lcoords = *gp->giveCoordinates();
        //zeta = giveLocalZetaCoord(gp);

        this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);
        
        this->computeStressVector(cartStressVector, genEps, gp, mat, tStep);
       
        this->transInitialCartesianToInitialContravar(gp, cartStressVector, contravarStressVector);
        this->computeStressResultantsAt(gp, contravarStressVector, S1g, S2g, S3g, tStep, genEps);

        fac1 = zeta + 0.5*gam*zeta*zeta;
        fac2 = 1. + gam*zeta;
        fac3 = 0.5*zeta*zeta;

        /* The following expressions are sectional forces if integrated over the thickness. However, 
           now they are integrated over the volume to produce f_int.*/
 
        // Membrane forces - N1, N2 
        f.at(1) = S1g.at(1); f.at(2) = S1g.at(2); f.at(3) = S1g.at(3);
        f.at(4) = S2g.at(1); f.at(5) = S2g.at(2); f.at(6) = S2g.at(3);

        // Moments - M1, M2
        f.at(7)  = fac1 * S1g.at(1);
        f.at(8)  = fac1 * S1g.at(2);
        f.at(9)  = fac1 * S1g.at(3);
        f.at(10) = fac1 * S2g.at(1);
        f.at(11) = fac1 * S2g.at(2);
        f.at(12) = fac1 * S2g.at(3);

        // Shear force - T
        f.at(13) = fac2 * S3g.at(1) + fac3 * ( S1g.at(1)*dg1 + S2g.at(1)*dg2) ;
        f.at(14) = fac2 * S3g.at(2) + fac3 * ( S1g.at(2)*dg1 + S2g.at(2)*dg2) ;
        f.at(15) = fac2 * S3g.at(3) + fac3 * ( S1g.at(3)*dg1 + S2g.at(3)*dg2) ;

        // Higher order force - Ms
        f.at(16) = fac3 * m.dotProduct(S1g);
        f.at(17) = fac3 * m.dotProduct(S2g);

        // Ts
        f.at(18) = fac3 * ( dm1.dotProduct(S1g) + dm2.dotProduct(S2g) ) + zeta * m.dotProduct(S3g) ;
        answer=f;
}

void 
TrDirShell :: computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep){
	// Analytically integrated over the thickness. Constant density assumed.
    
    IntegrationRule *iRule = integrationRulesArray [1]; // rule 2 for mid-plane integration only
	GaussPoint *gp;


    Material *mat = this->giveMaterial();
    FloatMatrix N, Nt;
    FloatArray lcoords, cartStressVector, contravarStressVector, sectionalForces, BF, a, da, unknowns, m, dm, aVec, daVec, fm(7), fM;
    double gam, dgam, dA;
    
    answer.resize(42); answer.zero();
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {

        gp = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(gp, N);
        this->giveUpdatedSolutionVector(aVec,tStep);
        this->computeVectorOf(EID_MomentumBalance, VM_Velocity, tStep, daVec);
	       
         a.beProductOf(N,aVec);  // [ x,  m,  gam]^T
        da.beProductOf(N,daVec); // [dx, dm, dgam]^T
	     m.setValues(3,  a.at(4),  a.at(5),  a.at(6) );  gam =  a.at(7);
        dm.setValues(3, da.at(4), da.at(5), da.at(6) ); dgam = da.at(7);

        double a1, a2, a3, h, h2, h3, h5, fac1, fac2, fac3, rho; 
        FloatArray coeff;


        rho = this->giveMaterial()->give('d', gp);
        h   = this->giveCrossSection()->give(CS_Thickness);
        h2  = h*h; h3 = h2*h; h5 = h2*h3; 
        this->computeThicknessMappingCoeff(gp, coeff);
        a1 = coeff.at(1); a2 = coeff.at(2); a3 = coeff.at(3);

        // Convective mass "matrix"
	    /*     3
  	       3 [ a*m
	       3   b*m
	       1  c*m.dm]*dgam
	    */
        
        fac1 = rho*h3*(20.*a3 + 3.*a1*h2)/240.;
        fac2 = rho*h5*(56.*a2 + 28.*a3*gam + 5.*a1*h2*gam)/4480.;
        fac3 = rho*h5*(28.*a3 + 5.*a1*h2)/4480.;
        fm.at(1) = fac1 * dm.at(1)*dgam;
        fm.at(2) = fac1 * dm.at(2)*dgam;
        fm.at(3) = fac1 * dm.at(3)*dgam;
        fm.at(4) = fac2 * dm.at(1)*dgam;
        fm.at(5) = fac2 * dm.at(2)*dgam;
        fm.at(6) = fac2 * dm.at(3)*dgam;
        fm.at(7) = fac3 * m.dotProduct(dm)*gam;

        dA = this->computeAreaAround(gp); 
        fM.beTProductOf(N,fm);
        fM.times(dA);
        answer.add(fM);
    }

}






void 
TrDirShell :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep){
	FloatMatrix mConsistent;

	//this->computeMassMatrix(mConsistent, tStep);
	// Reduce to lumped form 
	// Todo: add algorithm for this
	//answer.resize(mConsistent);
	//answer.resize(42,42); answer.zero();
	//for( int i = 1; i<=42; i++){
	//	answer.at(i,i) = mConsistent.at(i,i);
	//}
    OOFEM_ERROR("TrDirShell :: computeLumpedMassMatrix - No lumping algorithm implemented");
}


void 
TrDirShell :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep){


    //this->computeMassMatrixNum(answer, tStep);
    //return;
    
    // Analytically integrated over the thickness. Constant density assumed.
    // => integration rule #2
	IntegrationRule *iRule = integrationRulesArray [ 1 ];
	GaussPoint *gp;

	//------------------------------
    FloatMatrix N, Nt, Ntm, NtmN, mass;
    FloatArray a, unknowns, m(3);
    double gam, dA;
    this->giveUpdatedSolutionVector(a, tStep);  

	answer.resize(42,42); answer.zero();
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);

		this->computeNmatrixAt(gp, N);
		unknowns.beProductOf(N,a); // [x, m, gam]^T
		m.setValues(3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
		gam = unknowns.at(7);
		
		// Consistent mass matrix
		/*     3    3    1
  		   3 [a*I  b*I   c*m      [A  B  C
		   3       d*I   e*m    =     D  E
		   1  sym       f*m.m]     sym   F]
		*/
		double a1, a2, a3;
		FloatArray coeff;
		this->computeThicknessMappingCoeff(gp, coeff);
		a1 = coeff.at(1); a2 = coeff.at(2); a3 = coeff.at(3);


		double h, h2, h3, h5, fac1, fac2, fac3, fac4, fac5, fac6, gam2, rho; 
		rho = this->giveMaterial()->give('d', gp);
		h   = this->giveCrossSection()->give(CS_Thickness);
		h2  = h*h; h3 = h2*h; h5 = h2*h3;
		gam2 = gam*gam;

		mass.resize(7,7);
		fac1 = a3*h + (a1*h3)/12.;
		fac2 = h3*(40.*a2 + 20.*a3*gam + 3.*a1*h2*gam)/480. ;
		fac3 = h3*(20.*a3 + 3.*a1*h2)/480. *1.0;
		fac4 = (28.*a3*h3*(80. + 3.*h2*gam2) + 3.*h5*(112.*a2*gam + a1*(112. + 5.*h2*gam2)))/26880.;
		fac5 = h*(56.*a2 + 28.*a3*gam + 5.*a1*h2*gam)/8960. ;
		fac6 = h5*(28.*a3 + 5.*a1*h2)/8960.;
		mass.at(1,1) = mass.at(2,2) = mass.at(3,3) = fac1; // A
		mass.at(1,4) = mass.at(2,5) = mass.at(3,6) = fac2; // B
		mass.at(1,7) = fac3 * m.at(1);  
		mass.at(2,7) = fac3 * m.at(2); 
		mass.at(3,7) = fac3 * m.at(3); // C
		mass.at(4,4) = mass.at(5,5) = mass.at(6,6) = fac4; // D
		mass.at(4,7) = fac5 * m.at(1);  
		mass.at(5,7) = fac5 * m.at(2); 
		mass.at(6,7) = fac5 * m.at(3); // E
		mass.at(7,7) = fac6 * m.dotProduct(m); // F
        mass.symmetrized();

        dA = this->computeAreaAround(gp);
        Ntm.beTProductOf(N, mass);
		NtmN.beProductOf(Ntm, N);
        NtmN.times(dA*rho);
		answer.add(NtmN);
        
	}


}






// Edge 
void 
TrDirShell :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode)
{
    FloatArray fT, components, traction(3);

    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >(load);
    
    
    if ( edgeLoad ) {
      
        this->computeTractionForce(fT, iEdge, edgeLoad, tStep);

        IntArray mask;
        this->giveEdgeDofMapping(mask, iEdge);
        answer.resize( this->computeNumberOfDofs(EID_MomentumBalance) );
        answer.zero();

        answer.assemble(fT, mask);

        return;
    } else {
        _error("computeEdgeLoadVectorAt: incompatible load");
        return;
    }
}

void 
TrDirShell :: computeTractionForce(FloatArray &answer, const int iedge, BoundaryLoad *edgeLoad, TimeStep *tStep)
{
	
    IntegrationRule *iRule = integrationRulesArray [2]; // rule #3 for edge integration of distributed loads given in [*/m]
	GaussPoint *gp;

    FloatMatrix N, Q;
    FloatArray g1, g2, g3, FT, fT(7), m, aVec, a, T(3), M, G1, components, lcoords;
    double dA;
    answer.resize(21); answer.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        
        gp = iRule->getIntegrationPoint(i);
        lcoords = *gp->giveCoordinates();
    
        edgeLoad->computeValueAt(components,tStep, lcoords, VM_Total);
        this->edgeComputeNmatrixAt(gp, N);
        this->edgeEvalCovarBaseVectorsAt(gp, iedge, g2, g3, tStep);
        g2.normalize();
        g3.normalize();
        g1.beVectorProductOf(g2,g3);  g1.normalize();
        this->giveCoordTransMatrix(Q, g1, g2, g3);
             
        FloatArray c1(3), c2(3), t1, t2; 
        c1.at(1) = components.at(1); c1.at(2) = components.at(2); c1.at(3) = components.at(3);
        c2.at(1) = components.at(4); c2.at(2) = components.at(5); c2.at(3) = components.at(6);
        t1.beTProductOf(Q,c1);
        t2.beTProductOf(Q,c2);
        
        fT.at(1) = t1.at(1); fT.at(2) = t1.at(2); fT.at(3) = t1.at(3);
        fT.at(4) = t2.at(1); fT.at(5) = t2.at(2); fT.at(6) = t2.at(3);
        fT.at(7) = components.at(7);
        fT = components;
        dA = this->edgeComputeAreaAround(gp, iedge) ; 
        FT.beTProductOf(N,fT);
        FT.times(dA);
        answer.add(FT);

    }
}

void
TrDirShell :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */
    //IntArray edgeNodes;
    //this->interpolation.computeLocalEdgeMapping(edgeNodes, iEdge);
    //answer.resize(4);

    
    if ( iEdge == 1 ) { // edge between nodes 1-4-2
        answer.setValues(21, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,  22, 23, 24, 25, 26, 27, 28);

    } else if ( iEdge == 2 ) { // edge between nodes 2-5-3
        answer.setValues(21, 8, 9, 10, 11, 12, 13, 14,  15, 16, 17, 18, 19, 20, 21, 29, 30, 31, 32, 33, 34, 35 );

    } else if ( iEdge == 3 ) { // edge between nodes 3-6-1
        answer.setValues(21, 15, 16, 17, 18, 19, 20, 21,  1, 2, 3, 4, 5, 6, 7,  36, 37, 38, 39, 40, 41, 42);

    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
    

}




// Surface
void
TrDirShell :: computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                                int iSurf, TimeStep *tStep, ValueModeType mode)
{
    FloatMatrix T;
    FloatArray force, solVec;
    IntArray mask;

    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( surfLoad ) {
        this->giveUpdatedSolutionVector(solVec,tStep);
        this->computePressureForce(force, solVec, iSurf, surfLoad, tStep);
        
        this->giveSurfaceDofMapping(mask, 1);  // same dofs regardless of iSurf
        answer.resize( this->computeNumberOfDofs(EID_MomentumBalance) );
        answer.zero();
        answer.assemble(force, mask);


        /*
       // Numerical tangent 
            FloatMatrix Knum(42,42); Knum.zero(); 
            FloatArray f_hp(42), f_hm(42), solVec_hp, solVec_hm, f0;
    
            double h=1.0e-8;
            for( int i=1; i<=42; i++ ){
                solVec_hp = solVec; solVec_hm = solVec;
                solVec_hp.at(i) += h; solVec_hm.at(i) -= h;
                this->computePressureForce(f_hp, solVec_hp, iSurf, surfLoad, tStep);
                this->computePressureForce(f_hm, solVec_hm, iSurf, surfLoad, tStep);
                f_hm.negated();
                Knum.addSubVectorCol(f_hm,1,i);
                Knum.addSubVectorCol(f_hp,1,i);

            }

            
            Knum.times(1.0/(2.0*h));
            testmatrix = Knum;
    */
            /*
            FloatMatrix test;         
            printf("\n Knum(1:6,1:6) \n");
            test.beSubMatrixOf(Knum,1,6,1,6);
            test.printYourself();

            printf("\n Knum(7:12,1:6) \n");
            test.beSubMatrixOf(Knum,7,12,1,6);
            test.printYourself();

            printf("\n Knum(1:6,7:12) \n");
            test.beSubMatrixOf(Knum,1,6,7,12);
            test.printYourself();
            */

        return;
    } else {
        _error("computeSurfaceLoadVectorAt: incompatible load");
        return;
    }
}

void 
TrDirShell :: computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep){
	// Computes pressure loading. Acts normal to the current (deformed) surface.

    // Should be special integration rule for top and bottom surface!!
    IntegrationRule *iRule = integrationRulesArray [1]; // rule #2 for mid-plane integration only
	GaussPoint *gp;

    FloatMatrix N, B;
    FloatArray Fp, fp, temp, n(3), m, aVec, a, load, traction, genEps;;
    //this->giveUpdatedSolutionVector(solVec,tStep);

    
    answer.resize(42); answer.zero();
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        
        gp = iRule->getIntegrationPoint(i);
        this->computeBmatrixAt(gp, B);
        this->computeNmatrixAt(gp, N);
        genEps.beProductOf(B,solVec); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

        this->computePressureForceAt(gp, fp, iSurf, genEps, surfLoad, tStep);

        double dA = this->computeAreaAround(gp) ; 
        Fp.beTProductOf(N,fp);
        Fp.times(dA);
        answer.add(Fp);

    }


    /*   
    // -----------------------------------
    // Numerical tangent

    
    FloatMatrix Knum(7,18); Knum.zero(); 
    FloatArray f_hp(18), f_hm(18), genEps_hp, genEps_hm, f0;
    
    double h=1.0e-4;
    for( int i=1; i<=18; i++ ){
        genEps_hp = genEps; genEps_hm = genEps;
        genEps_hp.at(i) += h; genEps_hm.at(i) -= h;
        this->computePressureForceAt(gp, f_hp, iSurf, genEps_hp, surfLoad, tStep);
        this->computePressureForceAt(gp, f_hm, iSurf, genEps_hm, surfLoad, tStep);
        f_hm.negated();
        Knum.addSubVectorCol(f_hm,1,i);
        Knum.addSubVectorCol(f_hp,1,i);

    }


    Knum.times(1.0/(2.0*h));
    FloatMatrix test;
    
    

    
    printf("\n Knum() \n");
    test.beSubMatrixOf(Knum,4,6,13,18);
    test.printYourself();
    */
    
    //printf("\n Knum(7,1:6) \n");
    //test.beSubMatrixOf(Knum,7,7,1,6);
    //test.printYourself();
    
}

void 
TrDirShell :: computePressureForceAt(GaussPoint *gp, FloatArray &answer, const int iSurf, FloatArray genEps, BoundaryLoad *surfLoad, TimeStep *tStep){
	// Computes pressure loading. Acts normal to the current (deformed) surface.


    FloatArray g1, g2, g3, m, load, traction;
    double gam, zeta;
    answer.resize(7); answer.zero();

    if( iSurf == 1 ){
        zeta = this->giveCrossSection()->give(CS_Thickness)*(-0.5); // bottom surface -> iSurf = 1
    }else if( iSurf == 2 ){         
        zeta = 0.0;                                                 // midplane surface -> iSurf = 2
    }else if( iSurf == 3 ){
        zeta = this->giveCrossSection()->give(CS_Thickness)*0.5;    // top surface -> iSurf = 3
    }else{
        _error("computePressureForceAt: incompatible load surface must be 1, 2 or 3");
    }   

	    m.setValues(3,  genEps.at(13),  genEps.at(14),  genEps.at(15) );  gam =  genEps.at(18);

        if( surfLoad->giveClassID() == ConstantPressureLoadClass ){
            this->evalCovarBaseVectorsAt(gp, g1, g2, g3, tStep, genEps); // m=g3
            surfLoad->computeValueAt(load, tStep, * ( gp->giveCoordinates() ), VM_Total); // pressure components
            traction.beVectorProductOf(g1,g2); //traction.normalize(); // normal vector
            traction.times( -load.at(1) );
        }else if( surfLoad->giveClassID() == ConstantSurfaceLoadClass ){
            surfLoad->computeValueAt(traction, tStep, * ( gp->giveCoordinates() ), VM_Total); // traction vector

        }else{
            _error("computePressureForceAt: incompatible load type");
        }
        
        
        double fac1 = zeta + 0.5*zeta*zeta*gam;
        double fac2 = 0.5*zeta*zeta;
        /*  Force
        fp = [      n
                   f1*n
              f2*dotprod(n,m)]
        */
        answer.at(1) = traction.at(1);
        answer.at(2) = traction.at(2);
        answer.at(3) = traction.at(3);
        answer.at(4) = fac1 * traction.at(1);
        answer.at(5) = fac1 * traction.at(2);
        answer.at(6) = fac1 * traction.at(3);
        answer.at(7) = fac2 * m.dotProduct(traction);
    

}

void 
TrDirShell :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(42);
    for( int i = 1; i <= 42; i++){    
        answer.at(i) = i;
    }
    //answer.setValues(42,  1, 2, 3, 4, 5, 6, 7,        22, 23, 24, 25, 26, 27, 28,   8, 9, 10, 11, 12, 13, 14,
    //                     29, 30, 31, 32, 33, 34, 35,  15, 16, 17, 18, 19, 20, 21,  36, 37, 38, 39, 40, 41, 42);
}



// Mass matrix
void
TrDirShell :: computeThicknessMappingCoeff(GaussPoint *gp, FloatArray &answer){
	//thickness jacobian = ratio between volume and area: j0 = a3 + a2*zeta^2 + a1 * zeta
	// Returns array with a1-a3, used in expression for analytical integration of mass matrix.
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FloatArray lcoords = *gp->giveCoordinates();
	
    FloatMatrix dNdxi;
    this->interpolation.evaldNdxi(dNdxi, lcoords);

	FloatArray M, dM1(3), dM2(3), dX1(3), dX2(3);
    this->evalInitialDirectorAt(gp, M); 

	for (int i = 1; i <= 6; i++ ) {
		double x, y, z, Mx, My, Mz;
        FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
		x = nodeI->at(1); y = nodeI->at(2); z = nodeI->at(3);
        
		M=this->giveInitialNodeDirector(i);
		Mx = M.at(1); My = M.at(2); Mz = M.at(3);

        dX1.at(1) += dNdxi.at(i,1)* x; 
		dX1.at(2) += dNdxi.at(i,1)* y; 
		dX1.at(3) += dNdxi.at(i,1)* z;
        dM1.at(1) += dNdxi.at(i,1)* Mx; 
		dM1.at(2) += dNdxi.at(i,1)* My; 
		dM1.at(3) += dNdxi.at(i,1)* Mz; 

        dX2.at(1) += dNdxi.at(i,2)* x; 
		dX2.at(2) += dNdxi.at(i,2)* y; 
		dX2.at(3) += dNdxi.at(i,2)* z;
        dM2.at(1) += dNdxi.at(i,2)* Mx; 
		dM2.at(2) += dNdxi.at(i,2)* My; 
		dM2.at(3) += dNdxi.at(i,2)* Mz; 

    }

    double sc;
    FloatArray temp, temp2;
    temp.beVectorProductOf(dX1,dX2);
    sc = temp.computeNorm();
    answer.resize(3);
    answer.at(3) = M.dotProduct(temp)/sc;
    
    temp.beVectorProductOf(dX1,dM2); 
    temp2.beVectorProductOf(dX2,dM1); 
    temp.add(temp2);
    answer.at(2) = M.dotProduct(temp)/sc;

    temp.beVectorProductOf(dM1,dM2); 
    answer.at(1) = M.dotProduct(temp)/sc;

}




// Integration

double 
TrDirShell :: computeVolumeAround(GaussPoint *gp){
	FloatArray G1, G2, G3, temp;
	double detJ;
	this->evalInitialCovarBaseVectorsAt(gp,G1, G2, G3);
	temp.beVectorProductOf(G1, G2);
    detJ = temp.dotProduct(G3)*0.5*this->giveCrossSection()->give(CS_Thickness); 
    return detJ * gp->giveWeight();
}

double 
TrDirShell :: computeAreaAround(GaussPoint *gp){
	FloatArray G1, G2, G3, temp;
	double detJ;
	this->evalInitialCovarBaseVectorsAt(gp, G1, G2, G3);
	temp.beVectorProductOf(G1, G2);
    detJ = temp.computeNorm();
    return detJ * gp->giveWeight()*0.5 ;
}

double 
TrDirShell :: edgeComputeAreaAround(GaussPoint *gp, const int iedge){
	FloatArray G1, G3, temp;
	double detJ;
	this->edgeEvalInitialCovarBaseVectorsAt(gp, iedge,G1, G3);
    detJ = G1.computeNorm();
    return detJ * gp->giveWeight() ;
}

double 
TrDirShell :: computeThicknessAroundLayer(GaussPoint *gp, int layer){
	FloatArray G1, G2, G3, temp;
	double detJ;
	this->evalInitialCovarBaseVectorsAt(gp,G1, G2, G3);
	temp.beVectorProductOf(G1, G2);
    //detJ = temp.dotProduct(G3)*0.5*this->giveCrossSection()->give(CS_Thickness); 
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >(this->giveCrossSection());
    detJ = temp.dotProduct(G3)*0.5*layeredCS->giveLayerThickness(layer);
    return detJ * gp->giveWeight();

}

double 
TrDirShell :: computeVolumeAroundLayer(GaussPoint *gp, int layer){
	FloatArray G1, G2, G3, temp;
	double detJ;
	this->evalInitialCovarBaseVectorsAt(gp,G1, G2, G3);
	temp.beVectorProductOf(G1, G2);
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >(this->giveCrossSection());
    detJ = temp.dotProduct(G3)*0.5*layeredCS->giveLayerThickness(layer);
    return detJ * gp->giveWeight();

}

FEInterpolation *TrDirShell :: giveInterpolation()
{
    return &interpolation;
}


void 
TrDirShell :: computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep){
	// Analytically integrated over the thickness. Constant density assumed.
    // => integration rule #2
	IntegrationRule *iRule = integrationRulesArray [ 0 ];
	GaussPoint *gp;

	//------------------------------
    FloatMatrix N, Nt, Ntm, NtmN, mass;
    FloatArray a, unknowns, m(3);
    double gam, zeta;
    this->giveUpdatedSolutionVector(a, tStep);  

	answer.resize(42,42); answer.zero();
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {

        gp = iRule->getIntegrationPoint(i);
        FloatArray lcoords = *gp->giveCoordinates(); 
        zeta = giveLocalZetaCoord(gp);

		this->computeNmatrixAt(gp, N);
		unknowns.beProductOf(N,a); // [x, m, gam]^T
		m.setValues(3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
		gam = unknowns.at(7);
		
		// Consistent mass matrix
		/*     3    3    1
  		   3 [a*I  b*I   c*m      [A  B  C
		   3       d*I   e*m    =     D  E
		   1  sym       f*m.m]     sym   F]
		*/

		double h, h2, h3, h5, fac1, fac2, fac3, fac4, fac5, fac6, gam2, rho, dV; 
		rho = this->giveMaterial()->give('d', gp);
		h   = this->giveCrossSection()->give(CS_Thickness);
		h2  = h*h; h3 = h2*h; h5 = h2*h3;
		gam2 = gam*gam;

		mass.resize(7,7);
		fac1 = 4;
		fac2 = 2.0 * zeta*(2.0 + gam*zeta);
		fac3 = 2.0 * zeta*zeta;
		fac4 = zeta * zeta*(2.0 + gam*zeta) * (2.0 + gam*zeta);
		fac5 = zeta*zeta*zeta*(2.0 + gam*zeta);
		fac6 = zeta*zeta*zeta*zeta;
		mass.at(1,1) = mass.at(2,2) = mass.at(3,3) = fac1; // A
		mass.at(1,4) = mass.at(2,5) = mass.at(3,6) = fac2; // B
		mass.at(1,7) = fac3 * m.at(1);  
		mass.at(2,7) = fac3 * m.at(2); 
		mass.at(3,7) = fac3 * m.at(3); // C
		mass.at(4,4) = mass.at(5,5) = mass.at(6,6) = fac4; // D
		mass.at(4,7) = fac5 * m.at(1);  
		mass.at(5,7) = fac5 * m.at(2); 
		mass.at(6,7) = fac5 * m.at(3); // E
		mass.at(7,7) = fac6 * m.dotProduct(m); // F
        mass.symmetrized();

		dV = this->computeVolumeAround(gp) ; 
        Ntm.beTProductOf(N, mass);
		NtmN.beProductOf(Ntm, N);
        NtmN.times(dV*rho/4.0 );
		answer.add(NtmN);
        
	}


}




// More general purpose methods

void
TrDirShell :: giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3,
                                                 FloatArray &G1, FloatArray &G2, FloatArray &G3){
    answer.resize(3,3);
    // Q_ij = dotproduct(g_i, G_j) 
    answer.at(1,1) = g1.dotProduct(G1); answer.at(1,2) = g1.dotProduct(G2); answer.at(1,3) = g1.dotProduct(G3);
    answer.at(2,1) = g2.dotProduct(G1); answer.at(2,2) = g2.dotProduct(G2); answer.at(2,3) = g2.dotProduct(G3);
    answer.at(3,1) = g3.dotProduct(G1); answer.at(3,2) = g3.dotProduct(G2); answer.at(3,3) = g3.dotProduct(G3);
    
}

void
TrDirShell :: giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3){
    // Assumes a transformation from a cartesian system
    answer.resize(3,3);
    // Q_ij = dotproduct(g_i, E_j) with E1 = [1 0 0]^T etc
    answer.at(1,1) = g1.at(1); answer.at(1,2) = g1.at(2); answer.at(1,3) = g1.at(3);
    answer.at(2,1) = g2.at(1); answer.at(2,2) = g2.at(2); answer.at(2,3) = g2.at(3);
    answer.at(3,1) = g3.at(1); answer.at(3,2) = g3.at(2); answer.at(3,3) = g3.at(3);
    
}

void
TrDirShell :: giveBondTransMatrix(FloatMatrix &answer, FloatMatrix &Q){
/* Returns the Bond transformation matrix (M) between coordinate system 1 and 2, defined by the three base 
   vectors in each coordinate system g1i and g2i respectively. 
   Q is a transformation matrix defined as the direction cosines between the transformed system and the initial
   system.*/


    answer.resize(6,6);
    // block matrix 11 = Q^2
    answer.at(1,1) = Q.at(1,1)*Q.at(1,1); answer.at(1,2) = Q.at(1,2)*Q.at(1,2); answer.at(1,3) = Q.at(1,3)*Q.at(1,3);
    answer.at(2,1) = Q.at(2,1)*Q.at(2,1); answer.at(2,2) = Q.at(2,2)*Q.at(2,2); answer.at(2,3) = Q.at(2,3)*Q.at(2,3);
    answer.at(3,1) = Q.at(3,1)*Q.at(3,1); answer.at(3,2) = Q.at(3,2)*Q.at(3,2); answer.at(3,3) = Q.at(3,3)*Q.at(3,3);

    // block matrix 12 
    answer.at(1,4) = 2.*Q.at(1,2)*Q.at(1,3); answer.at(1,5) = 2.*Q.at(1,3)*Q.at(1,1); answer.at(1,6) = 2.*Q.at(1,1)*Q.at(1,2);
    answer.at(2,4) = 2.*Q.at(2,2)*Q.at(2,3); answer.at(2,5) = 2.*Q.at(2,3)*Q.at(2,1); answer.at(2,6) = 2.*Q.at(2,1)*Q.at(2,2);
    answer.at(3,4) = 2.*Q.at(3,2)*Q.at(3,3); answer.at(3,5) = 2.*Q.at(3,3)*Q.at(3,1); answer.at(3,6) = 2.*Q.at(3,1)*Q.at(3,2);

    // block matrix 21 
    answer.at(4,1) = Q.at(2,1)*Q.at(3,1); answer.at(4,2) = Q.at(2,2)*Q.at(3,2); answer.at(4,3) = Q.at(2,3)*Q.at(3,3);
    answer.at(5,1) = Q.at(3,1)*Q.at(1,1); answer.at(5,2) = Q.at(3,2)*Q.at(1,2); answer.at(5,3) = Q.at(3,3)*Q.at(1,3);
    answer.at(6,1) = Q.at(1,1)*Q.at(2,1); answer.at(6,2) = Q.at(1,2)*Q.at(2,2); answer.at(6,3) = Q.at(1,3)*Q.at(2,3);

    // block matrix 22 
    answer.at(4,4) = Q.at(2,2)*Q.at(3,3) + Q.at(2,3)*Q.at(3,2); 
    answer.at(4,5) = Q.at(2,1)*Q.at(3,3) + Q.at(2,3)*Q.at(3,1); 
    answer.at(4,6) = Q.at(2,2)*Q.at(3,1) + Q.at(2,1)*Q.at(3,2);
    answer.at(5,4) = Q.at(1,2)*Q.at(3,3) + Q.at(1,3)*Q.at(3,2);
    answer.at(5,5) = Q.at(1,3)*Q.at(3,1) + Q.at(1,1)*Q.at(3,3);
    answer.at(5,6) = Q.at(1,1)*Q.at(3,2) + Q.at(1,2)*Q.at(3,1);
    answer.at(6,4) = Q.at(1,2)*Q.at(2,3) + Q.at(1,3)*Q.at(2,2);
    answer.at(6,5) = Q.at(1,3)*Q.at(2,1) + Q.at(1,1)*Q.at(2,3);
    answer.at(6,6) = Q.at(1,1)*Q.at(2,2) + Q.at(1,2)*Q.at(2,1);

}



void 
TrDirShell :: computeLinearizedStiffness(GaussPoint *gp, Material *mat, TimeStep *tStep, 
                FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, FloatMatrix A[3][3], FloatArray &genEps){


    FloatArray lcoords, cartStressVector, contravarStressVector, sectionalForces, BF, a;
    FloatArray g1, g2, g3, m(3), dm1(3), dm2(3), temp1, temp2;

    //Material *mat = this->giveMaterial();
    FloatMatrix D, Dcart, Mmat, S;

    //A = L^iklj * (g_k x g_l) + S^ij*I
    mat->giveCharacteristicMatrix(Dcart, FullForm, TangentStiffness, gp, tStep); // L_ijkl - cartesian system (Voigt)
    this->transInitialCartesianToInitialContravar(gp, Dcart, D);      // L^ijkl - curvilinear system (Voigt)

         

    this->computeStressVector(cartStressVector, genEps, gp, mat, tStep);
    this->transInitialCartesianToInitialContravar(gp, cartStressVector, contravarStressVector);
    S.beMatrixForm(contravarStressVector);

    this->computeStressResultantsAt(gp, contravarStressVector, S1g, S2g, S3g, tStep, genEps);  
    this->evalCovarBaseVectorsAt(gp, g1, g2, g3, tStep, genEps);

    FloatMatrix gg11, gg12, gg13, gg21, gg22, gg23, gg31, gg32, gg33;


    gg11.beDyadicProductOf(g1,g1); gg12.beDyadicProductOf(g1,g2); gg13.beDyadicProductOf(g1,g3); 
    gg21.beDyadicProductOf(g2,g1); gg22.beDyadicProductOf(g2,g2); gg23.beDyadicProductOf(g2,g3); 
    gg31.beDyadicProductOf(g3,g1); gg32.beDyadicProductOf(g3,g2); gg33.beDyadicProductOf(g3,g3); 



    // position 11

     A[0][0].resize(3,3); A[0][0].beUnitMatrix();  A[0][0].times( S.at(1,1) );
     A[0][0].add(D.at(1,1),gg11);  A[0][0].add(D.at(1,6),gg12);  A[0][0].add(D.at(1,5),gg13); 
     A[0][0].add(D.at(6,1),gg21);  A[0][0].add(D.at(6,6),gg22);  A[0][0].add(D.at(6,5),gg23); 
     A[0][0].add(D.at(5,1),gg31);  A[0][0].add(D.at(5,6),gg32);  A[0][0].add(D.at(5,5),gg33); 
     
     
    // position 21
    A[1][0].resize(3,3); A[1][0].beUnitMatrix(); A[1][0].times( S.at(2,1) );
    A[1][0].add(D.at(6,1),gg11); A[1][0].add(D.at(6,6),gg12); A[1][0].add(D.at(6,5),gg13); 
    A[1][0].add(D.at(2,1),gg21); A[1][0].add(D.at(2,6),gg22); A[1][0].add(D.at(2,5),gg23); 
    A[1][0].add(D.at(4,1),gg31); A[1][0].add(D.at(4,6),gg32); A[1][0].add(D.at(4,5),gg33); 
    
    // position 31
    A[2][0].resize(3,3); A[2][0].beUnitMatrix(); A[2][0].times( S.at(3,1) );
    A[2][0].add(D.at(5,1),gg11); A[2][0].add(D.at(5,6),gg12); A[2][0].add(D.at(5,5),gg13); 
    A[2][0].add(D.at(4,1),gg21); A[2][0].add(D.at(4,6),gg22); A[2][0].add(D.at(4,5),gg23); 
    A[2][0].add(D.at(3,1),gg31); A[2][0].add(D.at(3,6),gg32); A[2][0].add(D.at(3,5),gg33); 
    
    // position 12 
    A[0][1].resize(3,3); A[0][1].beUnitMatrix(); A[0][1].times( S.at(1,2) );
    A[0][1].add(D.at(1,6),gg11); A[0][1].add(D.at(1,2),gg12); A[0][1].add(D.at(1,4),gg13); 
    A[0][1].add(D.at(6,6),gg21); A[0][1].add(D.at(6,2),gg22); A[0][1].add(D.at(6,4),gg23); 
    A[0][1].add(D.at(5,6),gg31); A[0][1].add(D.at(5,2),gg32); A[0][1].add(D.at(5,4),gg33); 
    
    // position 22
    A[1][1].resize(3,3); A[1][1].beUnitMatrix(); A[1][1].times( S.at(2,2) );
    A[1][1].add(D.at(6,6),gg11); A[1][1].add(D.at(6,2),gg12); A[1][1].add(D.at(6,4),gg13); 
    A[1][1].add(D.at(2,6),gg21); A[1][1].add(D.at(2,2),gg22); A[1][1].add(D.at(2,4),gg23); 
    A[1][1].add(D.at(4,6),gg31); A[1][1].add(D.at(4,2),gg32); A[1][1].add(D.at(4,4),gg33); 
    
    // position 32 
    A[2][1].resize(3,3); A[2][1].beUnitMatrix(); A[2][1].times( S.at(3,2) );
    A[2][1].add(D.at(5,6),gg11); A[2][1].add(D.at(5,2),gg12); A[2][1].add(D.at(5,4),gg13); 
    A[2][1].add(D.at(4,6),gg21); A[2][1].add(D.at(4,2),gg22); A[2][1].add(D.at(4,4),gg23); 
    A[2][1].add(D.at(3,6),gg31); A[2][1].add(D.at(3,2),gg32); A[2][1].add(D.at(3,4),gg33); 
    
    // position 13 
    A[0][2].resize(3,3); A[0][2].beUnitMatrix(); A[0][2].times( S.at(1,3) );
    A[0][2].add(D.at(1,5),gg11); A[0][2].add(D.at(1,4),gg12); A[0][2].add(D.at(1,3),gg13); 
    A[0][2].add(D.at(6,5),gg21); A[0][2].add(D.at(6,4),gg22); A[0][2].add(D.at(6,3),gg23); 
    A[0][2].add(D.at(5,5),gg31); A[0][2].add(D.at(5,4),gg32); A[0][2].add(D.at(5,3),gg33); 
    
    // position 23 
    A[1][2].resize(3,3); A[1][2].beUnitMatrix(); A[1][2].times( S.at(2,3) );
    A[1][2].add(D.at(6,5),gg11); A[1][2].add(D.at(6,4),gg12); A[1][2].add(D.at(6,3),gg13); 
    A[1][2].add(D.at(2,5),gg21); A[1][2].add(D.at(2,4),gg22); A[1][2].add(D.at(2,3),gg23); 
    A[1][2].add(D.at(4,5),gg31); A[1][2].add(D.at(4,4),gg32); A[1][2].add(D.at(4,3),gg33); 
    
    // position 33 
    A[2][2].resize(3,3); A[2][2].beUnitMatrix(); A[2][2].times( S.at(3,3) );
    A[2][2].add(D.at(5,5),gg11); A[2][2].add(D.at(5,4),gg12); A[2][2].add(D.at(5,3),gg13); 
    A[2][2].add(D.at(4,5),gg21); A[2][2].add(D.at(4,4),gg22); A[2][2].add(D.at(4,3),gg23); 
    A[2][2].add(D.at(3,5),gg31); A[2][2].add(D.at(3,4),gg32); A[2][2].add(D.at(3,3),gg33); 


   // printf("gg11 \n"); gg11.printYourself(); printf("gg12 \n"); gg12.printYourself(); printf("gg13 \n"); gg13.printYourself();
   // printf("gg21 \n"); gg21.printYourself(); printf("gg22 \n"); gg22.printYourself(); printf("gg23 \n"); gg23.printYourself();
   // printf("gg31 \n"); gg31.printYourself(); printf("gg32 \n"); gg32.printYourself(); printf("gg33 \n"); gg32.printYourself();


    //printf("A11 \n"); A[0][0].printYourself(); //printf("A12 \n"); A[0][1].printYourself(); printf("A13 \n"); A[0][2].printYourself();
    //printf("A21 \n"); A[1][0].printYourself(); printf("A22 \n"); A[1][1].printYourself(); printf("A23 \n"); A[1][2].printYourself();
    //printf("A31 \n"); A[2][0].printYourself(); printf("A32 \n"); A[2][1].printYourself(); printf("A33 \n"); A[2][2].printYourself();


    }




int 
TrDirShell :: giveVoigtIndex(const int ind1, const int ind2){
     
    if( ind1==1 && ind2==1 ){
        return 1;
    }else if( ind1==2 && ind2==2){
        return 2;
    }else if( ind1==3 && ind2==3){
        return 3;   
    }else if( (ind1==2 && ind2==3) || (ind1==3 && ind2==2) ){
        return 4;
    }else if( (ind1==1 && ind2==3) || (ind1==3 && ind2==1) ){
        return 5;
    }else if( (ind1==1 && ind2==2) || (ind1==2 && ind2==1) ){
        return 6;
    }else{
        OOFEM_ERROR("Error in giveVoigtIndex - bad indices");
        return -1;
    }




    };


void
TrDirShell :: giveTensorForm(const FloatMatrix &matrix, FloatArray &tensor){

        
    // convert 4th order tensor in Voigt to tensor form

    FloatArray L [1][1][1][1];
    int ind1, ind2;
    for( int i = 0; i<=2; i++ ){
        for( int j = 0; j<=2; j++ ){
            for( int k = 0; k<=2; k++ ){
                for( int l = 0; l<=2; l++ ){
                    
                    ind1 = giveVoigtIndex(i+1,j+1);
                    ind2 = giveVoigtIndex(k+1,l+1);
                    L[i][j][k][l].resize(1);
                    L[i][j][k][l] = matrix.at(ind1,ind2);

                }
            }
        }
    }

}



void 
TrDirShell :: compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer){

    
    
    answer.resize(18,18);
   for( int i = 1; i<=18; i++ ){
       for( int j = 1; j<=18; j++ ){

           if( abs(matrix1.at(i,j)) > 1.0e-12 ){
               double diff = ( matrix1.at(i,j)-matrix2.at(i,j) );
               double relDiff =  diff / matrix1.at(i,j);
               if( abs(relDiff)<1.0e-4){
                   answer.at(i,j) = 0.0;
               }else if( abs(diff)<1.0e3 ){
                   answer.at(i,j) = 0.0;
               }else
                   answer.at(i,j) = relDiff;
           }else{
               answer.at(i,j) = -1.0;
           }

       }
   }

}




Interface *TrDirShell :: giveInterface(InterfaceType it)
{
    switch (it) {
        case NodalAveragingRecoveryModelInterfaceType:
            return static_cast< NodalAveragingRecoveryModelInterface * >(this);//c++
        case LayeredCrossSectionInterfaceType:
            return ( LayeredCrossSectionInterface * ) this;
        default:
            return StructuralElement :: giveInterface(it);
    }
}


void TrDirShell :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

int TrDirShell :: NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( type == IST_DirectorField ) {
        return 3;
    }

    return 0;
}


void TrDirShell :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DirectorField ) {
        answer.resize(3);
        answer = this->giveInitialNodeDirector( node );
        answer.at(1) += this->giveNode(node)->giveDofWithID(w_u)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        answer.at(2) += this->giveNode(node)->giveDofWithID(w_v)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        answer.at(3) += this->giveNode(node)->giveDofWithID(w_w)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        answer.times(this->giveCrossSection()->give(CS_Thickness));
    } else   {
        answer.resize(0);
    }
}




} // end namespace oofem

