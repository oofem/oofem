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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "Shell7BasePhFi.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "feinterpol3d.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "constantpressureload.h"
#include "constantsurfaceload.h"
#include "vtkxmlexportmodule.h"
#include "fracturemanager.h"

#include "masterdof.h"

#include <fstream>

namespace oofem {

const int nLayers = 5; //@todo: Generalize!
const double disturB = 1e-8; //@todo: Generalize!


Shell7BasePhFi :: Shell7BasePhFi(int n, Domain *aDomain) : Shell7Base(n, aDomain), PhaseFieldElement(n, aDomain){
	this->numberOfLayers = nLayers;
}

IRResultType Shell7BasePhFi :: initializeFrom(InputRecord *ir)
{
    Shell7Base :: initializeFrom(ir);
    return IRRT_OK;
}



void
Shell7BasePhFi :: postInitialize()
{
    
    Shell7Base :: postInitialize();

}



void
Shell7BasePhFi :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // returns the total DofIDMask: [shell7Base IDMask, damage var. IdMask]
	IntArray answer_d;
    this->giveDofManDofIDMask_u(answer);
	this->giveDofManDofIDMask_d(answer_d);
	answer.followedBy(answer_d);

}

void
Shell7BasePhFi :: giveDofManDofIDMask_u(IntArray &answer) const
{
	answer.setValues(7, D_u, D_v, D_w, W_u, W_v, W_w, Gamma);
}

void
Shell7BasePhFi :: giveDofManDofIDMask_d(IntArray &answer) const
{
    // Returns [nextFreeDofID, nextFreeDofID+1, ..., nextFreeDofID + numberOfLayers]
    answer.resize(0);
    int sID = this->domain->giveNextFreeDofID(0);
	for (int i = 0; i < numberOfLayers; i++) {
		answer.followedBy(sID + i);
	}
}


double 
Shell7BasePhFi :: computeDamageInLayerAt(int layer, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN)
{
    // d = N_d * a_d
    FloatArray dVec, dVecRed;
	
    this->computeDamageUnknowns(dVec, valueMode, stepN);		// should be a vector with all damage nodal values for this element
																// ordered such that all damage dofs associated with node 1 comes first
    FloatArray Nvec, lcoords;					
	IntArray indx = computeDamageIndexArray(layer);         // Since dVec only contains damage vars this index vector should be 1 3 5 7 9 11 for the first layer (with 2 layers in the cs)

	dVecRed.beSubArrayOf(dVec, indx);
    this->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    
    // Make sure returned damage is always between [0,1]
    double d_temp = Nvec.dotProduct(dVecRed);
    //if ( d_temp < 0.0 ) {
    //    return 0.0;
    //} else if ( d_temp > 1.0 ) {
    //    return 1.0;
    //} else {
        return d_temp;
    //}
}

double 
Shell7BasePhFi :: computeOldDamageInLayerAt(int layer, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN)
{
    // d = N_d * a_d
    FloatArray dVec, dVecRed, ddVec;
	
    this->computeDamageUnknowns(dVec, valueMode, stepN);		// should be a vector with all damage nodal values for this element
																// ordered such that all damage dofs associated with node 1 comes first
	this->computeDamageUnknowns(ddVec, VM_Incremental, stepN);

	dVec.subtract(ddVec);

    FloatArray Nvec, lcoords;					
	IntArray indx = computeDamageIndexArray(layer);         // Since dVec only contains damage vars this index vector should be 1 3 5 7 9 11 for the first layer (with 2 layers in the cs)

	dVecRed.beSubArrayOf(dVec, indx);
    this->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    
    // Make sure returned damage is always between [0,1]
    double d_temp = Nvec.dotProduct(dVecRed);
    //if ( d_temp < 0.0 ) {
    //    return 0.0;
    //} else if ( d_temp > 1.0 ) {
    //    return 1.0;
    //} else {
        return d_temp;
    //}
}
double 
Shell7BasePhFi :: computeDamageInLayerAt_dist(int layer, int index, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN)
{
    // d = N_d * a_d
    FloatArray dVec, dVecRed;
	
    this->computeDamageUnknowns(dVec, valueMode, stepN);		// should be a vector with all damage nodal values for this element
																// ordered such that all damage dofs associated with node 1 comes first
	if (valueMode == VM_Total)
	{
		dVec.at(index) = dVec.at(index) + disturB;
	} else if (valueMode == VM_Incremental)
	{
		dVec.at(index) = dVec.at(index) + disturB/(stepN->giveTimeIncrement());
	}
	

    FloatArray Nvec, lcoords;					
	IntArray indx = computeDamageIndexArray(layer);         // Since dVec only contains damage vars this index vector should be 1 3 5 7 9 11 for the first layer (with 2 layers in the cs)

	dVecRed.beSubArrayOf(dVec, indx);
    this->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    
    // Make sure returned damage is always between [0,1]
    double d_temp = Nvec.dotProduct(dVecRed);
   // if ( d_temp < 0.0 ) {
   //     return 0.0;
   // } else if ( d_temp > 1.0 ) {
   //     return 1.0;
   // } else {
        return d_temp;
   // }
}

#if 0
double 
Shell7BasePhFi :: computeDamageInLayerAt_dist(int layer, int index, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN)
{
    // d = N_d * a_d
    FloatArray dVec, dVecRed;
	
    this->computeDamageUnknowns(dVec, valueMode, stepN);		// should be a vector with all damage nodal values for this element
																// ordered such that all damage dofs associated with node 1 comes first
	dVec.at(index) = dVec.at(index) + disturB;
    FloatArray Nvec, lcoords;					
	IntArray indx = computeDamageIndexArray(layer);         // Since dVec only contains damage vars this index vector should be 1 3 5 7 9 11 for the first layer (with 2 layers in the cs)

	dVecRed.beSubArrayOf(dVec, indx);
    this->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    
    // Make sure returned damage is always between [0,1]
    double d_temp = Nvec.dotProduct(dVecRed);
    if ( d_temp < 0.0 ) {
        return 0.0;
    } else if ( d_temp > 1.0 ) {
        return 1.0;
    } else {
        return d_temp;
    }
}
#endif

IntArray
Shell7BasePhFi :: computeDamageIndexArray(int layer)
{
	int numberOfNodes = this->giveNumberOfDofManagers();		
	
	IntArray indx(numberOfNodes);

	for (int i = 1; i <= numberOfNodes; i++)
	{
        indx.at(i) = layer + (i-1)*this->layeredCS->giveNumberOfLayers(); // NEW JB
	}
    return indx;
}

double
Shell7BasePhFi  :: computeGInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // computes g = (1-d)^2 + r0
    //double d = this->computeDamageInLayerAt(layer, gp, valueMode, stepN);
	double d = this->computeOldDamageInLayerAt(layer, gp, valueMode, stepN);
    //double d = this->computeDamageInLayerAt(layer, gp, valueMode, stepN);
	double r0 = 1.0e-10;

    return (1.0 - d) * (1.0 - d) + r0;  
    
}

double
Shell7BasePhFi  :: computeGInLayer_dist(int layer, int indx, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // computes g = (1-d)^2 + r0
    double d = this->computeDamageInLayerAt_dist(layer, indx, gp, valueMode, stepN);
    double r0 = 1.0e-10;

    return (1.0 - d) * (1.0 - d) + r0;  
    
}

double 
Shell7BasePhFi  :: computeGprimInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // compute g' =-2*(1-d)
    double d = this->computeDamageInLayerAt(layer, gp, valueMode, stepN);
    return -2.0 * (1.0 - d);
    
}

double 
Shell7BasePhFi  :: computeGprimInLayer_dist(int layer, int indx, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // compute g' =-2*(1-d)
    double d = this->computeDamageInLayerAt_dist(layer, indx, gp, valueMode, stepN);
    return -2.0 * (1.0 - d);
    
}

double 
Shell7BasePhFi  :: computeGbisInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // compute g'' = D(-2*(1-d))/Dd = 2 
    //@todo this method is trivial but here if one wants to change the expression for g
    //return 2.0;
    double d = this->computeDamageInLayerAt(layer, gp, valueMode, stepN);
    //if ( d < 0.0 ) {
    //    return 0.0;
    //} else if ( d > 1.0 ) {
    //    return 0.0;
    //} else {
        return 2.0;
    //}
}

int
Shell7BasePhFi :: giveNumberOfDofs() 
{
    return this->giveNumberOfuDofs() + this->giveNumberOfdDofs();
}

int
Shell7BasePhFi :: giveNumberOfuDofs() 
{
    return 7 * this->giveNumberOfDofManagers();
}

int
Shell7BasePhFi :: giveNumberOfdDofs() 
{
    return this->giveNumberOfDofManagers() * this->numberOfLayers;
}


// Tangent matrices

void
Shell7BasePhFi :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
	Shell7Base :: computeStiffnessMatrix(answer, rMode, tStep);
}


void
Shell7BasePhFi :: computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) {
         
    int ndofs   = Shell7BasePhFi :: giveNumberOfDofs();
	int ndofs_u = Shell7BasePhFi :: giveNumberOfuDofs(); 
	int ndofs_d = ndofs - ndofs_u;

    answer.resize( ndofs_u, ndofs_d);
    answer.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
  
	FloatArray fu(ndofs_u), fd(ndofs_d), ftemp_u, sectionalForces, sectionalForces_dv;
    FloatArray genEps, genEpsD, solVec, lCoords, Nd;
    FloatMatrix B, Bd, K_ud(ndofs_u, this->numberOfDofMans), K_temp;
    this->giveUpdatedSolutionVector(solVec, tStep);    // => x, m, gam

    IntArray ordering_temp(ndofs_u);
    for ( int i = 1; i <= ndofs_u; i++ ) {
        ordering_temp.at(i) = i;
    }

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        fu.zero();
        K_ud.zero();
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );
		IntArray indx_d = computeDamageIndexArray(layer);
        FloatArray dVecTotal, dVecLayer, Ddam_Dxi;
	
        this->computeDamageUnknowns(dVecTotal, VM_Total, tStep);		// vector with all damage nodal values for this element
	    dVecLayer.beSubArrayOf(dVecTotal, indx_d);                      // vector of damage variables for a given layer        
       
        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();
            this->computeBmatrixAt(lCoords, B);
			this->computeNdvectorAt(lCoords, Nd); // (1x6)

            this->computeGeneralizedStrainVectorNew(genEpsD, solVec, B);
            this->computeGeneralizedStrainVectorNew(genEps, solVec, B); // used for computing the stress

            double zeta = giveGlobalZcoord( *gp->giveCoordinates() ); 
            this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEpsD, zeta); // these are per unit volume
				
            // Computation of tangent: K_ud = \int Bu^t* (sec. forces) * g' * Nd
            ftemp_u.beTProductOf(B,sectionalForces);
            double dV = this->computeVolumeAroundLayer(gp, layer);
            double Gprim   = computeGprimInLayer(layer, gp, VM_Total, tStep);
			
			fu.add(dV*Gprim, ftemp_u);
            K_temp.beDyadicProductOf(fu, Nd);
            K_ud.add(K_temp);
        }
        
        // Just place the columns in the right order
        answer.assemble(K_ud, ordering_temp, indx_d);
		//ordering_temp.printYourself();
		//indx_d.printYourself();
    }

}

void
Shell7BasePhFi :: computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) {

    answer.resize( this->numberOfDofMans*this->layeredCS->giveNumberOfLayers(), 42 );
    answer.zero();

}

void
Shell7BasePhFi :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
{    
    // Computation of tangent: K_dd = \int Nd^t * ( -kp*neg_Mac'(alpha_dot)/delta_t + g_c/l + G''*Psibar) * Nd + 
    //                                \int Bd^t * (  g_c * l * [G^1 G^2]^t * [G^1 G^2] ) * Bd 
    //                              = K_dd1 + K_dd2
    int ndofs   = Shell7BasePhFi :: giveNumberOfDofs();
	int ndofs_u = Shell7BasePhFi :: giveNumberOfuDofs(); 
	int ndofs_d = ndofs - ndofs_u;

    answer.resize( ndofs_d, ndofs_d );
    answer.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();    
	
    FloatArray lCoords;
    FloatMatrix Bd, Nd, temp(this->giveNumberOfDofManagers(),this->giveNumberOfDofManagers());
 
    IntArray ordering_temp(ndofs_u); // [1, 2, ..., 42]
    for ( int i = 1; i <= ndofs_u; i++ ) {
        ordering_temp.at(i) = i;
    }

    double l       = this->giveInternalLength();
    double g_c     = this->giveCriticalEnergy();
    double kp      = this->givePenaltyParameter();
    double Delta_t = tStep->giveTimeIncrement();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        temp.zero();
  
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
		IntArray indx_d = computeDamageIndexArray(layer);       

        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();
			this->computeBdmatrixAt(lCoords, Bd);        // (2x6)
			this->computeNdMatrixAt(lCoords, Nd);        // (1x6)

            double Gbis   = computeGbisInLayer(layer, gp, VM_Total, tStep);
            double Psibar  = this->computeFreeEnergy( gp, tStep );            
            double dV = this->computeVolumeAroundLayer(gp, layer);
            

            // K_dd1 = ( -kp*neg_Mac'(d_dot) / Delta_t + g_c/ l + G'' * Psibar ) * N^t*N; 
            double Delta_d = computeDamageInLayerAt(layer, gp, VM_Incremental, tStep);
            double factorN = -kp * neg_MaCauleyPrime(Delta_d/Delta_t)/Delta_t +  g_c / l + Gbis * Psibar; 
            temp.plusProductSymmUpper(Nd, Nd, factorN * dV);
            

            // K_dd2 =  g_c * l * Bd^t * Bd;
            double factorB = g_c * l;
            temp.plusProductSymmUpper(Bd, Bd, factorB * dV);
        }
        temp.symmetrized();
        answer.assemble(temp, indx_d, indx_d);
    }

}


void
Shell7BasePhFi :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
{

	IntArray loc_u, loc_d;
	FloatMatrix ansLocal;

    loc_u = this->giveOrdering_Disp();
	loc_d = this->giveOrdering_Damage();

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs, nDofs );
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4, answer5, answer6;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);     
    //this->computeStiffnessMatrix_ud(answer2, rMode, tStep);	//@todo: check why this gives worse convergence than the numerical
	this->computeStiffnessMatrix_dd(answer4, rMode, tStep);
	//this->computeStiffnessMatrixNum_ud(answer5, rMode, tStep);
	//this->computeStiffnessMatrixNum_dd(answer6, rMode, tStep); // Yields same matrix as the analytical (answer 4)
	
	//answer2.printYourself();
	//answer4.printYourself();

	
    //this->computeStiffnessMatrix_du(answer3, rMode, tStep); //symmetric
    //answer3.beTranspositionOf(answer2);
	answer3.beTranspositionOf(answer5);
 	//loc_u.printYourself();
	//loc_d.printYourself();
	
    
    answer.assemble( answer1, loc_u, loc_u );
    //answer.assemble( answer5, loc_u, loc_d );
    //answer.assemble( answer3, loc_d, loc_u );
    answer.assemble( answer4, loc_d, loc_d );
	
	//answer4.printYourself();

}

void
Shell7BasePhFi :: computeStiffnessMatrixNum_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
{    

	int ndofs   = Shell7BasePhFi :: giveNumberOfDofs();
	int ndofs_u = Shell7BasePhFi :: giveNumberOfuDofs(); 
	int ndofs_d = ndofs - ndofs_u;
	int upD = 1;

    answer.resize( ndofs_u, ndofs_d);
    answer.zero();

	FloatArray force, force_dist, forceTemp, solVec(42), force_small;
	IntArray indxDisp;



	computeSectionalForces(force, tStep, solVec, upD);
	//force.printYourself();
	
	int numberOfLayers = this->layeredCS->giveNumberOfLayers(); 
	int numdofMans = this->giveNumberOfDofManagers();

	for ( int indx = 1; indx <= numberOfLayers*numdofMans; indx++ ) {
		computeSectionalForces_dist(force_dist, indx, tStep, solVec, upD);
		//force.printYourself();
		//force_dist.printYourself();
		//solVec.printYourself();
		//force_dist.printYourself();
		force_dist.subtract(force);
		
		//force_dist.printYourself();

		const IntArray &indxDisp = this->giveOrderingPhFi(Displacement);

		force_small.beSubArrayOf(force_dist,indxDisp);
		//force_small.printYourself();
		force_small.times(1/disturB);
		answer.addSubVectorCol(force_small,1,indx);
	}
}

void
Shell7BasePhFi :: computeStiffnessMatrixNum_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
{    

	int ndofs   = Shell7BasePhFi :: giveNumberOfDofs();
	int ndofs_u = Shell7BasePhFi :: giveNumberOfuDofs(); 
	int ndofs_d = ndofs - ndofs_u;
	int upD = 1;

    answer.resize( ndofs_d, ndofs_d);
    answer.zero();

	FloatArray force, force_dist, forceTemp, solVec(42), force_small;
	IntArray indxDam;



	computeSectionalForces(force, tStep, solVec, upD);
	//force.printYourself();
	
	int numberOfLayers = this->layeredCS->giveNumberOfLayers(); 
	int numdofMans = this->giveNumberOfDofManagers();

	for ( int indx = 1; indx <= numberOfLayers*numdofMans; indx++ ) {
		computeSectionalForces_dist(force_dist, indx, tStep, solVec, upD);
		//force_dist.printYourself();
		force_dist.subtract(force);
		
		//force_dist.printYourself();

		const IntArray &indxDam = this->giveOrderingPhFi(Damage);

		force_small.beSubArrayOf(force_dist,indxDam);
		//force_small.printYourself();
		force_small.times(1/disturB);
		answer.addSubVectorCol(force_small,1,indx);
	}
}

void
Shell7BasePhFi :: new_computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep)
{
    //Shell7Base :: new_computeBulkTangentMatrix(answer, solVec, solVecI, solVecJ,rMode, tStep);
    //return;
    //@todo optimize method - remove I,J-stuff since XFEM does not need this anymore
    FloatMatrix A [ 3 ] [ 3 ], lambdaI [ 3 ], lambdaJ [ 3 ];
    FloatMatrix L(18,18), B;
    FloatMatrix K;

    int ndofs = Shell7BasePhFi :: giveNumberOfuDofs();
    answer.resize(ndofs, ndofs); 
    answer.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     
    FloatMatrix temp;
    FloatArray genEpsI, genEpsJ, genEps, lCoords;

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
        StructuralMaterial *mat = static_cast< StructuralMaterial* >( domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) ) );

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            lCoords = *gp->giveCoordinates();

            this->computeBmatrixAt(lCoords, B, 0, 0);
            this->computeGeneralizedStrainVectorNew(genEpsI, solVecI, B);
            this->computeGeneralizedStrainVectorNew(genEpsJ, solVecJ, B);
            this->computeGeneralizedStrainVectorNew(genEps , solVec , B);
            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, A, genEps);

            double zeta = giveGlobalZcoord( *gp->giveCoordinates() );
            this->computeLambdaGMatrices(lambdaI, genEpsI, zeta);
            this->computeLambdaGMatrices(lambdaJ, genEpsJ, zeta);

            // L = sum_{i,j} (lambdaI_i)^T * A^ij * lambdaJ_j
            // @todo Naive implementation - should be optimized 
            // note: L will only be symmetric if lambdaI = lambdaJ
            L.zero();
            for ( int j = 0; j < 3; j++ ) {
                for ( int k = 0; k < 3; k++ ) {
                    this->computeTripleProduct(temp, lambdaI [ j ], A [ j ][ k ], lambdaJ [ k ]);
                    L.add(temp);
                }
            }
            
            this->computeTripleProduct(K, B, L, B);
            double dV = this->computeVolumeAroundLayer(gp, layer);	
			double g = this->computeGInLayer(layer, gp, VM_Total,  tStep);				// Check so that the correct computeDamageAt method is used (in Shell7BasePhFi)
            answer.add( g*dV, K );


        }
    }

    //@todo Note! Since this block matrix is assembled outside this method with ordering_disp no assembly should be made inside here (=> double assembly) JB
    //const IntArray &ordering = this->giveOrdering_All();
    //answer.assemble(tempAnswer, ordering, ordering);
}

#if 0
int 
	Shell7BasePhFi :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
	// Compute special IST quantities this element should support
	switch (type) {
	case IST_CauchyStressTensor:
		this->computeCauchyStressVector(answer, gp, tStep);			//@todo: add damage!!!
		return 1;
	default:
		return Element :: giveIPValue(answer, gp, type, tStep);
	}
}  
#endif // 0


// Internal forces

#if 0
void
Shell7BasePhFi :: giveUpdatedSolutionVector_d(FloatArray &answer, TimeStep *tStep)
{

	IntArray dofIdArray;
    
	//Shell7Base :: giveDofManDofIDMask(dummy, EID_MomentumBalance, dofIdArray);
    this->giveDofManDofIDMask_d(dofIdArray);

	this->computeVectorOfDofIDs(dofIdArray, VM_Total, tStep, temp);
    answer.assemble( temp, this->giveOrdering(AllInv) );
}
#endif

//


void
Shell7BasePhFi :: computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord)
{
     
    int ndofs = Shell7BasePhFi :: giveNumberOfDofs();
	int ndofs_u = Shell7BasePhFi :: giveNumberOfuDofs(); 
	int ndofs_d = ndofs - ndofs_u;

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    double sectionalForces_ds = 0;
	FloatArray fu(ndofs_u), fd(ndofs_d), ftemp_u, ftemp_d(this->giveNumberOfDofManagers()), sectionalForces, sectionalForces_dv;
    FloatArray genEps, genEpsD, totalSolVec, lCoords, Nd, f_debug(this->giveNumberOfDofManagers()), F_debug(ndofs_d);
    FloatMatrix B, Bd;

    totalSolVec = solVec;
	sectionalForces.resize(2);
    FloatArray dVecTotal, dVecLayer, Ddam_Dxi, gradd;
    this->computeDamageUnknowns(dVecTotal, VM_Total, tStep);		// vector with all damage nodal values for this element 

    fu.zero();
    fd.zero();

    FloatMatrix K_dd;
    this->computeStiffnessMatrix_dd(K_dd, TangentStiffness, tStep);
    F_debug.beProductOf(K_dd, dVecTotal);

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        ftemp_d.zero();
        f_debug.zero();
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );
		IntArray indx_d = computeDamageIndexArray(layer);
        dVecLayer.beSubArrayOf(dVecTotal, indx_d);                      // vector of damage variables for a given layer        
        
        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();
            this->computeBmatrixAt(lCoords, B);
			this->computeBdmatrixAt(lCoords, Bd);      // (2x6)
			this->computeNdvectorAt(lCoords, Nd); // (1x6)
            gradd.beProductOf(Bd,dVecLayer);        // [dalpha/dX1, dalpha/dX2]

            this->computeGeneralizedStrainVectorNew(genEps, totalSolVec, B); // used for computing the stress
            double zeta = giveGlobalZcoord( *gp->giveCoordinates() ); 
            	
            // Computation of sectional forces: fu = Bu^t*[N M T Ms Ts]^t
            this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEps, zeta); // these are per unit volume
            ftemp_u.beTProductOf(B,sectionalForces);
            double dV = this->computeVolumeAroundLayer(gp, layer);
			double g = this->computeGInLayer(layer, gp, VM_Total,  tStep);
			fu.add(dV*g, ftemp_u);

            // Computation of sectional forces: fd = Nd^t * sectionalForces_ds + Bd^t * sectionalForces_dv 
            this->computeSectionalForcesAt_d(sectionalForces_ds, sectionalForces_dv, gp, mat, tStep, zeta, layer, gradd); // these are per unit volume

            // "external force" for linear problem -(-2)*psibar
            double Psibar  = this->computeFreeEnergy( gp, tStep );
            f_debug.add( -2.0*Psibar*dV,  Nd );

			ftemp_d.add( sectionalForces_ds*dV,  Nd );
            ftemp_d.plusProduct(Bd, sectionalForces_dv, dV);

        }
        // Assemble layer contribution into the correct place
        F_debug.assemble(f_debug, indx_d);
        fd.assemble(ftemp_d, indx_d);
    }
    //F_debug.printYourself();
    //fd.printYourself();
    answer.resize( ndofs );
    answer.zero();
    const IntArray &ordering_disp   = this->giveOrderingPhFi(Displacement);
	const IntArray &ordering_damage = this->giveOrderingPhFi(Damage);
    answer.assemble(fu, ordering_disp);		//ordering_disp contains only displacement related dofs, not damage
	answer.assemble(fd, ordering_damage);
    //answer.assemble(F_debug, ordering_damage);
 
}

void
Shell7BasePhFi :: computeSectionalForces_dist(FloatArray &answer, int indx, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord)
{
     
    int ndofs = Shell7BasePhFi :: giveNumberOfDofs();
	int ndofs_u = Shell7BasePhFi :: giveNumberOfuDofs(); 
	int ndofs_d = ndofs - ndofs_u;

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    double sectionalForces_ds = 0;
	FloatArray fu(ndofs_u), fd(ndofs_d), ftemp_u, ftemp_d(this->giveNumberOfDofManagers()), sectionalForces, sectionalForces_dv;
    FloatArray genEps, genEpsD, totalSolVec, lCoords, Nd;
    FloatMatrix B, Bd;
    this->giveUpdatedSolutionVector(totalSolVec, tStep);    // => x, m, gam

    fu.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        ftemp_d.zero();
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );
		IntArray indx_d = computeDamageIndexArray(layer);
        FloatArray dVecTotal, dVecLayer, Ddam_Dxi, gradd;
	
        this->computeDamageUnknowns(dVecTotal, VM_Total, tStep);		// vector with all damage nodal values for this element
		dVecTotal.at(indx) = dVecTotal.at(indx) + disturB;
	    dVecLayer.beSubArrayOf(dVecTotal, indx_d);                      // vector of damage variables for a given layer     

        //dVecLayer.printYourself();
        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();
            this->computeBmatrixAt(lCoords, B);
			this->computeBdmatrixAt(lCoords, Bd);      // (2x6)
			this->computeNdvectorAt(lCoords, Nd); // (1x6)
            //Ddam_Dxi.beProductOf(Bd,dVecLayer);        // [dalpha/dxi1, dalpha/dxi2] OLD
            gradd.beProductOf(Bd,dVecLayer);        // [dalpha/dX1, dalpha/dX2]

            this->computeGeneralizedStrainVectorNew(genEpsD, solVec, B);
            this->computeGeneralizedStrainVectorNew(genEps, totalSolVec, B); // used for computing the stress

            double zeta = giveGlobalZcoord( *gp->giveCoordinates() ); 
            	
            // Computation of sectional forces: fu = Bu^t*[N M T Ms Ts]^t
//            this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEpsD, zeta); // these are per unit volume
            this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEps, zeta); // these are per unit volume
            ftemp_u.beTProductOf(B,sectionalForces);
            double dV = this->computeVolumeAroundLayer(gp, layer);
			double g = this->computeGInLayer_dist(layer, indx, gp, VM_Total,  tStep);
			fu.add(dV*g, ftemp_u);

            // Computation of sectional forces: fd = Nd^t * sectionalForces_ds + Bd^t * sectionalForces_dv 
            //this->computeSectionalForcesAt_d(sectionalForces_ds, sectionalForces_dv, gp, mat, tStep, zeta, layer, Ddam_Dxi); // these are per unit volume // OLD
            this->computeSectionalForcesAt_d_dist(sectionalForces_ds, sectionalForces_dv, indx, gp, mat, tStep, zeta, layer, gradd); // these are per unit volume

			ftemp_d.add( sectionalForces_ds*dV,  Nd );
			//ftemp_d.plusProduct(BdX, sectionalForces_dv, dV); // OLD
            ftemp_d.plusProduct(Bd, sectionalForces_dv, dV);

        }
        // Assemble layer contribution into the correct place
        fd.assemble(ftemp_d, indx_d);

    }

    answer.resize( ndofs );
    answer.zero();
    const IntArray &ordering_disp   = this->giveOrderingPhFi(Displacement);
	const IntArray &ordering_damage = this->giveOrderingPhFi(Damage);
    answer.assemble(fu, ordering_disp);		//ordering_disp contains only displacement related dofs, not damage
	answer.assemble(fd, ordering_damage);
   
}

double 
Shell7BasePhFi :: neg_MaCauley(double par)
{
    return 0.5 * ( abs(par) - par );
}

double 
Shell7BasePhFi :: neg_MaCauleyPrime(double par)
{
    return 0.5 * ( abs(par)/(par + 1.0e-12) - 1.0 ); // 0.5*(sign - 1) taken from Ragnars code
}

void
Shell7BasePhFi :: computeSectionalForcesAt_d(double &sectionalForcesScal, FloatArray &sectionalForcesVec, IntegrationPoint *gp,
                                            Material *mat, TimeStep *tStep, double zeta, int layer, FloatArray &gradd)
{
    //PhaseFieldCrossSection *cs = static_cast...  // in the future...
    sectionalForcesVec.zero();
    //sectionalForcesScal = -kp*neg_Mac(alpha_dot) + g_c/l*d + G'*Psibar 
    double kp      = this->givePenaltyParameter();
    double Delta_t = tStep->giveTimeIncrement();
    double d       = computeDamageInLayerAt(layer, gp, VM_Total, tStep);
    double Delta_d = computeDamageInLayerAt(layer, gp, VM_Incremental, tStep);
    double l       = this->giveInternalLength();
    double g_c     = this->giveCriticalEnergy();
    double Gprim   = computeGprimInLayer(layer, gp, VM_Total, tStep);
    double Psibar  = this->computeFreeEnergy( gp, tStep );

	if (Psibar < 0) {
		OOFEM_ERROR1("Shell7BasePhFi :: computeSectionalForcesAt_d - negative strain energy predicted")
	}
    
    sectionalForcesScal = -kp * neg_MaCauley(Delta_d/Delta_t) + g_c / l * d + Gprim * Psibar;
  
    //sectionalForcesVec = grad(alpha) * g_c * l
	sectionalForcesVec.beScaled(g_c*l,gradd);
    //sectionalForcesVec = gradd * g_c * l;
    //printf("scal %e  vec %e %e\n", sectionalForcesScal, sectionalForcesVec.at(1), sectionalForcesVec.at(2));
}

void
Shell7BasePhFi :: computeSectionalForcesAt_d_dist(double &sectionalForcesScal, FloatArray &sectionalForcesVec, int indx, IntegrationPoint *gp,
                                            Material *mat, TimeStep *tStep, double zeta, int layer, FloatArray &gradd)
{
    //PhaseFieldCrossSection *cs = static_cast...  // in the future...
    
    //sectionalForcesScal = -kp*neg_Mac(alpha_dot) + g_c/l*d + G'*Psibar 
    double kp      = this->givePenaltyParameter();
    double Delta_t = tStep->giveTimeIncrement();
    double d       = computeDamageInLayerAt_dist(layer, indx, gp, VM_Total, tStep);
    double Delta_d = computeDamageInLayerAt_dist(layer, indx, gp, VM_Incremental, tStep);
    double l       = this->giveInternalLength();
    double g_c     = this->giveCriticalEnergy();
    double Gprim   = computeGprimInLayer_dist(layer, indx, gp, VM_Total, tStep);
    double Psibar  = this->computeFreeEnergy( gp, tStep );
    
    sectionalForcesScal = -kp * neg_MaCauley(Delta_d/Delta_t) + g_c / l * d + Gprim * Psibar;
    sectionalForcesVec.zero();
    //sectionalForcesVec = grad(alpha) * g_c * l
    sectionalForcesVec.beScaled(g_c*l,gradd);
    //sectionalForcesVec = gradd * g_c * l;

}

#if 0
void
Shell7BasePhFi :: computeBdmatrixAt(FloatArray &lCoords, FloatMatrix &answer)
{
    // Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
    // B*a = [Dalpha/Dxi1, Dalpha/Dxi2] 
 
    FloatMatrix dNdxi;
    this->fei->evaldNdxi( dNdxi, lCoords, FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dNdxi);


}

#else

void
Shell7BasePhFi :: computeBdmatrixAt(FloatArray &lCoords, FloatMatrix &answer)
{
    // Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
    // B*a = [Dalpha/DX1, Dalpha/DX2] 
 
    FloatMatrix dNdxi;
    this->fei->evaldNdxi( dNdxi, lCoords, FEIElementGeometryWrapper(this) );

    FloatMatrix Gcon, Gcon_red, temp;
    Shell7Base :: evalInitialContravarBaseVectorsAt(lCoords, Gcon); //[G^1 G^2 G^3]
    Gcon_red.beSubMatrixOf(Gcon,1,2,1,2);
    
    temp.beTranspositionOf(dNdxi);
    answer.beProductOf(Gcon_red, temp); // dN/dX = [G^1 G^2] * dN/dxi 

}

#endif

void 
Shell7BasePhFi :: computeNdvectorAt(const FloatArray &iLocCoords, FloatArray &answer)
{
    // Returns the matrix {N} of the receiver, evaluated at aGaussPoint. Such that
    // N*a = alpha
    this->fei->evalN( answer, iLocCoords, FEIElementGeometryWrapper(this) );
}


void 
Shell7BasePhFi :: computeNdMatrixAt(const FloatArray &iLocCoords, FloatMatrix &answer)
{
    FloatArray Nvec;
    this->computeNdvectorAt(iLocCoords, Nvec);
    answer.resize(1,Nvec.giveSize());
    for ( int i = 1; i<= Nvec.giveSize(); i++ ){
        answer.at(1,i) = Nvec.at(i);   
    }

}


// Computation of solution vectors


void
Shell7BasePhFi :: computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType u, TimeStep *stepN, FloatArray &answer)
{
    // Routine to extract the solution vector for an element given an dofid array.
    // Size will be numberOfDofs and if a certain dofId does not exist a zero is used as value. 
    // This method may e.g. be used to obtain the enriched part of the solution vector
    ///@todo: generalize so it can be used by all XFEM elements
    answer.resize( Shell7BasePhFi ::giveNumberOfDofs());
    answer.zero();
    int k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);        
        for (int j = 1; j <= dofIdArray.giveSize(); j++ ) {
            
            if ( dMan->hasDofID( (DofIDItem) dofIdArray.at(j) ) ) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at(j) );
                ///@todo: will fail if any other dof then gamma is excluded from enrichment 
                /// since I only add "j". Instead I should skip certain dof numbers when incrementing
                answer.at(k+j) = d->giveUnknown(u, stepN);	// Martin: modification to be used also for rates
            }
        }
        k += 7;
    }
}





// VTK export

void 
Shell7BasePhFi :: giveShellExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )            
{   

    int numCells = this->layeredCS->giveNumberOfLayers();
    const int numCellNodes  = 15; // quadratic wedge
    int numTotalNodes = numCellNodes*numCells;

    vtkPiece.setNumberOfCells(numCells);
    vtkPiece.setNumberOfNodes(numTotalNodes);

    std::vector <FloatArray> nodeCoords;
    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);

    // Compute fictious node coords
    int nodeNum = 1;
    for ( int layer = 1; layer <= numCells; layer++ ) {
        

        // Node coordinates
        this->giveFictiousNodeCoordsForExport(nodeCoords, layer);       
        
        for ( int node = 1; node <= numCellNodes; node++ ) {    
            vtkPiece.setNodeCoords(nodeNum, nodeCoords[node-1] );
            nodeNum += 1;
        }

        // Connectivity       
        for ( int i = 1; i <= numCellNodes; i++ ) {            
            nodes.at(i) = val++;
        }
        vtkPiece.setConnectivity(layer, nodes);
        
        // Offset
        offset += numCellNodes;
        vtkPiece.setOffset(layer, offset);

        // Cell types
        vtkPiece.setCellType(layer, 26); // Quadratic wedge
    }


    // Export nodal variables from primary fields        
    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numTotalNodes);

    std::vector<FloatArray> updatedNodeCoords;
    FloatArray u(3), damage(1);
    std::vector<FloatArray> values;
    FloatArray dVecTotal, dVecLayer, damageArray;
	this->computeDamageUnknowns(dVecTotal, VM_Total, tStep);		// vector with all damage nodal values for this element
    
    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);
        nodeNum = 1;
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            
            if ( type == DisplacementVector ) { // compute displacement as u = x - X
                this->giveFictiousNodeCoordsForExport(nodeCoords, layer);
                this->giveFictiousUpdatedNodeCoordsForExport(updatedNodeCoords, layer, tStep);
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    u = updatedNodeCoords[j-1];
                    u.subtract(nodeCoords[j-1]);
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, u);
                    nodeNum += 1;        
                }
            } else if ( type == ScalarDamage ) {
                // compute damage in layer
                // set same damage in the nodes above and below the nodes on the midsurface
                IntArray indx_d = computeDamageIndexArray(layer);

	            dVecLayer.beSubArrayOf(dVecTotal, indx_d);                      // vector of damage variables for a given layer  
                damageArray.setValues(15, dVecLayer.at(1), dVecLayer.at(2), dVecLayer.at(3), dVecLayer.at(1), dVecLayer.at(2), dVecLayer.at(3),
                                          dVecLayer.at(4), dVecLayer.at(5), dVecLayer.at(6), dVecLayer.at(4), dVecLayer.at(5), dVecLayer.at(6),
                                          dVecLayer.at(1), dVecLayer.at(2), dVecLayer.at(3) );
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    damage.at(1) = damageArray.at(j);
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, damage);
                    nodeNum += 1;        
                }
            } else {
                NodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values[j-1]);
                    nodeNum += 1;
                }
            }
        }
    }

    // Export nodal variables from internal fields
    
    vtkPiece.setNumberOfInternalVarsToExport( internalVarsToExport.giveSize(), numTotalNodes );
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        //this->recoverShearStress(tStep);
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            recoverValuesFromIP(values, layer, type, tStep);        
            for ( int j = 1; j <= numCellNodes; j++ ) {
                vtkPiece.setInternalVarInNode( fieldNum, nodeNum, values[j-1] );
                //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);          
                nodeNum += 1;        
            }                                
        }  
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);;
      
        for ( int layer = 1; layer <= numCells; layer++ ) {     
            IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
            VTKXMLExportModule::computeIPAverage(average, iRuleL, this, type, tStep);
            
            if ( average.giveSize() == 6 ) {
                vtkPiece.setCellVar(i, layer, convV6ToV9Stress(average) );
            } else {
                vtkPiece.setCellVar(i, layer, average );
            }

        }

    }


    

}

#if 0
void 
Shell7BasePhFi :: recoverValuesFromIP(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{
    // recover nodal values by coosing the ip closest to the node

    //FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );

    // composite element interpolator
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);

    int numNodes = localNodeCoords.giveNumberOfColumns();
    recoveredValues.resize(numNodes);
    
    IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
    IntegrationPoint *ip;

    // Find closest ip to the nodes
    IntArray closestIPArray(numNodes);
    FloatArray nodeCoords, ipCoords, ipValues;

    for ( int i = 1; i <= numNodes; i++ ) {
        nodeCoords.beColumnOf(localNodeCoords, i);
        double distOld = 3.0; // should not be larger
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            ip = iRule->getIntegrationPoint(j);
            ipCoords = *ip->giveCoordinates();
            double dist = nodeCoords.distance(ipCoords);
            if ( dist < distOld ) {
                closestIPArray.at(i) = j;
                distOld = dist;
            }
        }
    }

   InternalStateValueType valueType =  giveInternalStateValueType(type);

    // recover ip values
    for ( int i = 1; i <= numNodes; i++ ) {
        ip = iRule->getIntegrationPoint( closestIPArray.at(i) );
        this->giveIPValue(ipValues, ip, type, tStep);
        if ( valueType == ISVT_TENSOR_S3 ) {
            recoveredValues[i-1].resize(9);
            recoveredValues[i-1] = convV6ToV9Stress(ipValues);
        } else if ( ipValues.giveSize() == 0 && type == IST_AbaqusStateVector) {
            recoveredValues[i-1].resize(23);
            recoveredValues[i-1].zero();
        } else if ( ipValues.giveSize() == 0 ) {
            recoveredValues[i-1].resize(giveInternalStateTypeSize(valueType));
            recoveredValues[i-1].zero();

        } else {
            recoveredValues[i-1] = ipValues;
        }
    }

}


void 
Shell7BasePhFi :: recoverShearStress(TimeStep *tStep)
{
    // Recover shear stresses at ip by numerical integration of the momentum balance through the thickness
    std::vector<FloatArray> recoveredValues;
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    IntegrationRule *iRuleThickness = specialIntegrationRulesArray[ 0 ];
    FloatArray dS, Sold;
    FloatMatrix B, Smat(2,6); // 2 stress components * num of in plane ip ///@todo generalize
    Smat.zero();
    FloatArray Tcon(6), Trec(6);  Tcon.zero(); Trec.zero();
    
     for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        this->recoverValuesFromIP(recoveredValues, layer, IST_StressTensor, tStep);
        //this->ZZNodalRecoveryMI_recoverValues(recoveredValues, layer, IST_StressTensor, tStep);
        double thickness = this->layeredCS->giveLayerThickness(layer);

        //set up vector of stresses in the ip's = [S1_xx, S1_yy, S1_xy, ..., Sn_xx, Sn_yy, Sn_xy]
        int numNodes = 15;
        FloatArray aS(numNodes*3); 
        for ( int j = 1, pos = 0; j <= numNodes; j++, pos+=3 ) {
            aS.at(pos + 1) = recoveredValues[j-1].at(1);   // S_xx
            aS.at(pos + 2) = recoveredValues[j-1].at(2);   // S_yy
            aS.at(pos + 3) = recoveredValues[j-1].at(6);   // S_xy
        }
        int numInPlaneIP = 6;

        for ( int i = 0; i < iRuleThickness->giveNumberOfIntegrationPoints(); i++ ) { 
            double  dz = thickness * iRuleThickness->getIntegrationPoint(i)->giveWeight();
            
            for ( int j = 0; j < numInPlaneIP; j++ ) { 

                int point = i*numInPlaneIP + j; // integration point number
                GaussPoint *gp = iRuleL->getIntegrationPoint(point);

                this->computeBmatrixForStressRecAt(*gp->giveCoordinates(), B, layer);
                dS.beProductOf(B,aS*(-dz)); // stress increment

                StructuralMaterialStatus* status = dynamic_cast< StructuralMaterialStatus* > ( gp->giveMaterialStatus() );
                Sold = status->giveStressVector();
                
                Smat.at(1,j+1) += dS.at(1); // add increment from each level
                Smat.at(2,j+1) += dS.at(2);

                //Tcon.at(j+1) += Sold.at(5)*dz;

                // Replace old stresses with  - this should probably not be done as it may affect the convergence in a nonlinear case
                Sold.at(5) = Smat.at(1,j+1); // S_xz
                Sold.at(4) = Smat.at(2,j+1); // S_yz


                status->letStressVectorBe(Sold);
                //Trec.at(j+1) += Sold.at(5)*dz;
            }
        }


    }

}


void
Shell7BasePhFi :: computeBmatrixForStressRecAt(FloatArray &lcoords, FloatMatrix &answer, int layer)
{
    // Returns the  special matrix {B} of the receiver, evaluated at aGaussPoint. Such that
    // B*a = [dS_xx/dx + dS_xy/dy, dS_yx/dx + dS_yy/dy ]^T, where a is the vector of in plane 
    // stresses [S_xx, S_yy, S_xy]
 
    // set up virtual cell geometry for an qwedge
    const int numNodes = 15;
    std::vector<FloatArray> nodes;
    giveFictiousNodeCoordsForExport(nodes, layer);

    int VTKWedge2EL [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    FloatArray *coords[numNodes];


    for ( int i = 1; i <= numNodes; i++ ) {
        int pos = VTKWedge2EL[ i-1 ];
        coords[ i - 1 ] = &nodes[ pos - 1];
        
    }
    
    FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );
    FloatMatrix dNdx;
    interpol->evaldNdx( dNdx, lcoords, FEIVertexListGeometryWrapper(numNodes, (const FloatArray **)coords ) );
    
    /*    
     * 1 [d/dx  0   d/dy
     * 1   0   d/dy d/dx]
     */
    int ndofs = numNodes*3;
    answer.resize(2, ndofs);
    for ( int i = 1, j = 0; i <= numNodes; i++, j += 3 ) {
        answer.at(1, j + 1) = dNdx.at(i, 1);
        answer.at(1, j + 3) = dNdx.at(i, 2);
        answer.at(2, j + 2) = dNdx.at(i, 2);
        answer.at(2, j + 3) = dNdx.at(i, 1);
    }
    
}





void 
Shell7BasePhFi :: giveFictiousNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}


void 
Shell7BasePhFi :: giveFictiousCZNodeCoordsForExport(std::vector<FloatArray> &nodes, int interface)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForCZExport.giveLocalNodeCoords(localNodeCoords);
    
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        localCoords.at(3) = 1.0;
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, interface, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}

void 
Shell7BasePhFi :: giveFictiousUpdatedNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, TimeStep *tStep)
{
    // compute fictious node coords

    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }
}


//void
//Shell7BasePhFi :: giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords) {
//    // Local coords for a quadratic wedge element (VTK cell type 26)
//    double z = 0.999;
//    nodeLocalXi1Coords.setValues(15, 1., 0., 0., 1., 0., 0., .5, 0., .5, .5, 0., .5, 1., 0., 0.);      
//    nodeLocalXi2Coords.setValues(15, 0., 1., 0., 0., 1., 0., .5, .5, 0., .5, .5, 0., 0., 1., 0.);
//    nodeLocalXi3Coords.setValues(15, -z, -z, -z,  z,  z,  z, -z, -z, -z,  z,  z,  z, 0., 0., 0.);
//}

#endif



// Misc functions



} // end namespace oofem
