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

#include "plprincipalstrain.h"
#include "xfem/propagationlaw.h"
#include "xfem/tipinfo.h"
#include "classfactory.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "spatiallocalizer.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "xfem/enrichmentitem.h"
#include "feinterpol.h"
#include "xfem/xfemmanager.h"

#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"

#include "xfem/XFEMDebugTools.h"

namespace oofem {
REGISTER_PropagationLaw(PLPrincipalStrain)

PLPrincipalStrain::PLPrincipalStrain():
mRadius(0.0),mIncrementLength(0.0),mStrainThreshold(0.0), mUseRadialBasisFunc(false)
{

}

PLPrincipalStrain::~PLPrincipalStrain() {

}

void PLPrincipalStrain :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, mRadius,                          _IFT_PLPrincipalStrain_Radius);
    IR_GIVE_FIELD(ir, mIncrementLength,         _IFT_PLPrincipalStrain_IncLength);
    IR_GIVE_FIELD(ir, mStrainThreshold, _IFT_PLPrincipalStrain_StrainThreshold);

    int useRadialBasisFunc = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, useRadialBasisFunc, _IFT_PLPrincipalStrain_RadialBasisFunc);
    if ( useRadialBasisFunc == 1 ) {
        mUseRadialBasisFunc = true;
    }
}

void PLPrincipalStrain :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mRadius,                             _IFT_PLPrincipalStrain_Radius);
    input.setField(mIncrementLength,            _IFT_PLPrincipalStrain_IncLength);
    input.setField(mStrainThreshold,        _IFT_PLPrincipalStrain_StrainThreshold);

    if ( mUseRadialBasisFunc ) {
        input.setField(1,       _IFT_PLPrincipalStrain_RadialBasisFunc);
    }
}


bool PLPrincipalStrain :: propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp)
{
    if ( !iEnrFront.propagationIsAllowed() ) {
        return false;
    }

    // Fetch crack tip data
    const TipInfo &tipInfo = iEnrFront.giveTipInfo();

    SpatialLocalizer *localizer = iDomain.giveSpatialLocalizer();


    const FloatArray &xT    = tipInfo.mGlobalCoord;
    const FloatArray &t     = tipInfo.mTangDir;
//    const FloatArray &n     = tipInfo.mNormalDir;

    // It is meaningless to propagate a tip that is not inside any element
    Element *el = localizer->giveElementContainingPoint(tipInfo.mGlobalCoord);
    if ( el != nullptr ) {



    	FloatArray x(xT);
    	x.add(mRadius, t);
//    	circPoints.push_back(x);


        std :: vector< double >sigTTArray, sigRTArray;

        FloatArray strainVec;

        if ( mUseRadialBasisFunc ) {
        	// Interpolate strain with radial basis functions

        	// Choose a cut-off length l:
        	// take the distance between two nodes in the element containing the
        	// crack tip multiplied by a constant factor.
        	// ( This choice implies that we hope that the element has reasonable
        	// aspect ratio.)
        	const auto &x1 = el->giveDofManager(1)->giveCoordinates();
        	const auto &x2 = el->giveDofManager(2)->giveCoordinates();
            const double l = 1.0 * distance(x1, x2);

        	// Use the octree to get all elements that have
        	// at least one Gauss point in a certain region around the tip.
        	const double searchRadius = 3.0 * l;
        	IntArray elIndices;
        	localizer->giveAllElementsWithIpWithinBox(elIndices, x, searchRadius);


        	// Loop over the elements and Gauss points obtained.
        	// Evaluate the interpolation.
        	FloatArray sumQiWiVi;
        	double sumWiVi = 0.0;
        	for ( int elIndex: elIndices ) {
        		Element *gpEl = iDomain.giveElement(elIndex);

        		for ( auto &gp_i: *gpEl->giveDefaultIntegrationRulePtr() ) {
        			////////////////////////////////////////
        			// Compute global gp coordinates
        			FloatArray N;
        			FEInterpolation *interp = gpEl->giveInterpolation();
        			interp->evalN( N, gp_i->giveNaturalCoordinates(), FEIElementGeometryWrapper(gpEl) );


        			// Compute global coordinates of Gauss point
        			FloatArray globalCoord(2);
        			globalCoord.zero();

        			for ( int i = 1; i <= gpEl->giveNumberOfDofManagers(); i++ ) {
        				DofManager *dMan = gpEl->giveDofManager(i);
        				globalCoord.at(1) += N.at(i) * dMan->giveCoordinate(1);
        				globalCoord.at(2) += N.at(i) * dMan->giveCoordinate(2);
        			}


        			////////////////////////////////////////
        			// Compute weight of kernel function

        			FloatArray tipToGP;
        			tipToGP.beDifferenceOf(globalCoord, xT);
        			bool inFrontOfCrack = true;
        			if ( tipToGP.dotProduct(t) < 0.0 ) {
        				inFrontOfCrack = false;
        			}

        			double r = distance(x, globalCoord);

        			if ( r < l && inFrontOfCrack ) {
        				double w = ( ( l - r ) / ( pow(2.0 * M_PI, 1.5) * pow(l, 3) ) ) * exp( -0.5 * pow(r, 2) / pow(l, 2) );

        				// Compute gp volume
        				double V = gpEl->computeVolumeAround(gp_i);

        				// Get stress
        				StructuralMaterialStatus *ms = dynamic_cast< StructuralMaterialStatus * >( gp_i->giveMaterialStatus() );
        				if ( ms == nullptr ) {
        					OOFEM_ERROR("failed to fetch MaterialStatus.");
        				}

        				FloatArray strainVecGP = ms->giveStrainVector();

        				if ( sumQiWiVi.giveSize() != strainVecGP.giveSize() ) {
        					sumQiWiVi.resize( strainVecGP.giveSize() );
        					sumQiWiVi.zero();
        				}

        				// Add to numerator
        				sumQiWiVi.add(w * V, strainVecGP);

        				// Add to denominator
        				sumWiVi += w * V;
        			}
        		}
        	}


        	if ( fabs(sumWiVi) > 1.0e-12 ) {
        		strainVec.beScaled(1.0 / sumWiVi, sumQiWiVi);
        	} else {
        		// Take strain from closest Gauss point
        		int region = 1;
        		bool useCZGP = false;
        		GaussPoint &gp = * ( localizer->giveClosestIP(x, region, useCZGP) );


        		// Compute strain
        		StructuralMaterialStatus *ms = dynamic_cast< StructuralMaterialStatus * >( gp.giveMaterialStatus() );
        		if ( ms == nullptr ) {
        			OOFEM_ERROR("failed to fetch MaterialStatus.");
        		}

        		strainVec = ms->giveStrainVector();
        	}
        } else {
        	// Take stress from closest Gauss point
        	int region = 1;
        	bool useCZGP = false;
        	GaussPoint &gp = * ( localizer->giveClosestIP(x, region, useCZGP) );


        	// Compute stresses
        	StructuralMaterialStatus *ms = dynamic_cast< StructuralMaterialStatus * >( gp.giveMaterialStatus() );
        	if ( ms == nullptr ) {
        		OOFEM_ERROR("failed to fetch MaterialStatus.");
        	}

        	strainVec = ms->giveStrainVector();
        }

        // Compute principal strain
		FloatArray principalVals;
		FloatMatrix principalDirs;
		StructuralMaterial::computePrincipalValDir(principalVals, principalDirs, strainVec, principal_strain);



        // Compare with threshold
        printf("Max principal strain: %e\n", principalVals[0]);
        if ( principalVals[0] > mStrainThreshold ) {

			FloatArray propNormal;
			propNormal.beColumnOf(principalDirs, 1);

			FloatArray propTangent = {-propNormal(1), propNormal(0)};

			if( propTangent.dotProduct(t) < 0.0 ) {
				propTangent.times(-1.0);
			}

            // Fill up struct
            oTipProp.mTipIndex = tipInfo.mTipIndex;
            oTipProp.mPropagationDir = propTangent;
            oTipProp.mPropagationLength = mIncrementLength;

            return true;
        }

    } // Tip is inside an element.


    return false;
}

} // end namespace oofem

