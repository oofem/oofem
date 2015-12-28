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

#include "prescribedmean.h"
#include "classfactory.h"
#include "masterdof.h"
#include "domain.h"
#include "feinterpol.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "function.h"

namespace oofem
{

double PrescribedMean :: domainSize;


REGISTER_BoundaryCondition(PrescribedMean);

IRResultType
PrescribedMean :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    GeneralBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, c, _IFT_PrescribedMean_Mean);
    IR_GIVE_FIELD(ir, dofid, _IFT_PrescribedMean_DofID);
    IR_GIVE_FIELD(ir, set, _IFT_GeneralBoundaryCondition_set);

    elementEdges = false;
    IR_GIVE_OPTIONAL_FIELD(ir, elementEdges, _IFT_PrescribedMean_Edge);

    int newdofid = this->domain->giveNextFreeDofID();
    lambdaIDs.clear();
    lambdaIDs.followedBy(newdofid);
    lambdaDman->appendDof( new MasterDof( lambdaDman, ( DofIDItem )newdofid ));

    domainSize=-1.;

    return IRRT_OK;

}

void
PrescribedMean :: assemble(SparseMtrx &answer, TimeStep *tStep, CharType type,
                           const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{

    if ( type != TangentStiffnessMatrix && type != StiffnessMatrix ) {
        return;
    }

    computeDomainSize();

    IntArray c_loc, r_loc;
    lambdaDman->giveLocationArray(lambdaIDs, r_loc, r_s);
    lambdaDman->giveLocationArray(lambdaIDs, c_loc, c_s);

    for ( int i = 1; i <= elements.giveSize(); i++ ) {
        int elementID = elements.at(i);
        Element *thisElement = this->giveDomain()->giveElement(elementID);
        FEInterpolation *interpolator = thisElement->giveInterpolation(DofIDItem(dofid));

        IntegrationRule *iRule = (elementEdges) ? (interpolator->giveBoundaryIntegrationRule(3, sides.at(i))) :
                                                  (interpolator->giveIntegrationRule(3));

        for ( GaussPoint * gp: * iRule ) {
            FloatArray lcoords = gp->giveNaturalCoordinates();
            FloatArray N; //, a;
            FloatMatrix temp, tempT;
            double detJ = 0.0;
            IntArray boundaryNodes, dofids={(DofIDItem) this->dofid}, r_Sideloc, c_Sideloc;

            if (elementEdges) {
                // Compute boundary integral
                interpolator->boundaryGiveNodes( boundaryNodes, sides.at(i) );
                interpolator->boundaryEvalN(N, sides.at(i), lcoords, FEIElementGeometryWrapper(thisElement));
                detJ = fabs ( interpolator->boundaryGiveTransformationJacobian(sides.at(i), lcoords, FEIElementGeometryWrapper(thisElement)) );
                // Retrieve locations for dofs on boundary
                thisElement->giveBoundaryLocationArray(r_Sideloc, boundaryNodes, dofids, r_s);
                thisElement->giveBoundaryLocationArray(c_Sideloc, boundaryNodes, dofids, c_s);
            } else {
                interpolator->evalN(N, lcoords, FEIElementGeometryWrapper(thisElement));
                detJ = fabs ( interpolator->giveTransformationJacobian(lcoords, FEIElementGeometryWrapper(thisElement) ) );
                IntArray DofIDStemp, rloc, cloc;

                thisElement->giveLocationArray(rloc, r_s, &DofIDStemp);
                thisElement->giveLocationArray(cloc, c_s, &DofIDStemp);

                r_Sideloc.clear();
                c_Sideloc.clear();
                for (int j=1; j<=DofIDStemp.giveSize(); j++) {
                    if ( DofIDStemp.at(j) == dofids.at(1) ) {
                        r_Sideloc.followedBy(rloc.at(j));
                        c_Sideloc.followedBy(cloc.at(j));
                    }
                }
            }

            // delta p part:
            temp = N*detJ*gp->giveWeight()*(1.0/domainSize);
            tempT.beTranspositionOf(temp);

            answer.assemble(r_Sideloc, c_loc, temp);
            answer.assemble(r_loc, c_Sideloc, tempT);
        }

        delete iRule;

    }

}

void
PrescribedMean :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                 CharType type, ValueModeType mode,
                                 const UnknownNumberingScheme &s, FloatArray *eNorm)
{

    if ( type == InternalForcesVector ) {
        giveInternalForcesVector(answer, tStep, type, mode, s);
    } else if ( type == ExternalForcesVector ) {
        giveExternalForcesVector(answer, tStep, type, mode, s);
    }
}

void
PrescribedMean :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep,
                                           CharType type, ValueModeType mode,
                                           const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    computeDomainSize();

    // Fetch unknowns of this boundary condition
    IntArray lambdaLoc;
    FloatArray lambda;
    lambdaDman->giveUnknownVector(lambda, lambdaIDs, mode, tStep);
    lambdaDman->giveLocationArray(lambdaIDs, lambdaLoc, s);

    for ( int i = 1; i <= elements.giveSize(); i++ ) {
        int elementID = elements.at(i);
        Element *thisElement = this->giveDomain()->giveElement(elementID);
        FEInterpolation *interpolator = thisElement->giveInterpolation(DofIDItem(dofid));

        IntegrationRule *iRule = (elementEdges) ? (interpolator->giveBoundaryIntegrationRule(3, sides.at(i))) :
                                                  (interpolator->giveIntegrationRule(3));

        for ( auto &gp: * iRule ) {
            FloatArray lcoords = gp->giveNaturalCoordinates();
            FloatArray a, N, pressureEqns, lambdaEqns;
            IntArray boundaryNodes, dofids={(DofIDItem) this->dofid}, locationArray;
            double detJ = 0.0;

            if (elementEdges) {
                // Compute integral
                interpolator->boundaryGiveNodes( boundaryNodes, sides.at(i) );
                thisElement->computeBoundaryVectorOf(boundaryNodes, dofids, VM_Total, tStep, a);
                interpolator->boundaryEvalN(N, sides.at(i), lcoords, FEIElementGeometryWrapper(thisElement));
                detJ = fabs ( interpolator->boundaryGiveTransformationJacobian(sides.at(i), lcoords, FEIElementGeometryWrapper(thisElement)) );

                // Retrieve locations for dofs with dofids
                thisElement->giveBoundaryLocationArray(locationArray, boundaryNodes, dofids, s);
            } else {
                thisElement->computeVectorOf(dofids, VM_Total, tStep, a);
                interpolator->evalN(N, lcoords, FEIElementGeometryWrapper(thisElement));
                detJ = fabs ( interpolator->giveTransformationJacobian(lcoords, FEIElementGeometryWrapper(thisElement)));

                IntArray DofIDStemp, loc;

                thisElement->giveLocationArray(loc, s, &DofIDStemp);

                locationArray.clear();
                for (int j=1; j<=DofIDStemp.giveSize(); j++) {
                    if ( DofIDStemp.at(j) == dofids.at(1) ) {
                        locationArray.followedBy(loc.at(j));
                    }
                }
            }

            // delta p part:
            pressureEqns = N*detJ*gp->giveWeight()*lambda.at(1)*(1.0/domainSize);

            // delta lambda part
            lambdaEqns.resize(1);
            lambdaEqns.at(1) = N.dotProduct(a);
            lambdaEqns.times(detJ*gp->giveWeight()*1.0/domainSize);
            lambdaEqns.at(1) = lambdaEqns.at(1);


            // delta p part
            answer.assemble(pressureEqns, locationArray);

            // delta lambda part
            answer.assemble(lambdaEqns, lambdaLoc);
        }
        delete iRule;
    }

}

void
PrescribedMean :: giveExternalForcesVector(FloatArray &answer, TimeStep *tStep,
                                           CharType type, ValueModeType mode,
                                           const UnknownNumberingScheme &s)
{
    computeDomainSize();

    FloatArray temp;
    IntArray lambdaLoc;

    temp.resize(1);
    temp.at(1) = c;

    lambdaDman->giveLocationArray(lambdaIDs, lambdaLoc, s);
    answer.assemble(temp, lambdaLoc);

    // Finally, compute value of loadtimefunction
    double factor;
    factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer.times(factor);
}

void
PrescribedMean :: computeDomainSize()
{
    if ( domainSize > 0.0 ) return;


    if ( elementEdges ) {
        IntArray setList = ((GeneralBoundaryCondition *)this)->giveDomain()->giveSet(set)->giveBoundaryList();

        elements.resize(setList.giveSize() / 2);
        sides.resize(setList.giveSize() / 2);

        for (int i = 1; i <= setList.giveSize(); i += 2) {
            elements.at(i/2+1) = setList.at(i);
            sides.at(i/2+1) = setList.at(i+1);
        }
    } else {
        IntArray setList = ((GeneralBoundaryCondition *)this)->giveDomain()->giveSet(set)->giveElementList();
        elements = setList;
    }

    domainSize = 0.0;

    for ( int i = 1; i <= elements.giveSize(); i++ ) {
        int elementID = elements.at(i);
        Element *thisElement = this->giveDomain()->giveElement(elementID);
        FEInterpolation *interpolator = thisElement->giveInterpolation(DofIDItem(dofid));

        IntegrationRule *iRule;

        if ( elementEdges ) {
            iRule = interpolator->giveBoundaryIntegrationRule(3, sides.at(i));
        } else {
            iRule = interpolator->giveIntegrationRule(3);
        }

        for ( GaussPoint * gp: * iRule ) {
            FloatArray lcoords = gp->giveNaturalCoordinates();

            double detJ;
            if ( elementEdges ) {
                detJ = fabs ( interpolator->boundaryGiveTransformationJacobian(sides.at(i), lcoords, FEIElementGeometryWrapper(thisElement)) );
            } else {
                detJ = fabs ( interpolator->giveTransformationJacobian(lcoords, FEIElementGeometryWrapper(thisElement)) );
            }
            domainSize = domainSize + detJ*gp->giveWeight();
        }

        delete iRule;
    }

    printf("%f\n", domainSize);

}

}
