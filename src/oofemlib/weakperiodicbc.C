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

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <memory>

#include "activebc.h"
#include "weakperiodicbc.h"
#include "inputrecord.h"
#include "element.h"
#include "node.h"
#include "masterdof.h"
#include "sparsemtrx.h"
#include "gausspoint.h"
#include "integrationrule.h"
#include "mathfem.h"
#include "fei2dtrlin.h"
#include "fei2dtrquad.h"
#include "classfactory.h"
#include "set.h"
#include "function.h"

#include "timestep.h"
// #include "../sm/Elements/tet21ghostsolid.h"
#include "../sm/Elements/nlstructuralelement.h"

namespace oofem {
REGISTER_BoundaryCondition(WeakPeriodicBoundaryCondition);

WeakPeriodicBoundaryCondition :: WeakPeriodicBoundaryCondition(int n, Domain *d) : ActiveBoundaryCondition(n, d), gammaDman( new Node(0, this->domain) )
{
    useBasisType = monomial;
    doUpdateSminmax = true;
}

WeakPeriodicBoundaryCondition :: ~WeakPeriodicBoundaryCondition()
{
}

IRResultType
WeakPeriodicBoundaryCondition :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    result = ActiveBoundaryCondition :: initializeFrom(ir);     ///@todo Carl, remove this line and use elementsidespositive/negative instead.
    if ( result != IRRT_OK ) {
        return result;
    }

    orderOfPolygon = 2;
    IR_GIVE_OPTIONAL_FIELD(ir, orderOfPolygon, _IFT_WeakPeriodicBoundaryCondition_order);

    int t = ( int ) monomial;
    IR_GIVE_OPTIONAL_FIELD(ir, t, _IFT_WeakPeriodicBoundaryCondition_descritizationType);
    useBasisType = ( basisType ) t;      // Fourierseries by default

    //dofids = P_f;        // Pressure as default
    IR_GIVE_OPTIONAL_FIELD(ir, dofids, _IFT_WeakPeriodicBoundaryCondition_dofids);
    ndofids = dofids.giveSize();

    ngp = -1;        // Pressure as default
    IR_GIVE_OPTIONAL_FIELD(ir, ngp, _IFT_WeakPeriodicBoundaryCondition_ngp);
    if ( ngp != -1 ) {
        OOFEM_WARNING("ngp isn't being used anymore! see how the interpolator constructs the integration rule automatically.");
        return IRRT_BAD_FORMAT;
    }


    nlgeo = false;
    IR_GIVE_OPTIONAL_FIELD(ir, nlgeo, _IFT_WeakPeriodicBoundaryCondition_nlgeo );


    g.resize(domain->giveNumberOfSpatialDimensions());
    g.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, g, _IFT_WeakPeriodicBoundaryCondition_gradient);

    IntArray temp;

    posSet = -1;
    negSet = -1;

    IR_GIVE_OPTIONAL_FIELD(ir, posSet, _IFT_WeakPeriodicBoundaryCondition_elementSidesPositiveSet);
    IR_GIVE_OPTIONAL_FIELD(ir, negSet, _IFT_WeakPeriodicBoundaryCondition_elementSidesNegativeSet);

    if ( posSet == -1 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, temp, _IFT_WeakPeriodicBoundaryCondition_elementSidesPositive);
        for ( int i = 0; i < temp.giveSize() / 2; i++ ) {
            side [ 0 ].push_back( temp.at(2 * i + 1) );
            element [ 0 ].push_back( temp.at(2 * i + 2) );
        }
    }

    if ( negSet == -1 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, temp, _IFT_WeakPeriodicBoundaryCondition_elementSidesNegative);
        for ( int i = 0; i < temp.giveSize() / 2; i++ ) {
            side [ 1 ].push_back( temp.at(2 * i + 1) );
            element [ 1 ].push_back( temp.at(2 * i + 2) );
        }
    }

    if ( this->domain->giveNumberOfSpatialDimensions() == 2 ) {
        ndof = (orderOfPolygon + 1) * dofids.giveSize();
        tcount = (orderOfPolygon + 1);
    } else if ( this->domain->giveNumberOfSpatialDimensions() == 3 ) {
        ndof = 1;
        for ( int i = 1; i <= orderOfPolygon; i++ ) {
            ndof = ndof + ( i + 1 );
        }
        tcount = ndof;
        ndof = ndof * ndofids;
    }

    // Create dofs for coefficients
    bcID = this->giveNumber();
    gamma_ids.clear();
    for ( int i = 0; i < ndof; i++ ) {
        int dofid = this->domain->giveNextFreeDofID();
        gamma_ids.followedBy(dofid);
        gammaDman->appendDof( new MasterDof( gammaDman.get(), ( DofIDItem )dofid ) );
    }

    return IRRT_OK;
}

void WeakPeriodicBoundaryCondition :: computeOrthogonalBasis()
{
    /* gsMatrix contains the coefficients for the orthogonal basis. The row represents the basis and the column the coefficients. */
    gsMatrix.resize(ndof, ndof);
    gsMatrix.zero();
    gsMatrix.at(1, 1) = 1;

    for ( int i = 2; i <= ndof; i++ ) { // Need ndof base functions. i indicates the row of gsMatrix, ie which polynomial is the current.
        gsMatrix.at(i, i) = 1;       // Copy from V

        // remove projection of v_i on all previous base functions
        for ( int j = 1; j < i; j++ ) {
            FloatArray uTemp;
            uTemp.resize(ndof);

            for ( int k = 1; k <= ndof; uTemp.at(k) = gsMatrix.at(j, k), k++ ) {
                ;
            }

            double thisValue = computeProjectionCoefficient(i, j);
            uTemp.times(-thisValue);

            for ( int k = 1; k <= ndof; gsMatrix.at(i, k) = gsMatrix.at(i, k) + uTemp.at(k), k++ ) {
                ;
            }

        }
    }
}

double WeakPeriodicBoundaryCondition :: computeProjectionCoefficient(int vIndex, int uIndex)
{
    /* Computes <v, u>/<u, u> where vIndex is the term in the polynomial and uIndex is the number of the base v being projected on. */
    int thisSide = 0;

    double value = 0.0, nom = 0.0, denom = 0.0;

    int A, B, thisOrder;
    getExponents(vIndex, A, B);
    thisOrder = A + B;
    getExponents(uIndex, A, B);
    thisOrder = thisOrder + A + B;

    // Loop over all elements
    for ( size_t ielement = 0; ielement < element [ thisSide ].size(); ielement++ ) {
        // Compute <v, u_i>/<u_i, u_i> on this element and store in nom/denom

        Element *thisElement = this->domain->giveElement( element [ thisSide ].at(ielement) );
        FEInterpolation *geoInterpolation = thisElement->giveInterpolation();

        std :: unique_ptr< IntegrationRule >iRule(geoInterpolation->giveBoundaryIntegrationRule(thisOrder, side [ thisSide ].at(ielement)));

        for ( GaussPoint *gp: *iRule ) {

            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            FloatArray gcoords;

            geoInterpolation->boundaryLocal2Global( gcoords, side [ thisSide ].at(ielement), lcoords, FEIElementGeometryWrapper(thisElement) );
            double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( side [ thisSide ].at(ielement), lcoords, FEIElementGeometryWrapper(thisElement) ) );

            int a, b;
            getExponents(vIndex, a, b);
            double vVal = pow(gcoords.at( surfaceIndexes.at(1) ), a) * pow(gcoords.at( surfaceIndexes.at(2) ), b);

            double uValue = 0.0;
            for ( int k = 1; k < ndof; k++ ) {
                int c, d;
                getExponents(k, c, d);
                uValue = uValue + gsMatrix.at(uIndex, k) * pow(gcoords.at( surfaceIndexes.at(1) ), c) * pow(gcoords.at( surfaceIndexes.at(2) ), d);
            }

            nom = nom + vVal *uValue *detJ *gp->giveWeight();
            denom = denom + uValue *uValue *detJ *gp->giveWeight();
        }
    }

    value = nom / denom;

    return value;
}

void WeakPeriodicBoundaryCondition :: giveEdgeNormal(FloatArray &answer, int element, int side)
{
    FloatArray xi;

    if ( this->domain->giveNumberOfSpatialDimensions() == 3 ) {
        xi.resize(2);
        xi(0) = 0.25;
        xi(1) = 0.25;
    } else {
        xi.resize(1);
        xi(0) = 0.5;
    }

    Element *thisElement = this->domain->giveElement(element);
    FEInterpolation *interpolation = thisElement->giveInterpolation( ( DofIDItem ) dofids(0) );

    interpolation->boundaryEvalNormal( answer, side, xi, FEIElementGeometryWrapper(thisElement) );
}

void WeakPeriodicBoundaryCondition :: updateDirection()
{
    // Check orientation for s
    FloatArray normal;

    if ( this->domain->giveNumberOfSpatialDimensions() == 2 ) {
        surfaceIndexes.resize(1);
        smin.resize(1);
        smax.resize(1);
    } else {
        surfaceIndexes.resize(2);
        smin.resize(2);
        smax.resize(2);
    }

    giveEdgeNormal( normal, element [ 0 ].at(0), side [ 0 ].at(0) );

    if ( fabs( normal.at(1) ) > 0.99999 ) {              // Normal points in X direction
        direction = 1;
        if ( this->domain->giveNumberOfSpatialDimensions() == 2 ) {
            surfaceIndexes.at(1) = 2;
        } else {
            surfaceIndexes.at(1) = 2;
            surfaceIndexes.at(2) = 3;
        }
    } else if ( fabs( normal.at(2) ) > 0.99999 ) {        // Normal points in Y direction
        direction = 2;
        if ( this->domain->giveNumberOfSpatialDimensions() == 2 ) {
            surfaceIndexes.at(1) = 1;
        } else {
            surfaceIndexes.at(1) = 1;
            surfaceIndexes.at(2) = 3;
        }
    } else if ( fabs( normal.at(3) ) > 0.99 ) {         // Normal points in Z direction
        direction = 3;
        if ( this->domain->giveNumberOfSpatialDimensions() == 2 ) {
            OOFEM_ERROR("3 dimensioal normal in a 2 dimensional problem.");
        } else {
            surfaceIndexes.at(1) = 1;
            surfaceIndexes.at(2) = 2;
        }
    } else {
        normal.printYourself();
        Element *thisElement = this->giveDomain()->giveElement( element [ 0 ].at(0) );
        OOFEM_ERROR("Only surfaces with normal in x, y or z direction supported. (el=%d, side=%d) \n", thisElement->giveLabel(), side [ 0 ].at(0) );
    }
}

void WeakPeriodicBoundaryCondition :: updateSminmax()
{
    if ( doUpdateSminmax ) {
        // If sets are used, now is the time to update lists of elements and sides since the sets are unknown in initializeFrom
        if ( posSet != -1 ) {
            IntArray posBoundary, negBoundary;

            posBoundary = this->giveDomain()->giveSet(posSet)->giveBoundaryList();
            for ( int i = 0; i < posBoundary.giveSize() / 2; i++ ) {
                side [ 0 ].push_back( posBoundary.at(2 * i + 2) );
                element [ 0 ].push_back( posBoundary.at(2 * i + 1) );
            }

            negBoundary = this->giveDomain()->giveSet(negSet)->giveBoundaryList();
            for ( int i = 0; i < negBoundary.giveSize() / 2; i++ ) {
                side [ 1 ].push_back( negBoundary.at(2 * i + 2) );
                element [ 1 ].push_back( negBoundary.at(2 * i + 1) );
            }
        }

        // Determine which is the positive and which is the negative side
        FloatArray normal;
        giveEdgeNormal( normal, element [ 0 ].at(0), side [ 0 ].at(0) );

        double normalSum = -1;
        ( this->giveDomain()->giveNumberOfSpatialDimensions() <= 2 ) ? normalSum = normal.at(1) + normal.at(2) : normalSum = normal.at(1) + normal.at(2) + normal.at(3);

        if ( normalSum > -0.000001 ) {        // No support for 3D yet
            sideSign [ 0 ] = 1;
            sideSign [ 1 ] = -1;
        } else {
            sideSign [ 0 ] = -1;
            sideSign [ 1 ] = 1;
        }

        updateDirection();

        for ( int i = 1; i <= surfaceIndexes.giveSize(); i++ ) {
            smin.at(i) = this->domain->giveDofManager(1)->giveCoordinate( surfaceIndexes.at(i) );
            smax.at(i) = this->domain->giveDofManager(1)->giveCoordinate( surfaceIndexes.at(i) );

            for ( int j = 1; j <= this->domain->giveNumberOfDofManagers(); j++ ) {
                double sValue = this->domain->giveDofManager(j)->giveCoordinate( surfaceIndexes.at(i) );
                smin.at(i) = std :: min(smin.at(i), sValue);
                smax.at(i) = std :: max(smax.at(i), sValue);
            }
        }
        doUpdateSminmax = false;

        if ( this->useBasisType == legendre ) {
            computeOrthogonalBasis();
        }
    }
}

void WeakPeriodicBoundaryCondition :: addElementSide(int newElement, int newSide)
{

    FloatArray normalNew, normal0;
    int addToList = 0;

    if ( element [ 0 ].size() > 0 ) {
        // If there are elements in the list, compare normals in order to determine which list to store them in
        giveEdgeNormal(normalNew, newElement, newSide);
        giveEdgeNormal( normal0, element [ 0 ].at(0), side [ 0 ].at(0) );
        double d = sqrt( pow(normalNew.at(1) - normal0.at(1), 2) + pow(normalNew.at(2) - normal0.at(2), 2) );
        if ( fabs(d) < 0.001 ) {
            addToList = 0;
        } else {
            addToList = 1;
        }
    } else {
        // Otherwise, check the normal in order to decide upon which direction the parameter runs (x or y)
        giveEdgeNormal(normalNew, newElement, newSide);
        if ( fabs(fabs( normalNew.at(1) ) - 1.) < 0.0001 ) {                       // Normal point in x direction, thus set direction to y
            direction = 2;
        } else {
            direction = 1;
        }
    }

    element [ addToList ].push_back(newElement);
    side [ addToList ].push_back(newSide);
}

void
WeakPeriodicBoundaryCondition :: computeDeformationGradient(FloatMatrix &answer, Element *e, FloatArray *lcoord, TimeStep *tStep)
{

    FloatArray F, u;
    FloatMatrix dNdx, BH, Fmatrix, Finv;
    FEInterpolation *interpolation = e->giveInterpolation( ( DofIDItem ) dofids(0) );

    // Fetch displacements
    e->computeVectorOf({1, 2, 3}, VM_Total, tStep, u);

    // Compute dNdx on boundary.
    interpolation->evaldNdx(dNdx, *lcoord, FEIElementGeometryWrapper (e) );

    // Compute displcement gradient BH
    BH.resize(9, dNdx.giveNumberOfRows() * 3);
    BH.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        BH.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        BH.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        BH.at(3, 3 * i - 0) = dNdx.at(i, 3);     // dw/dz
        BH.at(4, 3 * i - 1) = dNdx.at(i, 3);     // dv/dz
        BH.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        BH.at(5, 3 * i - 2) = dNdx.at(i, 3);     // du/dz
        BH.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        BH.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        BH.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
    }

    // Finally, compute deformation gradient F=BH*u+I
    F.beProductOf(BH, u);
    F.at(1)+=1.0;
    F.at(2)+=1.0;
    F.at(3)+=1.0;

    Fmatrix.beMatrixForm(F);
    Finv.beInverseOf(Fmatrix);
    answer.beTranspositionOf(Finv);

}

void WeakPeriodicBoundaryCondition :: computeElementTangent(FloatMatrix &B, Element *e, int boundary, TimeStep *tStep)
{

    OOFEM_ERROR("Function obsolete");

    FloatArray gcoords;
    IntArray bnodes;

    FEInterpolation *geoInterpolation = e->giveInterpolation();

    // Use correct interpolation for the dofid on which the condition is applied
    FEInterpolation *interpolation = e->giveInterpolation( ( DofIDItem ) dofids(0) );

    interpolation->boundaryGiveNodes(bnodes, boundary);

    B.resize(bnodes.giveSize(), ndof);
    B.zero();

    std :: unique_ptr< IntegrationRule >iRule(geoInterpolation->giveBoundaryIntegrationRule(orderOfPolygon, boundary));

    for ( GaussPoint *gp: *iRule ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        FloatArray N;

        // Find the value of parameter s which is the vert/horiz distance to 0
        geoInterpolation->boundaryLocal2Global( gcoords, boundary, lcoords, FEIElementGeometryWrapper(e) );

        FloatArray boundaryLocal;
        interpolation->global2local( boundaryLocal, gcoords, FEIElementGeometryWrapper(e) );

#if 0
        FloatMatrix FinvT;
        if (tStep->giveNumber()>0) {
            if ( dynamic_cast <NLStructuralElement *> (e) != NULL ) { // If finite strains, compute deformation gradient
                this->computeDeformationGradient(FinvT, e, &boundaryLocal, tStep);
                FinvT.printYourself();
            }
        }

#endif

        // Compute base function values
        interpolation->boundaryEvalN( N, boundary, lcoords, FEIElementGeometryWrapper(e) );
        // Compute Jacobian
        double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(e) ) );

        for ( int j = 0; j < B.giveNumberOfColumns(); j++ ) {

            double fVal = computeBaseFunctionValue(j, gcoords);

            for ( int k = 0; k < B.giveNumberOfRows(); k++ ) {
                B(k, j) += N(k) * fVal * detJ * gp->giveWeight();
            }
        }
    }
}

void WeakPeriodicBoundaryCondition :: assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( type != TangentStiffnessMatrix && type != StiffnessMatrix ) {
        return;
    }

    IntArray c_loc, r_loc;
    gammaDman->giveLocationArray(gamma_ids, r_loc, r_s);
    gammaDman->giveLocationArray(gamma_ids, c_loc, c_s);

    FloatMatrix B, BT;
    int normalSign;

    updateSminmax();

    // Assemble each side
    for ( int thisSide = 0; thisSide <= 1; thisSide++ ) {
        normalSign = sideSign [ thisSide ];

        for ( size_t ielement = 0; ielement < element [ thisSide ].size(); ielement++ ) {               // Loop over each element on this edge
            Element *thisElement = this->domain->giveElement( element [ thisSide ].at(ielement) );
            int boundary = side [ thisSide ].at(ielement);

            // Find dofs for this element side
            IntArray r_sideLoc, c_sideLoc;

            // Find dofs for this element which should be periodic
            IntArray bNodes;

            FEInterpolation *interpolation = thisElement->giveInterpolation( ( DofIDItem ) dofids(0) );
            FEInterpolation *geoInterpolation = thisElement->giveInterpolation();

            interpolation->boundaryGiveNodes( bNodes, side [ thisSide ].at(ielement) );

            thisElement->giveBoundaryLocationArray(r_sideLoc, bNodes, dofids, r_s);
            thisElement->giveBoundaryLocationArray(c_sideLoc, bNodes, dofids, c_s);

            B.resize(bNodes.giveSize()*ndofids, ndofids*tcount);
            B.zero();

            std :: unique_ptr< IntegrationRule >iRule(geoInterpolation->giveBoundaryIntegrationRule(orderOfPolygon, boundary));

            for ( GaussPoint *gp: *iRule ) {
                FloatArray lcoords = gp->giveNaturalCoordinates();
                FloatArray N, gcoords;

                geoInterpolation->boundaryLocal2Global( gcoords, boundary, lcoords, FEIElementGeometryWrapper(thisElement));

                interpolation->boundaryEvalN(N, boundary, lcoords, FEIElementGeometryWrapper(thisElement));

                double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(thisElement) ) );

                FloatMatrix Mbeta(ndofids, ndof), Mv(ndofids, bNodes.giveSize()*ndofids), NvTNbeta;

                for (int i=0; i<tcount; i++) {
                    for (int j=0; j<ndofids; j++) {
                        Mbeta.at(j+1, ndofids*i+j+1) = computeBaseFunctionValue(i, gcoords);
                    }
                }

                for (int i=0; i<N.giveSize(); i++) {
                    for (int j=0; j<ndofids; j++) {
                        Mv.at(j+1, ndofids*i+j+1) = N.at(i+1);
                    }
                }

                FloatMatrix defNv, F, Finv;
                double J=1.0;

                if (nlgeo) {
                    FloatArray elocal;
                    geoInterpolation->global2local(elocal, gcoords, FEIElementGeometryWrapper(thisElement));
                    computeDeformationGradient(F, thisElement, &elocal, tStep);
//                    J=F.giveDeterminant();
                    Finv.beInverseOf(F);
                    defNv.beProductOf(Finv, Mv);
                } else {
                    defNv = Mv;
                }

                NvTNbeta.beTProductOf(defNv, Mbeta);

                NvTNbeta.times(J * detJ * gp->giveWeight());
                B.add(NvTNbeta);


/*                for (int i=0; i<ndof; i++) {
                    double fVal = computeBaseFunctionValue(i, gcoords);
                    for ( int k = 0; k < B.giveNumberOfRows(); k++ ) {
                        B(k, i) += N(k) * fVal * detJ * gp->giveWeight();
                    }
                }*/


            }

            B.times(normalSign);
            BT.beTranspositionOf(B);

            answer.assemble(r_sideLoc, c_loc, B);
            answer.assemble(r_loc, c_sideLoc, BT);
        }
    }

}

double
WeakPeriodicBoundaryCondition :: computeBaseFunctionValue(int baseID, FloatArray coordinate)
{
    if (this->domain->giveNumberOfSpatialDimensions() == 2) {
        return computeBaseFunctionValue1D(baseID, coordinate.at(surfaceIndexes.at(1)));
    } else {
        return computeBaseFunctionValue2D(baseID, {coordinate.at(surfaceIndexes.at(1)), coordinate.at(surfaceIndexes.at(2))});
    }
}

double WeakPeriodicBoundaryCondition :: computeBaseFunctionValue2D(int baseID, FloatArray coordinate)
{
    double fVal = 0.0;

    if ( useBasisType == monomial ) {
        int a, b;
        getExponents(baseID+1, a, b);
        fVal = pow(coordinate.at(1), a) * pow(coordinate.at(2), b);
    } else if ( useBasisType == legendre ) {
        for ( int i = 1; i <= ndof; i++ ) {
            int a, b;
            getExponents(i, a, b);
            fVal = fVal + gsMatrix.at(baseID+1, i) * pow(coordinate.at(1), a) * pow(coordinate.at(2), b);
        }
    }
    // printf("Value for u_%u att coordinate %f, %f is %f\n", baseID, coordinate.at(1), coordinate.at(2), fVal);
    return fVal;
}

double WeakPeriodicBoundaryCondition :: computeBaseFunctionValue1D(int baseID, double coordinate)
{
    double fVal = 0.0;
    FloatArray sideLength;

    // compute side lengths
    sideLength.resize( smax.giveSize() );
    for ( int i = 1; i <= smax.giveSize(); i++ ) {
        sideLength.at(i) = smax.at(i) - smin.at(i);
    }

    if ( useBasisType == monomial ) {
        fVal = pow(coordinate, baseID);
    } else if ( useBasisType == trigonometric ) {
        if ( baseID % 2 == 0 ) {           // Even (does not yet work in 3D)
            fVal = cos( ( ( double ) baseID ) / 2. * ( coordinate * 2. * M_PI / sideLength.at(1) ) );
        } else {
            fVal = sin( ( ( double ) baseID + 1 ) / 2. * ( coordinate * 2. * M_PI / sideLength.at(1) ) );
        }
    } else if ( useBasisType == legendre ) {
        double n = ( double ) baseID;
        coordinate = 2.0 * coordinate - 1.0;
        for ( int k = 0; k <= baseID; k++ ) {
            fVal = fVal + binomial(n, k) * binomial(-n - 1.0, k) * pow( ( 1.0 - coordinate ) / 2.0, ( double ) k );
        }
    }

    return fVal;
}

void WeakPeriodicBoundaryCondition :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                     CharType type, ValueModeType mode,
                                                     const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    if ( type == InternalForcesVector ) {
        giveInternalForcesVector(answer, tStep, type, mode, s);
    } else if ( type == ExternalForcesVector ) {
        giveExternalForcesVector(answer, tStep, type, mode, s);
    }

}

void
WeakPeriodicBoundaryCondition :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep,
                                                          CharType type, ValueModeType mode,
                                                          const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    // Fetch unknowns of this boundary condition
    IntArray gammaLoc;
    FloatArray gamma;
    gammaDman->giveUnknownVector(gamma, gamma_ids, mode, tStep);
    gammaDman->giveLocationArray(gamma_ids, gammaLoc, s);

    // Values from solution
    FloatArray a;
    // Find dofs for this element side
    IntArray sideLocation, masterDofIDs;

    FloatMatrix B;

    int normalSign;

    updateSminmax();

    // Assemble each side
    for ( int thisSide = 0; thisSide <= 1; thisSide++ ) {
        normalSign = sideSign [ thisSide ];

        for ( size_t ielement = 0; ielement < element [ thisSide ].size(); ielement++ ) {                   // Loop over each element on this edge
            Element *thisElement = this->domain->giveElement( element [ thisSide ].at(ielement) );
            int boundary = side [ thisSide ].at(ielement);

            // Find dofs for this element which should be periodic
            IntArray bNodes;

            FEInterpolation *interpolation = thisElement->giveInterpolation( ( DofIDItem ) dofids(0) );
            FEInterpolation *geoInterpolation = thisElement->giveInterpolation();

            interpolation->boundaryGiveNodes( bNodes, boundary );

            thisElement->giveBoundaryLocationArray(sideLocation, bNodes, dofids, s, &masterDofIDs);
            thisElement->computeBoundaryVectorOf(bNodes, dofids, VM_Total, tStep, a);

            FloatArray vProd, gammaProd;

            B.resize(bNodes.giveSize(), ndof);
            B.zero();

            std :: unique_ptr< IntegrationRule >iRule(geoInterpolation->giveBoundaryIntegrationRule(orderOfPolygon, boundary));

            // Where we test with velocity
            vProd.resize(bNodes.giveSize()*dofids.giveSize());
            vProd.zero();

            // Where we test with gamma
            gammaProd.resize(ndof);
            gammaProd.zero();

            // For test:
            FloatArray myProd(gammaProd.giveSize()), myProdGamma(vProd.giveSize());

            for ( GaussPoint *gp: *iRule ) {
                FloatArray lcoords = gp->giveNaturalCoordinates();
                FloatArray N, gcoords;
                FloatMatrix C, D, Nbeta, Nv;

                geoInterpolation->boundaryLocal2Global( gcoords, boundary , lcoords, FEIElementGeometryWrapper(thisElement));
                interpolation->boundaryEvalN(N, boundary, lcoords, FEIElementGeometryWrapper(thisElement));
                double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(thisElement) ) );

                Nbeta.resize(ndofids, ndof);
                Nbeta.zero();
                for (int i=0; i<tcount; i++) {
                    double val = computeBaseFunctionValue(i, gcoords);
                    for (int j=0; j<ndofids; j++) {
                        Nbeta.at(j+1, i*ndofids+j+1) = val;
                    }
                }

                Nv.resize(ndofids, N.giveSize()*ndofids);
                Nv.zero();
                for (int i=0; i<ndofids; i++) {
                    for (int j=0; j<N.giveSize(); j++) {
                        Nv.at(i+1, ndofids*j+i+1) = N(j);
                    }
                }

                FloatMatrix defNv, F, Finv;
                double J=1.0;

                if (nlgeo) {
                    FloatArray elocal;
                    geoInterpolation->global2local(elocal, gcoords, FEIElementGeometryWrapper(thisElement));
                    computeDeformationGradient(F, thisElement, &elocal, tStep);
//                    J = F.giveDeterminant();
                    Finv.beInverseOf(F);
                    defNv.beProductOf(Finv, Nv);
                } else {
                    J=1.0;
                    defNv = Nv;
                }

                C.beTProductOf(Nbeta, defNv);
                D.beTranspositionOf(C);

                gammaProd.plusProduct(D, a, J*detJ*gp->giveWeight()*normalSign);
                vProd.plusProduct(C, gamma, J*detJ*gp->giveWeight()*normalSign);

                // Old:
/*                FloatArray BaseFunctionValues;
                BaseFunctionValues.resize(ndof);
                BaseFunctionValues.zero();
                for (int i=0; i<ndof; i++) {
                    BaseFunctionValues.at(i+1) = computeBaseFunctionValue(i, gcoords);
                }

                C.beDyadicProductOf(N, BaseFunctionValues);
                D.beTranspositionOf(C);

                myProd.plusProduct(C, a, detJ*gp->giveWeight()*normalSign);
                myProdGamma.plusProduct(D, gamma, detJ*gp->giveWeight()*normalSign); */

            }

            if ( eNorms ) {
                eNorms->assembleSquared( vProd, gamma_ids );
                eNorms->assembleSquared( gammaProd, masterDofIDs);
            }

            answer.assemble(gammaProd, gammaLoc);
            answer.assemble(vProd, sideLocation);
//            answer.assemble(myProd, gammaLoc);
//            answer.assemble(myProdGamma, sideLocation);
        }
    }
}

void
WeakPeriodicBoundaryCondition :: giveExternalForcesVector(FloatArray &answer, TimeStep *tStep,
                                                          CharType type, ValueModeType mode,
                                                          const UnknownNumberingScheme &s)
{

    updateSminmax();

    IntArray gammaLoc;
    gammaDman->giveLocationArray(gamma_ids, gammaLoc, s);

    FloatArray temp;
    temp.resize(ndof);
    temp.zero();

    int normalSign;

    for ( int thisSide = 0; thisSide <= 1; thisSide++ ) {
        normalSign = sideSign [ thisSide ];
        for ( size_t ielement = 0; ielement < element [ thisSide ].size(); ielement++ ) {                   // Loop over each element on this edge

            Element *thisElement = this->domain->giveElement( element [ thisSide ].at(ielement) );

            FEInterpolation *geoInterpolation = thisElement->giveInterpolation();

            std :: unique_ptr< IntegrationRule >iRule(geoInterpolation->giveBoundaryIntegrationRule(orderOfPolygon, side [ thisSide ].at(ielement) ));

            for ( GaussPoint *gp: *iRule ) {

                FloatArray gcoords;
                FloatArray lcoords = gp->giveNaturalCoordinates();

                // Find the value of parameter s which is the vert/horiz distance to 0
                geoInterpolation->boundaryLocal2Global( gcoords, side [ thisSide ].at(ielement), lcoords, FEIElementGeometryWrapper( thisElement ) );
                // Compute Jacobian
                double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( side [ thisSide ].at(ielement), lcoords, FEIElementGeometryWrapper(thisElement) ) );

                for (int j=0; j<ndof; j++) {
                    FloatArray coord;

                    double fVal = computeBaseFunctionValue(j, gcoords );

                    temp.at(j+1)=temp.at(j+1) + normalSign*this->g.dotProduct(gcoords)*fVal*gp->giveWeight()*detJ;
                }
            }
        }
    }

    answer.assemble(temp, gammaLoc);

    // Finally, compute value of loadtimefunction
    double factor;
    factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer.times(factor);
}

int WeakPeriodicBoundaryCondition :: giveNumberOfInternalDofManagers()
{
    return 1;
}

DofManager *WeakPeriodicBoundaryCondition :: giveInternalDofManager(int i)
{
    if ( i == 1 ) {
        return gammaDman.get();
    } else {
        return NULL;
    }
}

double WeakPeriodicBoundaryCondition :: factorial(int n)
{
    int x = 1;
    for ( int i = 1; i <= n; i++ ) {
        x = x * i;
    }
    return x;
}

double WeakPeriodicBoundaryCondition :: binomial(double n, int k)
{
    double f = 1.0;
    for ( int i = 1; i <= k; i++ ) {
        f = f * ( n - ( k - i ) ) / i;
    }
    return f;
}

void WeakPeriodicBoundaryCondition :: getExponents(int term, int &i, int &j)
{
    bool doContinue = true;

    // c is the number of the current term
    int c = 0;

    // n is the order of the current polynomial (row in Pascals triangle)
    int n = 0;

    while ( doContinue ) {
        for ( int t = 0; t <= n; t++ ) {
            c++;
            if ( c == term ) {
                i = n - t;
                j = t;
                doContinue = false;
                break;
            }
        }
        n++;
    }
}
}
