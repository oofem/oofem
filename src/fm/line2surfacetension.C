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

#include "line2surfacetension.h"
#include "node.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "material.h"
#include "fei2dlinequad.h"

namespace oofem {
FEI2dLineQuad Line2SurfaceTension :: fei(1, 2);

Line2SurfaceTension :: Line2SurfaceTension(int n, Domain *aDomain) : LineSurfaceTension(n, aDomain)
{
    numberOfDofMans  = 3;
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ 1 ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 2, _Unknown);
}

Line2SurfaceTension :: ~Line2SurfaceTension()
{
}

void Line2SurfaceTension :: computeN(FloatArray &answer, const FloatArray &lcoords) const
{
    this->fei.evalN(answer, lcoords, FEIElementGeometryWrapper(this));
}

FEInterpolation *Line2SurfaceTension :: giveInterpolation()
{
    return &this->fei;
}

double Line2SurfaceTension :: computeNXIntegral() const
{
    Node *node;
    double x1, x2, x3, y1, y2, y3;

    node = this->giveNode(1);
    x1 = node->giveCoordinate(1);
    y1 = node->giveCoordinate(2);

    node = this->giveNode(2);
    x2 = node->giveCoordinate(1);
    y2 = node->giveCoordinate(2);

    node = this->giveNode(3);
    x3 = node->giveCoordinate(1);
    y3 = node->giveCoordinate(2);

    return (x3*(8*y1 - y2) + 2*x1*(y2 - 4*y3) + x2*(-2*y1 + y3))/6;
}

void Line2SurfaceTension :: computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep)
{
    ///@todo Support axisymm.
    //domainType dt = this->giveDomain()->giveDomainType();
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    double t = 1, gamma_s;
    ///@todo Should i use this? Not used in FM module (but perhaps it should?) / Mikael.
    //t = this->giveDomain()->giveCrossSection(1)->give(CS_Thickness);
    gamma_s = this->giveMaterial()->give('g', NULL);

    FloatMatrix xy(2, 3);
    Node *node;
    for ( int i = 1; i <= 3; i++ ) {
        node = giveNode(i);
        xy.at(1, i) = node->giveCoordinate(1);
        xy.at(2, i) = node->giveCoordinate(2);
    }

    FloatArray A;
    FloatArray dNdxi(3);
    FloatArray es(2); // tangent vector to curve
    FloatMatrix BJ(2, 6);
    BJ.zero();

    answer.resize(6);
    answer.zero();

    for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        //interpolation.evaldNdx(dN, domain, dofManArray, * gp->giveCoordinates(), 0.0);
        double xi = gp->giveCoordinate(1);

        // Some simplifications can be performed, since the mapping J is a scalar.
        dNdxi.at(1) = -0.5 + xi;
        dNdxi.at(2) =  0.5 + xi;
        dNdxi.at(3) = -2.0 * xi;

        es.beProductOf(xy, dNdxi);
        double J = es.computeNorm();
        es.times(1 / J); //es.normalize();

        // dNds = dNdxi/J
        // B.at(1,1) = dNds.at(1); and so on.

        BJ.at(1, 1) = BJ.at(2, 2) = dNdxi.at(1);
        BJ.at(1, 3) = BJ.at(2, 4) = dNdxi.at(2);
        BJ.at(1, 5) = BJ.at(2, 6) = dNdxi.at(3);

        A.beTProductOf(BJ, es);
        answer.add( - gamma_s * t * gp->giveWeight(), A); // Note! Negative sign!
    }
}

void Line2SurfaceTension :: computeTangent(FloatMatrix &answer, TimeStep *tStep)
{
#if 1
    answer.resize(6, 6);
    answer.zero();
#else
    ///@todo Support axisymm.
    domainType dt = this->giveDomain()->giveDomainType();
    if (dt == _3dAxisymmMode) {
        OOFEM_ERROR("Line2SurfaceTension :: computeTangent - Axisymm not implemented");
    }
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    double t = 1, gamma_s;
    ///@todo Should i use this? Not that meaningful for flow problems.
    //t = this->giveDomain()->giveCrossSection(1)->give(CS_Thickness);
    gamma_s = this->giveMaterial()->give('g', NULL);

    FloatMatrix xy(2,3);
    Node *node;
    for (int i = 1; i <= 3; i++) {
        node = giveNode(i);
        xy.at(1,i) = node->giveCoordinate(1);
        xy.at(2,i) = node->giveCoordinate(2);
    }

    FloatArray A;
    FloatArray dNdxi(3);
    FloatArray es(2); // tangent vector to curve
    FloatMatrix BJ(2,6);
    BJ.zero();
    FloatMatrix temp1,temp2;

    answer.resize(6,6);
    answer.zero();
    for (int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);

        double xi = gp->giveCoordinate(1);

        dNdxi.at(1) = -0.5+xi;
        dNdxi.at(2) =  0.5+xi;
        dNdxi.at(3) = -2.0*xi;

        es.beProductOf(xy,dNdxi);
        double J = es.computeNorm();
        es.times(1/J); //es.normalize();

        BJ.at(1,1) = BJ.at(2,2) = dNdxi.at(1);
        BJ.at(1,3) = BJ.at(2,4) = dNdxi.at(2);
        BJ.at(1,5) = BJ.at(2,6) = dNdxi.at(3);

        A.beTProductOf(BJ,es);

        temp1.beTProductOf(BJ,BJ);
        temp2.beDyadicProductOf(A,A);
        temp1.subtract(temp2);
        temp1.times(t*gp->giveWeight()/J*(tStep->giveTimeIncrement()));
        answer.add(temp1);
    }
    answer.times(gamma_s);
#endif
}

double Line2SurfaceTension :: SpatialLocalizerI_giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords)
{
    if (!this->computeLocalCoordinates(lcoords, gcoords)) {
        lcoords.resize(0);
        closest.resize(0);
        return -1.0;
    }
    // compute local coordinates already gives closest point.
    this->computeGlobalCoordinates(closest, lcoords);
    return closest.distance(gcoords);
}

double Line2SurfaceTension :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &gcoords)
{
    return this->giveNode(3)->giveCoordinates()->distance(gcoords);
}


int Line2SurfaceTension :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
        TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    FloatArray lcoords, closest;
    double dist = this->SpatialLocalizerI_giveClosestPoint(lcoords, closest, gcoords);
    if (dist < 0) {
        answer.resize(0);
        return false;
    }
    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void Line2SurfaceTension :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n;
    this->computeN(n, lcoords);

    answer.resize(2);
    answer.zero();
    for (int i = 1; i <= n.giveSize(); i++) {
        answer(0) += n.at(i)*this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(EID_MomentumBalance, mode, tStep);
        answer(1) += n.at(i)*this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(EID_MomentumBalance, mode, tStep);
    }
}

} // end namespace oofem
