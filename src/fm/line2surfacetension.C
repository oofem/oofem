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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "timestep.h"
#include "mathfem.h"

#include <stdlib.h>
#include <math.h>

namespace oofem {
Line2SurfaceTension :: Line2SurfaceTension(int n, Domain *aDomain) : LineSurfaceTension(n, aDomain)
{
    numberOfDofMans  = 3;
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ 1 ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 2, _Unknown);
}

Line2SurfaceTension :: ~Line2SurfaceTension() {}

void Line2SurfaceTension :: computeN(FloatArray &answer, const FloatArray &lcoords) const
{
    double xi = lcoords(0);
    answer.resize(3);
    answer(0) = 0.5*(xi-1.0)*xi;
    answer(1) = 0.5*(xi+1.0)*xi;
    answer(2) = 1.0-xi*xi;
}

int Line2SurfaceTension :: computeLocalCoordinates(FloatArray &lcoords, const FloatArray &gcoords)
{
    FloatArray x1_x2, p_x3, x3_x2_x1;
    x1_x2 = *this->giveNode(1)->giveCoordinates();
    x1_x2.subtract(*this->giveNode(2)->giveCoordinates());
    p_x3 = gcoords;
    p_x3.subtract(*this->giveNode(3)->giveCoordinates());
    x3_x2_x1 = *this->giveNode(3)->giveCoordinates();
    x3_x2_x1.add(*this->giveNode(3)->giveCoordinates());
    x3_x2_x1.subtract(*this->giveNode(2)->giveCoordinates());
    x3_x2_x1.subtract(*this->giveNode(1)->giveCoordinates());

    double b0 = 0.50*x1_x2.dotProduct(p_x3); // 1/2*(x1 - x2).(p - x3)
    double b1 = 0.25*x1_x2.computeSquaredNorm() + x3_x2_x1.dotProduct(p_x3); // 1/4*(x1 - x2)^2 + (2*x3 - x2 - x1).(p - x3)
    double b2 = 0.75*x1_x2.dotProduct(x3_x2_x1); // 3/4*(x1-x2).(2x3 - x2 - x1)
    double b3 = 0.50*x3_x2_x1.computeSquaredNorm(); // 1/2*(2*x3 - x2 - x3)^2

    double r[3];
    int roots;
    cubic(b3, b2, b1, b0, &r[0], &r[1], &r[2], &roots);

    // Copy them over, along with boundary cases.
    int points = 2;
    double p[5] = {-1.0, 1.0};
    for (int i = 0; i < roots; i++) {
        if (r[i] > -1.0 && r[i] < 1.0) {
            // The cubic solver has pretty bad accuracy, performing a single newton iteration
            // which typically improves the solution by many many orders of magnitude.
            // TODO: Move this to the cubic solver.
            r[i] -= (b0 + b1*r[i] + b2*r[i]*r[i] + b3*r[i]*r[i]*r[i])/(b1 + 2*b2*r[i] + 3*b3*r[i]*r[i]);
            //printf("r[%d] -> %e\n",i, (b0 + b1*r[i] + b2*r[i]*r[i] + b3*r[i]*r[i]*r[i]));
            p[points] = r[i];
            points++;
        }
    }
/*
    int points = 101;
    double p[points];
    for (int i = 0; i < points; i++) {
        p[i] = -1.0 + 2.0*i/(points-1);
    }
*/
    double min_distance2 = 0.0, min_xi, distance2;
    FloatArray f(2), xi(1);

    for (int i = 0; i < points; i++) {
        xi(0) = p[i];
        this->computeGlobalCoordinates(f,xi);
        distance2 = f.distance_square(gcoords);
        if ( i == 0 || distance2 < min_distance2 ) {
            min_distance2 = distance2;
            min_xi = xi(0);
        }
    }
    // Local coordinates would be preferable..
    lcoords.resize(1);
    lcoords(0) = min_xi;

    return true;
}

void Line2SurfaceTension :: computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep)
{
    //domainType dt = this->giveDomain()->giveDomainType(); // TODO, support axisymm
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    double t = 1;
    //t = this->giveDomain()->giveCrossSection(1)->give(CS_Thickness); // TODO: Should i use this?

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
        answer.add( -this->gamma_s * t * gp->giveWeight(), A); // Note! Negative sign!
    }

    // This part is not verified yet.
    if ( this->bflag1 ) {
        dNdxi.at(1) = -1.5;
        dNdxi.at(2) = -0.5;
        dNdxi.at(3) =  2.0;
        es.beProductOf(xy, dNdxi);
        es.normalize();
        answer.at(1) -= es.at(1) * this->gamma_s;
        answer.at(2) -= es.at(2) * this->gamma_s;
    }

    if ( this->bflag2 ) {
        dNdxi.at(1) =  0.5;
        dNdxi.at(2) =  1.5;
        dNdxi.at(3) = -2.0;
        es.beProductOf(xy, dNdxi);
        es.normalize();
        answer.at(3) -= es.at(1) * this->gamma_s;
        answer.at(4) -= es.at(2) * this->gamma_s;
    }
}

void Line2SurfaceTension :: computeTangent(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();
    /*
     *  domainType dt = this->giveDomain()->giveDomainType(); // TODO, support axisymm
     *  IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
     *  double thickness = 1;
     *
     *  FloatMatrix xy(2,3);
     *  Node *node;
     *  for (int i = 1; i <= 3; i++) {
     *      node = giveNode(i);
     *      xy.at(1,i) = node->giveCoordinate(1);
     *      xy.at(2,i) = node->giveCoordinate(2);
     *  }
     *
     *  FloatArray A;
     *  FloatArray dNdxi(3);
     *  FloatArray es(2); // tangent vector to curve
     *  FloatMatrix BJ(2,6);
     *  BJ.zero();
     *  FloatMatrix temp1,temp2;
     *
     *  answer.resize(6,6);
     *  answer.zero();
     *  for (int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
     *      GaussPoint *gp = iRule->getIntegrationPoint(k);
     *
     *      double xi = gp->giveCoordinate(1);
     *
     *      dNdxi.at(1) = -0.5+xi;
     *      dNdxi.at(2) =  0.5+xi;
     *      dNdxi.at(3) =   -2*xi;
     *
     *      es.beProductOf(xy,dNdxi);
     *      double J = es.computeNorm();
     *      es.times(1/J); //es.normalize();
     *
     *      BJ.at(1,1) = BJ.at(2,2) = dNdxi.at(1);
     *      BJ.at(1,3) = BJ.at(2,4) = dNdxi.at(2);
     *      BJ.at(1,5) = BJ.at(2,6) = dNdxi.at(3);
     *
     *      A.beTProductOf(BJ,es);
     *
     *      temp1.beTProductOf(BJ,BJ);
     *      temp2.beDyadicProductOf(A,A);
     *      temp1.subtract(temp2);
     *      temp1.times(this->gamma_s*thickness*gp->giveWeight()/J*(tStep->giveTimeIncrement()));
     *      answer.add(temp1);
     *  }
     *  if (dt ==
     */
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

int Line2SurfaceTension :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
        TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    FloatArray lcoords, closest;
    int ok = this->SpatialLocalizerI_giveClosestPoint(lcoords, closest, gcoords);
    if (!ok) {
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
