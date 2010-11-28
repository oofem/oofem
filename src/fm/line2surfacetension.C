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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

void Line2SurfaceTension :: computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep)
{
    //domainType dt = this->giveDomain()->giveDomainType(); // TODO, support axisymm
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    double t;
    //t = this->giveDomain()->giveCrossSection(1)->give(CS_Thickness); // TODO: Should i use this?
    t = 1;

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
        dNdxi.at(3) =   -2 * xi;

        es.beProductOf(xy, dNdxi);
        double J = es.computeNorm();
        es.times(1 / J); //es.normalize();

        // dNds = dNdxi/J
        // B.at(1,1) = dNds.at(1); and so on.

        BJ.at(1, 1) = BJ.at(2, 2) = dNdxi.at(1);
        BJ.at(1, 3) = BJ.at(2, 4) = dNdxi.at(2);
        BJ.at(1, 5) = BJ.at(2, 6) = dNdxi.at(3);

        A.beTProductOf(BJ, es);
        A.times( -this->gamma_s * t * gp->giveWeight() ); // Note! Negative sign!

        answer.add(A);
    }

    // This part is not verified yet.
    if ( this->bflag1 ) {
        dNdxi.at(1) = -1.5;
        dNdxi.at(2) = -0.5;
        dNdxi.at(3) = 2;
        es.beProductOf(xy, dNdxi);
        es.normalize();
        answer.at(1) -= es.at(1) * this->gamma_s;
        answer.at(2) -= es.at(2) * this->gamma_s;
    }

    if ( this->bflag2 ) {
        dNdxi.at(1) = 0.5;
        dNdxi.at(2) = 1.5;
        dNdxi.at(3) = -2;
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
} // end namespace oofem
