/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
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

#include "stokesflow.h"
#include "stokesflowvelocityhomogenization.h"
#include "primaryfield.h"
#include "classfactory.h"
#include "deadweight.h"
#include "tr21stokes.h"
#include "dofmanager.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_EngngModel(StokesFlowVelocityHomogenization);

StokesFlowVelocityHomogenization :: StokesFlowVelocityHomogenization(int i, EngngModel *_master) : StokesFlow(i, _master)
{
}

StokesFlowVelocityHomogenization :: ~StokesFlowVelocityHomogenization()
{ }


double
StokesFlowVelocityHomogenization :: giveAreaOfRVE()
{
    FloatArray min, max;

    min = * this->giveDomain(1)->giveDofManager(1)->giveCoordinates();
    max = * this->giveDomain(1)->giveDofManager(1)->giveCoordinates();

    for ( auto &node : this->giveDomain(1)->giveDofManagers() ) {
        min.beMinOf( min, * node->giveCoordinates() );
        max.beMaxOf( max, * node->giveCoordinates() );
    }

    max.subtract(min);
    return max.product();
}


void
StokesFlowVelocityHomogenization :: getMeans(FloatArray &v, TimeStep *tStep)
{
    FloatMatrix N;
    FloatArray v_hatTemp, unknowns;

    v.clear();

    for ( auto &elem : this->giveDomain(1)->giveElements() ) {
        this->integrateNMatrix(N, *elem, tStep);
        elem->computeVectorOf({V_u, V_v, V_w}, VM_Total, tStep, unknowns);
        v_hatTemp.beProductOf(N, unknowns);
        v.add(v_hatTemp);
    }

    v.times( 1. / this->giveAreaOfRVE() );
}


void
StokesFlowVelocityHomogenization :: computeTangent(FloatMatrix &answer, TimeStep *tStep)
{
    IntArray loc, col;

    Domain *domain = this->giveDomain(1);
    int nsd = domain->giveNumberOfSpatialDimensions();
    int ndof = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    // Build F matrix
    IntArray dofs = {V_u, V_v}; ///@todo Should not be hardcoded to just V_u and V_v
    FloatMatrix F(ndof, nsd), Fe, N;
    col.enumerate(nsd);

    for ( auto &elem : domain->giveElements() ) {

        this->integrateNMatrix(N, *elem, tStep);
        
        elem->giveLocationArray( loc, dofs, EModelDefaultEquationNumbering() );
        Fe.beTranspositionOf(N);
        F.assemble(Fe, loc, col);
    }

    FloatMatrix H;

    std :: unique_ptr< SparseLinearSystemNM > linMethod( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) );

    H.resize( F.giveNumberOfRows(), F.giveNumberOfColumns() );
    H.zero();

    linMethod->solve(*stiffnessMatrix, F, H);

    answer.beTProductOf(H, F);
    answer.times( 1. / this->giveAreaOfRVE() );
}


void
StokesFlowVelocityHomogenization :: integrateNMatrix(FloatMatrix &N, Element &elem, TimeStep *tStep)
{
    FloatArray n, n2;

    for ( GaussPoint *gp: *elem.giveDefaultIntegrationRulePtr() ) {
        FloatArray &lcoords = * gp->giveNaturalCoordinates();

        ///@todo Ask the element for the N-matrix instead
        elem.giveInterpolation()->evalN( n, lcoords, FEIElementGeometryWrapper(&elem) );
        double detJ = elem.giveInterpolation()->giveTransformationJacobian( lcoords, FEIElementGeometryWrapper(&elem) );
        n2.add(gp->giveWeight() * detJ, n);
    }

    N.beNMatrixOf(n2, this->giveDomain(1)->giveNumberOfSpatialDimensions());
}

}
