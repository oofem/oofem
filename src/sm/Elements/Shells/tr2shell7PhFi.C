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

#include "tr2shell7phfi.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "vtkxmlexportmodule.h"
#include "classfactory.h"
#include "tr2shell7.h"

namespace oofem {
REGISTER_Element(Tr2Shell7PhFi);

FEI3dTrQuad Tr2Shell7PhFi :: interpolation;

IntArray Tr2Shell7PhFi :: ordering_disp(42);
IntArray Tr2Shell7PhFi :: ordering_disp_inv(42);
IntArray Tr2Shell7PhFi :: ordering_base(42);
IntArray Tr2Shell7PhFi :: ordering_base_inv(42); 
IntArray Tr2Shell7PhFi:: ordering_all(1);
IntArray Tr2Shell7PhFi:: ordering_damage(1);
IntArray Tr2Shell7PhFi :: ordering_gr(1);
IntArray Tr2Shell7PhFi :: ordering_gr_edge(21);
bool Tr2Shell7PhFi :: __initialized = Tr2Shell7PhFi :: initOrdering();


Tr2Shell7PhFi:: Tr2Shell7PhFi(int n, Domain *aDomain) : Shell7BasePhFi(n, aDomain)
{
    this->numberOfDofMans = 6;

    this->ordering_damage.resize(0);
    this->ordering_gr.resize(0);
    this->ordering_disp_inv.resize(0);

    IntArray localDam, localDisp(7);	// hard coded for 7 parameter shell!!

    localDam.resize(0);
    localDisp.setValues(7, 1, 2, 3, 19, 20, 21, 37);

    for (int i = 1; i <= numberOfLayers; i++) {
        //localDam.followedBy(i + this->giveNumberOfuDofs());
        localDam.followedBy(7 + i);
    }

    int numDofsPerNode = 7 + numberOfLayers;

    for (int i = 1; i <= this->numberOfDofMans; i++) {
        this->ordering_damage.followedBy(localDam);
        this->ordering_gr.followedBy(localDisp);
        this->ordering_gr.followedBy(localDam);
        this->ordering_disp_inv.followedBy(localDisp);

        localDisp.at(1) += 3;
        localDisp.at(2) += 3;
        localDisp.at(3) += 3;
        localDisp.at(4) += 3;
        localDisp.at(5) += 3;
        localDisp.at(6) += 3;
        localDisp.at(7) += 1;

        localDam.add(numberOfLayers);
    }


    this->ordering_all.resize(0);
    this->ordering_all = this->ordering_disp;
    this->ordering_all.followedBy(ordering_damage);

    // JB - New

    IntArray temp_x(0), temp_m(0), temp_dam(0), local_temp_x(0), local_temp_m(0), local_temp_gam(0), local_temp_dam(0);
    temp_x.setValues(3, 1, 2, 3);
    temp_m.setValues(3, 4, 5, 6);
    int temp_gam = 7;

    for (int i = 1; i <= numberOfLayers; i++) {
        temp_dam.followedBy(7 + i);
    }

    for (int i = 1; i <= this->numberOfDofMans; i++) {
        local_temp_x.followedBy(temp_x);
        local_temp_m.followedBy(temp_m);
        local_temp_gam.followedBy(temp_gam);
        local_temp_dam.followedBy(temp_dam);
        temp_x.add(numDofsPerNode);
        temp_m.add(numDofsPerNode);
        temp_gam +=numDofsPerNode;
        temp_dam.add(numDofsPerNode);
    }
    //local_temp_x.printYourself();
    //local_temp_m.printYourself();
    //local_temp_gam.printYourself();
    //local_temp_dam.printYourself();

    // Construct the two orddering arrays
    ordering_disp.resize(0);
    ordering_damage.resize(0);
    this->ordering_disp.followedBy(local_temp_x); // xbar field
    this->ordering_disp.followedBy(local_temp_m); // director field
    this->ordering_disp.followedBy(local_temp_gam); // gamma field
    this->ordering_damage = local_temp_dam;

    //ordering_disp.printYourself();
    //ordering_damage.printYourself();
}

void
Tr2Shell7PhFi :: giveDofManDofIDMask_u(IntArray &answer)
{
    Shell7BasePhFi :: giveDofManDofIDMask_u(answer);
}

void
Tr2Shell7PhFi :: giveDofManDofIDMask_d(IntArray &answer)
{
    Shell7BasePhFi :: giveDofManDofIDMask_d(answer);
}

const IntArray &
Tr2Shell7PhFi:: giveOrdering(SolutionField fieldType) const
{
    OOFEM_ERROR("Tr2Shell7PhFi :: giveOrdering not implemented: Use Tr2Shell7PhFi :: giveOrderingPhFi instead");
    return 0;
}


const IntArray &
Tr2Shell7PhFi:: giveOrderingPhFi(SolutionFieldPhFi fieldType) const
{
    if ( fieldType == All ) {
        return this->ordering_all;
    } else if ( fieldType == Displacement ) {
        return this->ordering_disp;
    } else if ( fieldType == Damage ) {
        return this->ordering_damage;
    } else if ( fieldType == AllInv ) {
        return this->ordering_gr;
    } else { /*if ( fieldType == EdgeInv )*/
        OOFEM_ERROR("Tr2Shell7PhFi :: giveOrdering: the requested ordering is not implemented")
        //return this->ordering_gr_edge;
    }
}

const IntArray &
Tr2Shell7PhFi :: giveOrdering_All() const
{
    return this->ordering_base; // When shell7base methods are called this is what is wanted
}

const IntArray &
Tr2Shell7PhFi :: giveOrdering_AllInv() const
{
    // Same as in Shell7Base 
    return this->ordering_base_inv;
}

const IntArray &
Tr2Shell7PhFi :: giveOrdering_Disp() const
{
    return this->ordering_disp; 
}

const IntArray &
Tr2Shell7PhFi :: giveOrdering_Damage() const
{
    return this->ordering_damage;
}


void
Tr2Shell7PhFi:: giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords)
{
    nodeLocalXiCoords.setValues(6, 1., 0., 0., .5, 0., .5);      // corner nodes then midnodes, uncertain of node numbering
    nodeLocalEtaCoords.setValues(6, 0., 1., 0., .5, .5, 0.);
}


FEInterpolation *Tr2Shell7PhFi:: giveInterpolation() const { return & interpolation; }


void
Tr2Shell7PhFi:: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        int nPointsTri  = 6;   // points in the plane
        int nPointsEdge = 2;   // edge integration
        specialIntegrationRulesArray = new IntegrationRule * [ 3 ];

        // Midplane and thickness

        // Midplane (Mass matrix integrated analytically through the thickness)
        specialIntegrationRulesArray [ 1 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 1 ]->SetUpPointsOnWedge(nPointsTri, 1, _3dMat); //@todo replce with triangle which has a xi3-coord


        // Edge
        specialIntegrationRulesArray [ 2 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 2 ]->SetUpPointsOnLine(nPointsEdge, _3dMat);


        // Layered cross section for bulk integration
        //@todo - must use a cast here since check consistency has not been called yet
        LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( Tr2Shell7PhFi:: giveCrossSection() );
        if ( layeredCS == NULL ) {
            OOFEM_ERROR("Tr2Shell7PhFionly supports layered cross section");
        }
        this->numberOfIntegrationRules = layeredCS->giveNumberOfLayers();
        this->numberOfGaussPoints = layeredCS->giveNumberOfLayers() * nPointsTri * layeredCS->giveNumIntegrationPointsInLayer();
        layeredCS->setupLayeredIntegrationRule(integrationRulesArray, this, nPointsTri);


        // Thickness integration for stress recovery
        specialIntegrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 0 ]->SetUpPointsOnLine(layeredCS->giveNumIntegrationPointsInLayer(), _3dMat);
    }
}


void
Tr2Shell7PhFi:: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge == 1 ) {        // edge between nodes 1-4-2
        //answer.setValues(21, 1, 2, 3, 4, 5, 6, 7,   8, 9, 10, 11, 12, 13, 14,   22, 23, 24, 25, 26, 27, 28);
        answer.setValues(21, 1, 2, 3, 8, 9, 10, 22, 23, 24,  4, 5, 6, 11, 12, 13, 25, 26, 27,   7, 14, 28);
    } else if ( iEdge == 2 ) { // edge between nodes 2-5-3
        //answer.setValues(21,   8, 9, 10, 11, 12, 13, 14,   15, 16, 17, 18, 19, 20, 21,   29, 30, 31, 32, 33, 34, 35 );
        answer.setValues(21,   8, 9, 10, 15, 16, 17, 29, 30, 31,   11, 12, 13, 18, 19, 20, 32, 33, 34,   14, 21, 35);
    } else if ( iEdge == 3 ) { // edge between nodes 3-6-1
        //answer.setValues(21,   15, 16, 17, 18, 19, 20, 21,   1, 2, 3, 4, 5, 6, 7,   36, 37, 38, 39, 40, 41, 42);
        answer.setValues(21,   15, 16, 17, 1, 2, 3, 36, 37, 38,   18, 19, 20, 4, 5, 6, 39, 40, 41,   21, 7, 42);
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}


void
Tr2Shell7PhFi::giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(42);
    for ( int i = 1; i <= 42; i++ ) {
        answer.at(i) = i;
    }
}


double
Tr2Shell7PhFi::computeAreaAround(GaussPoint *gp, double xi)
{
    FloatArray G1, G2, temp;
    FloatMatrix Gcov;
    FloatArray lcoords(3);
    lcoords.at(1) = gp->giveCoordinate(1);
    lcoords.at(2) = gp->giveCoordinate(2);
    lcoords.at(3) = xi;
    this->evalInitialCovarBaseVectorsAt(lcoords, Gcov);
    G1.beColumnOf(Gcov, 1);
    G2.beColumnOf(Gcov, 2);
    temp.beVectorProductOf(G1, G2);
    double detJ = temp.computeNorm();
    return detJ * gp->giveWeight();
}


double
Tr2Shell7PhFi::computeVolumeAroundLayer(GaussPoint *gp, int layer)
{
    double detJ;
    FloatMatrix Gcov;
    FloatArray lcoords;
    lcoords = * gp->giveCoordinates();
    this->evalInitialCovarBaseVectorsAt(lcoords, Gcov);
    detJ = Gcov.giveDeterminant() * 0.5 * this->layeredCS->giveLayerThickness(layer);
    return detJ * gp->giveWeight();
}


void
Tr2Shell7PhFi:: compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer)
{
    int ndofs = 42;
    answer.resize(ndofs, ndofs);
    for ( int i = 1; i <= ndofs; i++ ) {
        for ( int j = 1; j <= 18; j++ ) {
            if ( fabs( matrix1.at(i, j) ) > 1.0e-12 ) {
                double diff = ( matrix1.at(i, j) - matrix2.at(i, j) );
                double relDiff =  diff / matrix1.at(i, j);
                if ( fabs(relDiff) < 1.0e-4 ) {
                    answer.at(i, j) = 0.0;
                } else if ( fabs(diff) < 1.0e3 ) {
                    answer.at(i, j) = 0.0;
                } else {
                    answer.at(i, j) = relDiff;
                }
            } else {
                answer.at(i, j) = -1.0;
            }
        }
    }
}

} // end namespace oofem
