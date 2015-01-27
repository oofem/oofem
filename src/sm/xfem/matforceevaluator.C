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

#include "matforceevaluator.h"
#include "xfem/tipinfo.h"
#include "domain.h"
#include "spatiallocalizer.h"
#include "Elements/nlstructuralelement.h"
#include "dofmanager.h"
#include "integrationrule.h"
#include "feinterpol.h"
#include "gausspoint.h"
#include "CrossSections/structuralcrosssection.h"

#include <fstream>
#include <set>

namespace oofem {

MaterialForceEvaluator::MaterialForceEvaluator() {

}

MaterialForceEvaluator::~MaterialForceEvaluator() {

}

void MaterialForceEvaluator::computeMaterialForce(FloatArray &oMatForce, Domain &iDomain, const TipInfo &iTipInfo, TimeStep *tStep, const double &iRadius)
{
    oMatForce.resize(2);
    oMatForce.clear();


    // Fetch elements within a specified volume around the crack tip
    SpatialLocalizer *localizer = iDomain.giveSpatialLocalizer();
    std :: set< int > elList;
    localizer->giveAllElementsWithNodesWithinBox(elList, iTipInfo.mGlobalCoord, iRadius);

    if( elList.empty() ) {
        FloatArray lcoords, closest;
        Element *closestEl = localizer->giveElementClosestToPoint(lcoords, closest, iTipInfo.mGlobalCoord);

        if(closestEl != NULL) {
            if( closest.distance(iTipInfo.mGlobalCoord) < 1.0e-9 ) {
                int elInd = closestEl->giveGlobalNumber();
                int elPlaceInArray = iDomain.giveElementPlaceInArray(elInd);
                elList.insert( elPlaceInArray );
            }
        }
        else {
            printf("Could not find closest element.\n");
        }
    }

    // Compute phi in nodes
    std::unordered_map<int, double> weightInNodes;

    for( int elIndex : elList ) {
//        int elPlaceInArray = iDomain.giveElementPlaceInArray(elIndex);
        Element *el = iDomain.giveElement(elIndex);

        const IntArray &elNodes = el->giveDofManArray();
        for(int nodeInd : elNodes) {
            DofManager *dMan = iDomain.giveDofManager(nodeInd);

            weightInNodes[nodeInd] = computeWeightFunctionInPoint( *(dMan->giveCoordinates()), iTipInfo.mGlobalCoord, iRadius);
        }

    }

    // Loop over elements and Gauss points
    for( int elIndex : elList ) {
//        int elPlaceInArray = iDomain.giveElementPlaceInArray(elIndex);
        NLStructuralElement *el = dynamic_cast<NLStructuralElement*>( iDomain.giveElement(elIndex) );

        if(el == NULL) {
            OOFEM_ERROR("Could not cast to NLStructuralElement.")
        }

//        const IntArray &elNodes = el->giveDofManArray();

        for(GaussPoint *gp : *(el->giveDefaultIntegrationRulePtr()) ) {

            FloatArray gradWeight;
            const FloatArray &pos = gp->giveGlobalCoordinates();

#if 0
            // Compute grad(phi)
            FloatMatrix dNdx;
            FEInterpolation *interp = el->giveInterpolation();
            const FEIElementGeometryWrapper geomWrapper(el);
            interp->evaldNdx(dNdx, gp->giveNaturalCoordinates(), geomWrapper);

            FloatArray weightInElNodes;

            weightInElNodes.resize( elNodes.giveSize() );
            for(int i = 0; i < weightInElNodes.giveSize(); i++) {
                weightInElNodes[i] = weightInNodes[elNodes[i]];
            }

            gradWeight.beTProductOf(dNdx, weightInElNodes);
#else
            FloatArray q;
//            q.beDifferenceOf(pos, iTipInfo.mGlobalCoord);
//            printf("q: "); q.printYourself();
            q = {pos[0]-iTipInfo.mGlobalCoord[0],pos[1]-iTipInfo.mGlobalCoord[1]};

            gradWeight.beScaled(-1.0/(iRadius*q.computeNorm()),q);
#endif

            // Compute Eshelby stress
            StructuralCrossSection *cs = dynamic_cast<StructuralCrossSection*>( gp->giveCrossSection() );
            if(cs != NULL) {

                FloatArray Fv;
                el->computeDeformationGradientVector(Fv, gp, tStep);

                FloatArray eshelbyStressV;
                cs->giveEshelbyStresses(eshelbyStressV, gp, Fv, tStep);

                FloatArray fullEshelbyStressV;
                StructuralMaterial :: giveFullVectorForm( fullEshelbyStressV, eshelbyStressV, _PlaneStrain );

                FloatMatrix eshelbyStress3D;
                eshelbyStress3D.beMatrixForm(fullEshelbyStressV);

                FloatMatrix eshelbyStress2D;
                eshelbyStress2D.beSubMatrixOf(eshelbyStress3D, 1, 2, 1, 2);


                if( computeWeightFunctionInPoint( pos, iTipInfo.mGlobalCoord, iRadius) > 0.0 ) {
                    // Add contribution
                    FloatArray contrib(2);
                    contrib.beProductOf(eshelbyStress2D, gradWeight);
                    contrib.times( -el->computeVolumeAround(gp) );

//                    if( contrib.computeNorm() > 1.0e10 ) {
//                        printf("Fv: "); Fv.printYourself();
//                    }

                    oMatForce.add(contrib);
                }
            }
        }

    }

//    printf("oMatForce: "); oMatForce.printYourself();

}

double MaterialForceEvaluator::computeWeightFunctionInPoint(const FloatArray &iCoord, const FloatArray &iTipCoord, const double &iRadius) const
{
    double weight = 0.0;

    double r = iTipCoord.distance(iCoord);

//    if(r <= iRadius) {
        weight = 1.0 - r/iRadius;
//    }

//    printf("r: %e weight: %e\n", r, weight);

    return weight;
}

} /* namespace oofem */
