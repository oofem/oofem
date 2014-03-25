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

#include "surfacetensionbc.h"
#include "element.h"
#include "node.h"
#include "crosssection.h"
#include "error.h"
#include "alist.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "classfactory.h"
#include "set.h"
#include "sparsemtrx.h"

#include "integrationdomain.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "mathfem.h"

#include <utility>
#include <list>

namespace oofem {
REGISTER_BoundaryCondition(SurfaceTensionBoundaryCondition);

IRResultType SurfaceTensionBoundaryCondition :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    IR_GIVE_FIELD(ir, this->gamma, _IFT_SurfaceTensionBoundaryCondition_gamma);

    this->useTangent = ir->hasField(_IFT_SurfaceTensionBoundaryCondition_useTangent);

    return ActiveBoundaryCondition :: initializeFrom(ir);
}

void SurfaceTensionBoundaryCondition :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type,
                                                           const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        eid = EID_MomentumBalance;
    }
    if ( !this->useTangent || eid != EID_MomentumBalance || type != TangentStiffnessMatrix ) {
        return;
    }

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    IntArray bNodes;

    rows.resize(boundaries.giveSize() / 2);
    cols.resize(boundaries.giveSize() / 2);

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        ElementGeometry *eG = this->giveDomain()->giveElementGeometry( boundaries.at(pos * 2 - 1) );
		ElementEvaluator *eE = this->giveDomain()->giveElementEvaluator(boundaries.at(pos * 2 - 1));

        int boundary = boundaries.at(pos * 2);

        eG->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        eE->giveBoundaryLocationArray(rows [ pos ], bNodes, eid, r_s, eG);
        eE->giveBoundaryLocationArray(cols [ pos ], bNodes, eid, c_s, eG);
    }
}

void SurfaceTensionBoundaryCondition :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                                                 CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        eid = EID_MomentumBalance;
    }
    if ( !this->useTangent || eid != EID_MomentumBalance || type != TangentStiffnessMatrix ) {
        return;
    }

    OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assemble - Not implemented yet.");

    FloatMatrix Ke;
    IntArray r_loc, c_loc, bNodes;

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        ElementGeometry *eG = this->giveDomain()->giveElementGeometry( boundaries.at(pos * 2 - 1) );
		ElementEvaluator *eE = this->giveDomain()->giveElementEvaluator(boundaries.at(pos * 2 - 1));

        int boundary = boundaries.at(pos * 2);

        eG->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        eE->giveBoundaryLocationArray(r_loc, bNodes, eid, r_s, eG);
        eE->giveBoundaryLocationArray(c_loc, bNodes, eid, c_s, eG);
        this->computeTangentFromElement(Ke, eG, boundary, tStep);
        answer->assemble(r_loc, c_loc, Ke);
    }
}

void SurfaceTensionBoundaryCondition :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                                       CharType type, ValueModeType mode,
                                                       const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    if ( type != InternalForcesVector ) {
        return;
    }
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        eid = EID_MomentumBalance;
    }
    if ( eid != EID_MomentumBalance || mode != VM_Total ) {
        return;
    }

    FloatArray fe;
    IntArray loc, dofids, masterdofids, bNodes;

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        ElementGeometry *eG = this->giveDomain()->giveElementGeometry( boundaries.at(pos * 2 - 1) );
		ElementEvaluator *eE = this->giveDomain()->giveElementEvaluator(boundaries.at(pos * 2 - 1));
		
        int boundary = boundaries.at(pos * 2);

        eG->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

		eE->giveBoundaryLocationArray(loc, bNodes, eid, s, eG, & masterdofids);
        this->computeLoadVectorFromElement(fe, eG, boundary, tStep);
        answer.assemble(fe, loc);
        if ( eNorms ) {
            for ( int i = 1; i <= loc.giveSize(); ++i ) {
                if ( loc.at(i) > 0 ) {
                    eNorms->at( masterdofids.at(1) ) += fe.at(i) * fe.at(i);
                }
            }
        }
    }
}

void SurfaceTensionBoundaryCondition :: computeTangentFromElement(FloatMatrix &answer, ElementGeometry *e, int side, TimeStep *tStep)
{
    FEInterpolation *fei = e->giveInterpolation();
    IntegrationRule *iRule = e->giveDefaultIntegrationRulePtr();
    if ( !fei || !iRule ) {
        OOFEM_ERROR("SurfaceTensionBoundaryCondition :: computeTangentFromElement - No interpolation available for element.");
    }

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    integrationDomain id = e->giveIntegrationDomain();
    int nodes = e->giveNumberOfNodes();

    if ( nsd == 2 ) {
        if ( !( ( side == -1 && id == _Line ) || ( side > 0 && ( id == _Triangle || id == _Square ) ) ) ) {
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Not a surface element.");
        }
        if ( side == -1 ) {
            side = 1;
        }

        FEInterpolation2d *fei2d = static_cast< FEInterpolation2d * >(fei);

        ///@todo More of this grunt work should be moved to the interpolation classes
        FloatMatrix xy(2, nodes);
        Node *node;
        for ( int i = 1; i <= nodes; i++ ) {
            node = e->giveNode(i);
            xy.at(1, i) = node->giveCoordinate(1);
            xy.at(2, i) = node->giveCoordinate(2);
        }

        FloatArray tmpA(2 *nodes);
        FloatArray es; // Tangent vector to curve
        FloatArray dNds;
        FloatMatrix B(2, 2 *nodes);
        B.zero();
        answer.resize(2 * nodes, 2 * nodes);
        answer.zero();

        if ( e->giveDomain()->isAxisymmetric() ) {
            FloatArray N;
            FloatArray gcoords;
            FloatArray tmpB(2 *nodes);
            for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                fei2d->edgeEvaldNds( dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryEvalN( N, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryLocal2Global( gcoords, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                double r = gcoords(0); // First coordinate is the radial coord.

                es.beProductOf(xy, dNds);

                // Construct the different matrices in the integrand;
                for ( int i = 0; i < nodes; i++ ) {
                    tmpA(i * 2 + 0) = dNds(i) * es(0);
                    tmpA(i * 2 + 1) = dNds(i) * es(1);
                    tmpB(i * 2 + 0) = N(i);
                    tmpB(i * 2 + 1) = 0;
                    B(i * 2, 0) = B(i * 2 + 1, 1) = dNds(i);
                }

                double dV = 2 *M_PI *gamma *J *gp->giveWeight();
                answer.plusDyadUnsym(tmpA, tmpB, dV);
                answer.plusDyadUnsym(tmpB, tmpA, dV);
                answer.plusProductSymmUpper(B, B, r * dV);
                answer.plusDyadUnsym(tmpA, tmpA, -r * dV);
            }
        } else {
            for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                double t = e->giveCrossSection()->give(CS_Thickness, gp); ///@todo The thickness is not often relevant or used in FM.
                fei2d->edgeEvaldNds( dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );

                es.beProductOf(xy, dNds);

                // Construct the different matrices in the integrand;
                for ( int i = 0; i < nodes; i++ ) {
                    tmpA(i * 2 + 0) = dNds(i) * es(0);
                    tmpA(i * 2 + 1) = dNds(i) * es(1);
                    B(i * 2, 0) = B(i * 2 + 1, 1) = dNds(i);
                }

                double dV = t * gamma * J * gp->giveWeight();
                answer.plusProductSymmUpper(B, B, dV);
                answer.plusDyadSymmUpper(tmpA, -dV);
            }
        }

        answer.symmetrized();
    }  else if ( nsd ==  3 ) {
        if ( !( ( ( id == _Triangle || id == _Square ) && side == -1 ) || ( ( id == _Tetrahedra || id == _Cube ) && side > 0 ) ) ) {
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Not a surface element.");
        }
        if ( side == -1 ) {
            side = 1;
        }

        FEInterpolation3d *fei3d = static_cast< FEInterpolation3d * >(fei);

        OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - 3D tangents not implemented yet.");

        FloatMatrix tmp(3 *nodes, 3 *nodes);
        FloatMatrix dNdx;
        FloatArray n;
        answer.resize(3 * nodes, 3 * nodes);
        answer.zero();
        for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(k);
            fei3d->surfaceEvaldNdx( dNdx, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
            /*double J = */ fei->boundaryEvalNormal( n, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
            //double dV = gamma * J * gp->giveWeight();

            for ( int i = 0; i < nodes; i++ ) {
                //tmp(3*i+0) = dNdx(i,0) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(0);
                //tmp(3*i+1) = dNdx(i,1) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(1);
                //tmp(3*i+2) = dNdx(i,2) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(2);
            }
            //answer.plusProductSymmUpper(A,B, dV);
            ///@todo  Derive expressions for this.
        }
    } else {
        OOFEM_WARNING("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Only 2D or 3D is possible!");
    }
}

void SurfaceTensionBoundaryCondition :: computeLoadVectorFromElement(FloatArray &answer, ElementGeometry *e, int side, TimeStep *tStep)
{
    FEInterpolation *fei = e->giveInterpolation();
    IntegrationRule *iRule = e->giveDefaultIntegrationRulePtr();
    if ( !fei || !iRule ) {
        OOFEM_ERROR("SurfaceTensionBoundaryCondition :: computeLoadVectorFromElement - No interpolation or default integration available for element.");
    }

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    integrationDomain id = iRule->giveIntegrationDomain();
    int nodes = e->giveNumberOfNodes();

    if ( nsd == 2 ) {
        if ( !( ( side == -1 && id == _Line ) || ( side > 0 && ( id == _Triangle || id == _Square ) ) ) ) {
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: computeLoadVectorFromElement - Not a surface element.");
        }
        if ( side == -1 ) {
            side = 1;
        }

        FEInterpolation2d *fei2d = static_cast< FEInterpolation2d * >(fei);

        ///@todo More of this grunt work should be moved to the interpolation classes
        FloatMatrix xy(2, nodes);
        Node *node;
        for ( int i = 1; i <= nodes; i++ ) {
            node = e->giveNode(i);
            xy.at(1, i) = node->giveCoordinate(1);
            xy.at(2, i) = node->giveCoordinate(2);
        }

        FloatArray tmp(2 *nodes); // Integrand
        FloatArray es; // Tangent vector to curve
        FloatArray dNds;
        answer.resize(2 * nodes);
        answer.zero();

        if ( e->giveDomain()->isAxisymmetric() ) {
            FloatArray N;
            FloatArray gcoords;
            for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                fei2d->edgeEvaldNds( dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryEvalN( N, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryLocal2Global( gcoords, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                double r = gcoords(0); // First coordinate is the radial coord.

                es.beProductOf(xy, dNds);
                for ( int i = 0; i < nodes; i++ ) {
                    tmp(2 * i)   = dNds(i) * es(0) * r + N(i);
                    tmp(2 * i + 1) = dNds(i) * es(1) * r;
                }

                double dA = 2 *M_PI *gamma *J *gp->giveWeight();
                answer.add(dA, tmp);
            }
        } else {
            for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                double t = e->giveCrossSection()->give(CS_Thickness, gp);
                fei2d->edgeEvaldNds( dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
                es.beProductOf(xy, dNds);

                for ( int i = 0; i < nodes; i++ ) {
                    tmp(2 * i)   = dNds(i) * es(0);
                    tmp(2 * i + 1) = dNds(i) * es(1);
                    //B.at(1, 1+i*2) = B.at(2, 2+i*2) = dNds(i);
                }
                //tmp.beTProductOf(B, es);

                double dA = t * gamma * J * gp->giveWeight();
                answer.add(dA, tmp);
            }
        }
    } else if ( nsd ==  3 ) {
        if ( !( ( ( id == _Triangle || id == _Square ) && side == -1 ) || ( ( id == _Tetrahedra || id == _Cube ) && side > 0 ) ) ) {
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Not a surface element.");
        }
        if ( side == -1 ) {
            side = 1;
        }

        FEInterpolation3d *fei3d = static_cast< FEInterpolation3d * >(fei);
        FloatArray tmp(3 *nodes);
        FloatMatrix dNdx;
        FloatArray n;
        answer.resize(3 * nodes);
        answer.zero();
        for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(k);
            fei3d->surfaceEvaldNdx( dNdx, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );
            double J = fei->boundaryEvalNormal( n, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e) );

            for ( int i = 0; i < nodes; i++ ) {
                tmp(3 * i + 0) = dNdx(i, 0) - ( dNdx(i, 0) * n(0) * +dNdx(i, 1) * n(1) + dNdx(i, 2) * n(2) ) * n(0);
                tmp(3 * i + 1) = dNdx(i, 1) - ( dNdx(i, 0) * n(0) * +dNdx(i, 1) * n(1) + dNdx(i, 2) * n(2) ) * n(1);
                tmp(3 * i + 2) = dNdx(i, 2) - ( dNdx(i, 0) * n(0) * +dNdx(i, 1) * n(1) + dNdx(i, 2) * n(2) ) * n(2);
            }
            answer.add(-gamma * J * gp->giveWeight(), tmp);
        }
        OOFEM_WARNING("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - 3D Completely untested!");
    }
}
} // end namespace oofem
