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
#include <memory>

namespace oofem {
REGISTER_BoundaryCondition(SurfaceTensionBoundaryCondition);

IRResultType SurfaceTensionBoundaryCondition :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, this->gamma, _IFT_SurfaceTensionBoundaryCondition_gamma);

    this->useTangent = ir->hasField(_IFT_SurfaceTensionBoundaryCondition_useTangent);

    return ActiveBoundaryCondition :: initializeFrom(ir);
}

void SurfaceTensionBoundaryCondition :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                           const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    IntArray bNodes;

    rows.resize(boundaries.giveSize() / 2);
    cols.resize(boundaries.giveSize() / 2);

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        e->giveBoundaryLocationArray(rows [ pos ], bNodes, this->dofs, r_s);
        e->giveBoundaryLocationArray(cols [ pos ], bNodes, this->dofs, c_s);
    }
}

void SurfaceTensionBoundaryCondition :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                                 CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }

    OOFEM_ERROR("Not implemented yet.");

    FloatMatrix Ke;
    IntArray r_loc, c_loc, bNodes;

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        e->giveBoundaryLocationArray(r_loc, bNodes, this->dofs, r_s);
        e->giveBoundaryLocationArray(c_loc, bNodes, this->dofs, c_s);
        this->computeTangentFromElement(Ke, e, boundary, tStep);
        answer.assemble(r_loc, c_loc, Ke);
    }
}

void SurfaceTensionBoundaryCondition :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                       CharType type, ValueModeType mode,
                                                       const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    if ( type != ExternalForcesVector ) {
        return;
    }

    FloatArray fe;
    IntArray loc, masterdofids, bNodes;

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        e->giveBoundaryLocationArray(loc, bNodes, this->dofs, s, & masterdofids);
        this->computeLoadVectorFromElement(fe, e, boundary, tStep);
        answer.assemble(fe, loc);
        if ( eNorms ) {
            eNorms->assembleSquared(fe, masterdofids);
        }
    }
}

void SurfaceTensionBoundaryCondition :: computeTangentFromElement(FloatMatrix &answer, Element *e, int side, TimeStep *tStep)
{
    FEInterpolation *fei = e->giveInterpolation();
    if ( !fei ) {
        OOFEM_ERROR("No interpolation available for element.");
    }
    std :: unique_ptr< IntegrationRule > iRule( fei->giveBoundaryIntegrationRule(fei->giveInterpolationOrder()-1, side) );

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    int nodes = e->giveNumberOfNodes();
    if ( side == -1 ) {
        side = 1;
    }

    answer.clear();

    if ( nsd == 2 ) {
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

        if ( e->giveDomain()->isAxisymmetric() ) {
            FloatArray N;
            FloatArray gcoords;
            FloatArray tmpB(2 *nodes);
            for ( GaussPoint *gp: *iRule ) {
                fei2d->edgeEvaldNds( dNds, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryEvalN( N, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryLocal2Global( gcoords, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
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
            for ( GaussPoint *gp: *iRule ) {
                double t = e->giveCrossSection()->give(CS_Thickness, gp); ///@todo The thickness is not often relevant or used in FM.
                fei2d->edgeEvaldNds( dNds, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );

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

        FEInterpolation3d *fei3d = static_cast< FEInterpolation3d * >(fei);

        OOFEM_ERROR("3D tangents not implemented yet.");

        //FloatMatrix tmp(3 *nodes, 3 *nodes);
        FloatMatrix dNdx;
        FloatArray n;
        for ( GaussPoint *gp: *iRule ) {
            fei3d->surfaceEvaldNdx( dNdx, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
            /*double J = */ fei->boundaryEvalNormal( n, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
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
        OOFEM_WARNING("Only 2D or 3D is possible!");
    }
}

void SurfaceTensionBoundaryCondition :: computeLoadVectorFromElement(FloatArray &answer, Element *e, int side, TimeStep *tStep)
{
    FEInterpolation *fei = e->giveInterpolation();
    if ( !fei ) {
        OOFEM_ERROR("No interpolation or default integration available for element.");
    }
    std :: unique_ptr< IntegrationRule > iRule( fei->giveBoundaryIntegrationRule(fei->giveInterpolationOrder()-1, side) );

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();

    if ( side == -1 ) {
        side = 1;
    }

    answer.clear();

    if ( nsd == 2 ) {

        FEInterpolation2d *fei2d = static_cast< FEInterpolation2d * >(fei);

        ///@todo More of this grunt work should be moved to the interpolation classes
        IntArray bnodes;
        fei2d->boundaryGiveNodes(bnodes, side);
        int nodes = bnodes.giveSize();
        FloatMatrix xy(2, nodes);
        for ( int i = 1; i <= nodes; i++ ) {
            Node *node = e->giveNode(bnodes.at(i));
            ///@todo This should rely on the xindex and yindex in the interpolator..
            xy.at(1, i) = node->giveCoordinate(1);
            xy.at(2, i) = node->giveCoordinate(2);
        }

        FloatArray tmp; // Integrand
        FloatArray es; // Tangent vector to curve
        FloatArray dNds;

        if ( e->giveDomain()->isAxisymmetric() ) {
            FloatArray N;
            FloatArray gcoords;
            for ( GaussPoint *gp: *iRule ) {
                fei2d->edgeEvaldNds( dNds, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryEvalN( N, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                fei->boundaryLocal2Global( gcoords, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                double r = gcoords(0); // First coordinate is the radial coord.

                es.beProductOf(xy, dNds);

                tmp.resize( 2 * nodes);
                for ( int i = 0; i < nodes; i++ ) {
                    tmp(2 * i)   = dNds(i) * es(0) * r + N(i);
                    tmp(2 * i + 1) = dNds(i) * es(1) * r;
                }

                answer.add(- 2 * M_PI * gamma * J * gp->giveWeight(), tmp);
            }
        } else {
            for ( GaussPoint *gp: *iRule ) {
                double t = e->giveCrossSection()->give(CS_Thickness, gp);
                fei2d->edgeEvaldNds( dNds, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                double J = fei->boundaryGiveTransformationJacobian( side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
                es.beProductOf(xy, dNds);

                tmp.resize( 2 * nodes);
                for ( int i = 0; i < nodes; i++ ) {
                    tmp(2 * i)   = dNds(i) * es(0);
                    tmp(2 * i + 1) = dNds(i) * es(1);
                    //B.at(1, 1+i*2) = B.at(2, 2+i*2) = dNds(i);
                }
                //tmp.beTProductOf(B, es);

                answer.add(- t * gamma * J * gp->giveWeight(), tmp);
            }
        }
    } else if ( nsd ==  3 ) {

        FEInterpolation3d *fei3d = static_cast< FEInterpolation3d * >(fei);
        FloatArray n, surfProj;
        FloatMatrix dNdx, B;
        for ( GaussPoint *gp: *iRule ) {
            fei3d->surfaceEvaldNdx( dNdx, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );
            double J = fei->boundaryEvalNormal( n, side, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e) );

            // [I - n(x)n]  in voigt form:
            surfProj = {1. - n(0)*n(0), 1. - n(1)*n(1), 1. - n(2)*n(2),
                        - n(1)*n(2), - n(0)*n(2), - n(0)*n(1),
            };

            // Construct B matrix of the surface nodes
            B.resize(6, 3 * dNdx.giveNumberOfRows());
            for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
                B.at(1, 3 * i - 2) = dNdx.at(i, 1);
                B.at(2, 3 * i - 1) = dNdx.at(i, 2);
                B.at(3, 3 * i - 0) = dNdx.at(i, 3);

                B.at(5, 3 * i - 2) = B.at(4, 3 * i - 1) = dNdx.at(i, 3);
                B.at(6, 3 * i - 2) = B.at(4, 3 * i - 0) = dNdx.at(i, 2);
                B.at(6, 3 * i - 1) = B.at(5, 3 * i - 0) = dNdx.at(i, 1);
            }

            answer.plusProduct(B, surfProj, -gamma * J * gp->giveWeight() );
        }
    }
}
} // end namespace oofem
