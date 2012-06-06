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

#include "surfacetensionbc.h"
#include "elementside.h"
#include "element.h"
#include "node.h"
#include "crosssection.h"
#include "error.h"
#include "alist.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "sparsemtrx.h"

#include "integrationdomain.h"
#include "integrationrule.h"
#include "gausspnt.h"
#include "mathfem.h"

#include <utility>
#include <list>

namespace oofem {

IRResultType SurfaceTensionBoundaryCondition :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    result = ActiveBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->gamma, IFT_SurfaceTensionBoundaryCondition_gamma, "gamma");

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_SurfaceTensionBoundaryCondition_useTangent, "usetangent");
    this->useTangent = val;

    return IRRT_OK;
}

void SurfaceTensionBoundaryCondition :: giveLocationArrays(AList<IntArray> &rows, AList<IntArray> &cols, EquationID eid, CharType type,
                                const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain)
{
    if (!this->useTangent || type != TangentStiffnessMatrix)
        return;

    IntArray r_loc, c_loc, dofids;
    rows.growTo(this->elements.size() + this->sides.size());
    cols.growTo(this->elements.size() + this->sides.size());

    int i = 1;
    for (std :: list<int> :: const_iterator it = elements.begin(); it != elements.end(); ++it ) {
        Element *e = this->giveDomain()->giveElement(*it);
        e->giveLocationArray(r_loc, eid, r_s);
        e->giveLocationArray(c_loc, eid, c_s);
        rows.put(i, new IntArray(r_loc));
        cols.put(i, new IntArray(c_loc));
        i++;
    }
    for (std :: list< std::pair<int, int> > :: const_iterator it = sides.begin(); it != sides.end(); ++it ) {
        Element *e = this->giveDomain()->giveElement(it->first);
        int side = it->second;
        ElementSide *es = e->giveSide(side);
        e->giveDofManDofIDMask(1, eid, dofids); // NOTE! Assumes that the first node contains the relevant dof ID's.
        es->giveLocationArray(dofids, r_loc, r_s);
        es->giveLocationArray(dofids, c_loc, c_s);
        rows.put(i, new IntArray(r_loc));
        cols.put(i, new IntArray(c_loc));
        i++;
    }
}

void SurfaceTensionBoundaryCondition :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain)
{
    if (!this->useTangent || type != TangentStiffnessMatrix)
        return;

    OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assemble - Not implemented yet.");

    FloatMatrix Ke;
    IntArray r_loc, c_loc;
    IntArray dofids;

    for (std :: list<int> :: const_iterator it = elements.begin(); it != elements.end(); ++it ) {
        Element *e = this->giveDomain()->giveElement(*it);
        e->giveLocationArray(r_loc, eid, r_s); // Assumes the element has a sane/normal location vector.
        e->giveLocationArray(c_loc, eid, c_s);
        this->computeTangentFromElement(Ke, e, -1, tStep);
        answer->assemble(r_loc, c_loc, Ke);
    }
    for (std :: list< std::pair<int, int> > :: const_iterator it = sides.begin(); it != sides.end(); ++it ) {
        Element *e = this->giveDomain()->giveElement(it->first);
        int side = it->second;
        ElementSide *es = e->giveSide(side);
        e->giveDofManDofIDMask(1, eid, dofids); // NOTE! Assumes that the first node contains the relevant dof ID's.
        es->giveLocationArray(dofids, r_loc, r_s); // Assumes the element has a sane/normal location vector.
        es->giveLocationArray(dofids, c_loc, c_s);
        this->computeTangentFromElement(Ke, e, side, tStep);
        answer->assemble(r_loc, c_loc, Ke);
    }
}

double SurfaceTensionBoundaryCondition :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain)
{
    if (type != InternalForcesVector) {
        return 0.0;
    }
    if (eid == EID_MomentumBalance_ConservationEquation) {
        eid = EID_MomentumBalance;
    }
    if (eid != EID_MomentumBalance || mode != VM_Total) {
        return 0.0;
    }

    FloatArray fe;
    IntArray loc;
    IntArray dofids;
    double norm = 0.0;

    for (std :: list<int> :: const_iterator it = elements.begin(); it != elements.end(); ++it ) {
        Element *e = this->giveDomain()->giveElement(*it);
        e->giveLocationArray(loc, eid, s); // Assumes the element has a sane/normal location vector.
        this->computeLoadVectorFromElement(fe, e, -1, tStep);
        answer.assemble(fe, loc);
        norm += fe.computeSquaredNorm();
    }
    for (std :: list< std::pair<int, int> > :: const_iterator it = sides.begin(); it != sides.end(); ++it ) {
        Element *e = this->giveDomain()->giveElement(it->first);
        int side = it->second;
        ElementSide *es = e->giveSide(side);
        e->giveDofManDofIDMask(1, eid, dofids); // NOTE! Assumes that the first node contains the relevant dof ID's.
        es->giveLocationArray(dofids, loc, s); // Assumes the element has a sane/normal location vector.
        this->computeLoadVectorFromElement(fe, e, side, tStep);
        answer.assemble(fe, loc);
        norm += fe.computeSquaredNorm();
    }
    return norm;
}

void SurfaceTensionBoundaryCondition :: computeTangentFromElement(FloatMatrix &answer, Element *e, int side, TimeStep *tStep)
{
    FEInterpolation *fei = e->giveInterpolation();
    IntegrationRule *iRule = e->giveDefaultIntegrationRulePtr();
    if (!fei || !iRule)
        OOFEM_ERROR("SurfaceTensionBoundaryCondition :: computeTangentFromElement - No interpolation available for element.");

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    integrationDomain id = e->giveIntegrationDomain();
    int nodes = e->giveNumberOfNodes();

    if (nsd == 2) {
        if (!((side == -1 && id == _Line) || (side > 0 && (id == _Triangle || id == _Square))))
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Not a surface element.");
        if (side == -1)
            side = 1;

        FEInterpolation2d *fei2d = static_cast<FEInterpolation2d*>(fei);

        ///@todo More of this grunt work should be moved to the interpolation classes
        FloatMatrix xy(2, nodes);
        Node *node;
        for ( int i = 1; i <= nodes; i++ ) {
            node = e->giveNode(i);
            xy.at(1, i) = node->giveCoordinate(1);
            xy.at(2, i) = node->giveCoordinate(2);
        }

        FloatArray tmpA(2*nodes);
        FloatArray es; // Tangent vector to curve
        FloatArray dNds;
        FloatMatrix B(2,2*nodes);
        B.zero();
        answer.resize(2*nodes, 2*nodes);
        answer.zero();

        if (e->giveDomain()->giveDomainType() == _3dAxisymmMode) {
            FloatArray N;
            FloatArray gcoords;
            FloatArray tmpB(2*nodes);
            for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                fei2d->edgeEvaldNds(dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                fei2d->edgeEvalN(N, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                double J = fei2d->edgeGiveTransformationJacobian(side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                fei2d->edgeLocal2global(gcoords, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                double r = gcoords(0); // First coordinate is the radial coord.

                es.beProductOf(xy, dNds);

                // Construct the different matrices in the integrand;
                for (int i = 0; i < nodes; i++) {
                    tmpA(i*2+0) = dNds(i)*es(0);
                    tmpA(i*2+1) = dNds(i)*es(1);
                    tmpB(i*2+0) = N(i);
                    tmpB(i*2+1) = 0;
                    B(i*2,0) = B(i*2+1, 1) = dNds(i);
                }

                double dV = 2 * M_PI * gamma * J * gp->giveWeight();
                answer.plusDyadSymmUpper(tmpA, tmpB, dV);
                answer.plusDyadSymmUpper(tmpB, tmpA, dV);
                answer.plusProductSymmUpper(B, B, r*dV);
                answer.plusDyadSymmUpper(tmpA, tmpA, -r*dV);
            }

        } else {
            double t = e->giveCrossSection()->give(CS_Thickness); ///@todo The thickness is not often relevant or used in FM.
            for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                fei2d->edgeEvaldNds(dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                double J = fei2d->edgeGiveTransformationJacobian(side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));

                es.beProductOf(xy, dNds);

                // Construct the different matrices in the integrand;
                for (int i = 0; i < nodes; i++) {
                    tmpA(i*2+0) = dNds(i)*es(0);
                    tmpA(i*2+1) = dNds(i)*es(1);
                    B(i*2,0) = B(i*2+1, 1) = dNds(i);
                }

                double dV = t * gamma * J * gp->giveWeight();
                answer.plusProductSymmUpper(B, B, dV);
                answer.plusDyadSymmUpper(tmpA, tmpA, -dV);
            }
        }

        answer.symmetrized();

    }  else if (nsd ==  3) {
        if ( !(  ((id == _Triangle || id == _Square) && side == -1) || ((id == _Tetrahedra || id == _Cube) && side > 0)  ) )
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Not a surface element.");
        if (side == -1)
            side = 1;

        FEInterpolation3d *fei3d = static_cast<FEInterpolation3d*>(fei);

        OOFEM_ERROR("SurfaceTensionBoundaryCondition :: computeTangentFromElement - 3D tangents not implemented yet.");

        FloatMatrix tmp(3*nodes, 3*nodes);
        FloatMatrix dNdx;
        FloatArray n;
        answer.resize(3*nodes, 3*nodes);
        answer.zero();
        for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(k);
            fei3d->surfaceEvaldNdx(dNdx, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
            fei3d->surfaceEvalNormal(n, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
            //double J = fei3d->surfaceGiveTransformationJacobian(side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e), 0.0);
            //double dV = gamma * J * gp->giveWeight();

            for (int i = 0; i < nodes; i++) {
                //tmp(3*i+0) = dNdx(i,0) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(0);
                //tmp(3*i+1) = dNdx(i,1) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(1);
                //tmp(3*i+2) = dNdx(i,2) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(2);
            }
            //answer.plusProductSymmUpper(A,B, dV);
            ///@todo  Derive expressions for this.
        }
        OOFEM_WARNING("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - 3D Completely untested!");
    }
}

void SurfaceTensionBoundaryCondition :: computeLoadVectorFromElement(FloatArray &answer, Element *e, int side, TimeStep *tStep)
{
    FEInterpolation *fei = e->giveInterpolation();
    IntegrationRule *iRule = e->giveDefaultIntegrationRulePtr();
    if (!fei || !iRule)
        OOFEM_ERROR("SurfaceTensionBoundaryCondition :: computeLoadVectorFromElement - No interpolation or default integration available for element.");

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    integrationDomain id = iRule->giveIntegrationDomain();
    int nodes = e->giveNumberOfNodes();

    if (nsd == 2) {
        if ( !((side == -1 && id == _Line) || (side > 0 && (id == _Triangle || id == _Square))) )
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Not a surface element.");
        if (side == -1)
            side = 1;

        FEInterpolation2d *fei2d = static_cast<FEInterpolation2d*>(fei);

        ///@todo More of this grunt work should be moved to the interpolation classes
        FloatMatrix xy(2, nodes);
        Node *node;
        for ( int i = 1; i <= nodes; i++ ) {
            node = e->giveNode(i);
            xy.at(1, i) = node->giveCoordinate(1);
            xy.at(2, i) = node->giveCoordinate(2);
        }

        FloatArray tmp(2*nodes); // Integrand
        FloatArray es; // Tangent vector to curve
        FloatArray dNds;
        answer.resize(2*nodes);
        answer.zero();

        if (e->giveDomain()->giveDomainType() == _3dAxisymmMode) {
            FloatArray N;
            FloatArray gcoords;
            for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                fei2d->edgeEvaldNds(dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                fei2d->edgeEvalN(N, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                double J = fei2d->edgeGiveTransformationJacobian(side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                fei2d->edgeLocal2global(gcoords, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                double r = gcoords(0); // First coordinate is the radial coord.

                es.beProductOf(xy, dNds);
                for (int i = 0; i < nodes; i++) {
                    tmp(2*i)   = dNds(i)*es(0)*r + N(i);
                    tmp(2*i+1) = dNds(i)*es(1)*r;
                }

                double dA = 2 * M_PI * gamma * J * gp->giveWeight();
                answer.add( dA, tmp);
            }

        } else {
            double t = e->giveCrossSection()->give(CS_Thickness);
            for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(k);
                fei2d->edgeEvaldNds(dNds, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                double J = fei2d->edgeGiveTransformationJacobian(side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
                es.beProductOf(xy, dNds);

                for (int i = 0; i < nodes; i++) {
                    tmp(2*i)   = dNds(i)*es(0);
                    tmp(2*i+1) = dNds(i)*es(1);
                    //B.at(1, 1+i*2) = B.at(2, 2+i*2) = dNds(i);
                }
                //tmp.beTProductOf(B, es);

                double dA = t * gamma * J * gp->giveWeight();
                answer.add( dA, tmp);
            }
        }

    } else if (nsd ==  3) {
        if ( !(  ((id == _Triangle || id == _Square) && side == -1) || ((id == _Tetrahedra || id == _Cube) && side > 0)  ) )
            OOFEM_ERROR("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - Not a surface element.");
        if (side == -1)
            side = 1;

        FEInterpolation3d *fei3d = static_cast<FEInterpolation3d*>(fei);
        FloatArray tmp(3*nodes);
        FloatMatrix dNdx;
        FloatArray n;
        answer.resize(3*nodes);
        answer.zero();
        for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(k);
            fei3d->surfaceEvaldNdx(dNdx, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
            fei3d->surfaceEvalNormal(n, side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));
            double J = fei3d->surfaceGiveTransformationJacobian(side, * gp->giveCoordinates(), FEIElementGeometryWrapper(e));

            for (int i = 0; i < nodes; i++) {
                tmp(3*i+0) = dNdx(i,0) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(0);
                tmp(3*i+1) = dNdx(i,1) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(1);
                tmp(3*i+2) = dNdx(i,2) - (dNdx(i,0)*n(0)* + dNdx(i,1)*n(1) + dNdx(i,2)*n(2))*n(2);
            }
            answer.add( - gamma * J * gp->giveWeight(), tmp);
        }
        OOFEM_WARNING("SurfaceTensionBoundaryCondition :: assembleVectorFromElement - 3D Completely untested!");
    }
}

} // end namespace oofem
