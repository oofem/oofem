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

#include "transportelement.h"
#include "domain.h"
#include "transportmaterial.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "verbose.h"
#include "mathfem.h"
#include "crosssection.h"
#include "transportcrosssection.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "dof.h"
#include "stationarytransportproblem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

const double TransportElement :: stefanBoltzmann = 5.67e-8; //W/m2/K4

TransportElement :: TransportElement(int n, Domain *aDomain, ElementMode em) :
    Element(n, aDomain), emode( em )
{
}


TransportElement :: ~TransportElement()
{ }


void
TransportElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( emode == HeatTransferEM ) {
        answer = {T_f};
    } else if ( emode == HeatMass1TransferEM ) {
        answer = {T_f, C_1};
    } else if ( emode == Mass1TransferEM ) {
        answer = {C_1};
    } else {
        OOFEM_ERROR("Unknown ElementMode");
    }
}


void
TransportElement :: giveCharacteristicMatrix(FloatMatrix &answer,
                                             CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == TangentStiffnessMatrix ) {
        FloatMatrix tmp;
        this->computeConductivityMatrix(answer, Conductivity, tStep);
        this->computeBCMtrxAt(tmp, tStep, VM_Total);
        answer.add(tmp);
    } else if ( mtrx == ConductivityMatrix ) {
        this->computeConductivityMatrix(answer, Conductivity, tStep);
    } else if ( mtrx == CapacityMatrix || mtrx == MassMatrix ) {
        this->computeCapacityMatrix(answer, tStep);
    } else if ( mtrx == LumpedMassMatrix ) {
        this->computeCapacityMatrix(answer, tStep);
        for ( int i = 1; i <= answer.giveNumberOfRows(); i++ ) {
            double s = 0.0;
            for ( int j = 1; j <= answer.giveNumberOfColumns(); j++ ) {
                s += answer.at(i, j);
                answer.at(i, j) = 0.0;
            }
            answer.at(i, i) = s;
        }
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx));
    }
}


void
TransportElement :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                             TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == InternalForcesVector ) {
        this->computeInternalForcesVector(answer, tStep);
    } else if ( mtrx == ExternalForcesVector ) {
        this->computeExternalForcesVector(answer, tStep, mode);
    } else if ( mtrx == InertiaForcesVector ) {
        this->computeInertiaForcesVector(answer, tStep);
    } else if ( mtrx == LumpedMassMatrix ) {
        this->computeLumpedCapacityVector(answer, tStep);
    } else {
        OOFEM_ERROR( "Unknown Type of characteristic mtrx (%s)",
                __CharTypeToString(mtrx) );
    }
}


int
TransportElement :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    if ( this->material > 0 && !dynamic_cast< TransportMaterial* >( this->giveMaterial() ) ) {
        OOFEM_WARNING("cross section without support for transport problems");
        result = 0;
    }

#if 0
    if ( !this->giveCrossSection()->testCrossSectionExtension(CS_TransportCapability) ) {
        OOFEM_WARNING("cross-section without support for transport problems", 1);
        result = 0;
    }

#endif
    return result;
}

void
TransportElement :: computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize( this->computeNumberOfDofs(), this->computeNumberOfDofs() );
    answer.zero();

    if ( emode == HeatTransferEM || emode == Mass1TransferEM ) {
        this->computeCapacitySubMatrix(answer, Capacity, 0, tStep);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatMatrix subAnswer;
        MatResponseMode rmode [ 2 ] = {
            Capacity_hh, Capacity_ww
        };
        //double coeff = 1.0; //this->giveMaterial()->give('d');

        for ( int i = 1; i <= 2; i++ ) {
            this->computeCapacitySubMatrix(subAnswer, rmode [ i - 1 ], 0, tStep);
            this->assembleLocalContribution(answer, subAnswer, 2, i, i);
        }
    } else {
        OOFEM_ERROR("Unknown ElementMode");
    }
}

void
TransportElement :: computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    answer.resize( this->computeNumberOfDofs(), this->computeNumberOfDofs() );
    answer.zero();
    if ( emode == HeatTransferEM || emode == Mass1TransferEM ) {
        this->computeConductivitySubMatrix(answer, 0, Conductivity_hh, tStep);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatMatrix subAnswer;
        MatResponseMode rmode [ 2 ] [ 2 ] = { { Conductivity_hh, Conductivity_hw }, { Conductivity_wh, Conductivity_ww } };

        for ( int i = 1; i <= 2; i++ ) {
            for ( int j = 1; j <= 2; j++ ) {
                this->computeConductivitySubMatrix(subAnswer, 0, rmode [ i - 1 ] [ j - 1 ], tStep);
                this->assembleLocalContribution(answer, subAnswer, 2, i, j);
            }
        }
    } else {
        OOFEM_ERROR("Unknown ElementMode");
    }
}


void
TransportElement :: computeNAt(FloatArray &answer, const FloatArray &lcoord)
{
    this->giveInterpolation()->evalN( answer, lcoord, FEIElementGeometryWrapper(this) );
}


void
TransportElement :: computeNmatrixAt(FloatMatrix &answer, const FloatArray &lcoords)
{
    FloatArray n;
    this->computeNAt(n, lcoords);
    if ( this->emode == HeatTransferEM  || this->emode == Mass1TransferEM ) {
        answer.beNMatrixOf(n, 1);
    } else {
        answer.beNMatrixOf(n, 2);
    }
}


void
TransportElement :: computeBmatrixAt(FloatMatrix &answer, const FloatArray &lcoords)
{
    FloatMatrix dnx;
    ///@todo We should change the transposition in evaldNdx;
    this->giveInterpolation()->evaldNdx( dnx, lcoords, FEIElementGeometryWrapper(this) );
    if ( emode == HeatTransferEM || emode == Mass1TransferEM ) {
        answer.beTranspositionOf(dnx);
    } else if ( this->emode == HeatMass1TransferEM ) {
        int nodes = dnx.giveNumberOfRows();
        int nsd = dnx.giveNumberOfColumns();
        answer.resize(nsd*2, nodes*2);
        answer.zero();
        for (int i = 0; i < nodes; ++i) {
            for (int j = 0; j < nsd; ++j) {
                answer(j, i*2) = dnx(i, j);
                answer(j+nsd, i*2+1) = dnx(i, j);
            }
        }
    }
}


void
TransportElement :: computeGradientMatrixAt(FloatMatrix &answer, const FloatArray &lcoords)
{
    FloatMatrix dnx;
    this->giveInterpolation()->evaldNdx( dnx, lcoords, FEIElementGeometryWrapper(this) );
    ///@todo We should change the transposition in evaldNdx;
    answer.beTranspositionOf(dnx);
}


void
TransportElement :: computeEgdeNAt(FloatArray &answer, int iedge, const FloatArray &lcoords)
{
    FEInterpolation *interp = this->giveInterpolation();
    if ( dynamic_cast< FEInterpolation2d * >(interp) ) {
        dynamic_cast< FEInterpolation2d * >(interp)->edgeEvalN( answer, iedge, lcoords, FEIElementGeometryWrapper(this) );
    } else if ( dynamic_cast< FEInterpolation3d * >(interp) ) {
        dynamic_cast< FEInterpolation3d * >(interp)->edgeEvalN( answer, iedge, lcoords, FEIElementGeometryWrapper(this) );
    }
}


void
TransportElement :: giveEdgeDofMapping(IntArray &answer, int iEdge)
{
    FEInterpolation *interp = this->giveInterpolation();
    if ( dynamic_cast< FEInterpolation2d * >(interp) ) {
        dynamic_cast< FEInterpolation2d * >(interp)->computeLocalEdgeMapping(answer, iEdge);
    } else if ( dynamic_cast< FEInterpolation3d * >(interp) ) {
        dynamic_cast< FEInterpolation3d * >(interp)->computeLocalEdgeMapping(answer, iEdge);
    }
}


void
TransportElement :: computeEdgeIpGlobalCoords(FloatArray &answer, const FloatArray &lcoords, int iEdge)
{
    FEInterpolation *interp = this->giveInterpolation();
    if ( dynamic_cast< FEInterpolation2d * >(interp) ) {
        dynamic_cast< FEInterpolation2d * >(interp)->edgeLocal2global( answer, iEdge, lcoords, FEIElementGeometryWrapper(this) );
    } else if ( dynamic_cast< FEInterpolation3d * >(interp) ) {
        dynamic_cast< FEInterpolation3d * >(interp)->edgeLocal2global( answer, iEdge, lcoords, FEIElementGeometryWrapper(this) );
    }
}

void
TransportElement :: computeSurfaceNAt(FloatArray &answer, int iSurf, const FloatArray &lcoord)
{
    FEInterpolation *interp = this->giveInterpolation();
    if ( dynamic_cast< FEInterpolation3d * >(interp) ) {
        dynamic_cast< FEInterpolation3d * >(interp)->surfaceEvalN( answer, iSurf, lcoord, FEIElementGeometryWrapper(this) );
    }
}

void
TransportElement :: giveSurfaceDofMapping(IntArray &answer, int iSurf)
{
    FEInterpolation *interp = this->giveInterpolation();
    if ( dynamic_cast< FEInterpolation3d * >(interp) ) {
        dynamic_cast< FEInterpolation3d * >(interp)->computeLocalSurfaceMapping(answer, iSurf);
    }
}

void
TransportElement :: computeSurfIpGlobalCoords(FloatArray &answer, const FloatArray &lcoord, int iSurf)
{
    FEInterpolation *interp = this->giveInterpolation();
    if ( dynamic_cast< FEInterpolation3d * >(interp) ) {
        dynamic_cast< FEInterpolation3d * >(interp)->surfaceLocal2global( answer, iSurf, lcoord, FEIElementGeometryWrapper(this) );
    }
}

int
TransportElement :: giveApproxOrder(int unknownIndx)
{
    return this->giveInterpolation()->giveInterpolationOrder();
}

void
TransportElement :: computeCapacitySubMatrix(FloatMatrix &answer, MatResponseMode rmode, int iri, TimeStep *tStep)
{
    FloatArray n;
    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    answer.clear();
    for ( GaussPoint *gp: *integrationRulesArray [ iri ] ) {
        this->computeNAt( n, gp->giveNaturalCoordinates() );
        // ask for capacity coefficient. In basic units [J/K/m3]
        double c = mat->giveCharacteristicValue(rmode, gp, tStep);
        double dV = this->computeVolumeAround(gp);
        answer.plusDyadSymmUpper(n, dV * c);
    }

    answer.symmetrized();
}

void
TransportElement :: computeConductivitySubMatrix(FloatMatrix &answer, int iri, MatResponseMode rmode, TimeStep *tStep)
{
    double dV;
    FloatMatrix b, d, db;

    answer.resize( this->giveNumberOfDofManagers(), this->giveNumberOfDofManagers() );
    answer.zero();
    for ( GaussPoint *gp: *integrationRulesArray [ iri ] ) {
        this->computeConstitutiveMatrixAt(d, rmode, gp, tStep);
        this->computeGradientMatrixAt(b, gp->giveNaturalCoordinates());
        dV = this->computeVolumeAround(gp);

        db.beProductOf(d, b);
        answer.plusProductSymmUpper(b, db, dV);
    }

    answer.symmetrized();
}

void
TransportElement :: computeInternalSourceRhsSubVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, int indx)
{
    // Computes numerically the generator Rhs vector of the receiver due to the generator
    //  at tStep.
    // // load is first transformed to local cs.
    // // load vector is then transformed to coordinate system in each node.
    // // (should be global coordinate system, but there may be defined
    // //  different coordinate system in each node)
    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    FloatArray val, n;
    answer.clear();

    // add internal source produced by material (if any)
    if ( mat->hasInternalSource() ) {
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            this->computeNAt( n, gp->giveNaturalCoordinates() );
            double dV = this->computeVolumeAround(gp);
            mat->computeInternalSourceVector(val, gp, tStep, mode);

            answer.add(val.at(indx) * dV, n);
        }
    }
}


void
TransportElement :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    if ( emode == HeatTransferEM || emode == Mass1TransferEM ) {
        this->computeInternalSourceRhsSubVectorAt(answer, tStep, mode, 1);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatArray subAnswer;

        for ( int i = 1; i <= 2; i++ ) {
            this->computeInternalSourceRhsSubVectorAt(subAnswer, tStep, mode, i);
            if ( subAnswer.isNotEmpty() ) {
                if ( answer.isEmpty() ) {
                    answer.resize( 2 * subAnswer.giveSize() );
                    answer.zero();
                }

                this->assembleLocalContribution(answer, subAnswer, 2, i);
            }
        }
    } else {
        OOFEM_ERROR("Unknown ElementMode");
    }
}


void
TransportElement :: computeIntSourceLHSMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );
    if ( mat->hasInternalSource() ) {
        answer.resize( this->computeNumberOfDofs(), this->computeNumberOfDofs() );
        answer.zero();

        if ( emode == HeatTransferEM || emode == Mass1TransferEM ) {
            this->computeIntSourceLHSSubMatrix(answer, IntSource, 0, tStep);
        } else if ( emode == HeatMass1TransferEM ) {
            FloatMatrix subAnswer;
            MatResponseMode rmode [ 2 ] = {
                IntSource_hh, IntSource_ww
            };
            //double coeff = 1.0; //this->giveMaterial()->give('d');

            for ( int i = 1; i <= 2; i++ ) {
                this->computeIntSourceLHSSubMatrix(subAnswer, rmode [ i - 1 ], 0, tStep);
                this->assembleLocalContribution(answer, subAnswer, 2, i, i);
            }
        } else {
            OOFEM_ERROR("Unknown ElementMode");
        }
    } else {
        answer.clear();
    }
}

void
TransportElement :: computeIntSourceLHSSubMatrix(FloatMatrix &answer, MatResponseMode rmode, int iri, TimeStep *tStep)
// computes LHS matrix due to material internal source (dHeat/dTemperature, dWaterSource/dw)
// IntSource - Heat transfer
// IntSource_hh - HeMo heat source
// IntSource_ww - HeMo water content source
// hw, wh is not used now
{
    FloatArray n;
    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    answer.clear();
    for ( GaussPoint *gp: *integrationRulesArray [ iri ] ) {
        this->computeNAt( n, gp->giveNaturalCoordinates() );
        // ask for coefficient from material
        double c = mat->giveCharacteristicValue(rmode, gp, tStep);
        double dV = this->computeVolumeAround(gp);
        answer.plusDyadSymmUpper(n, dV * c);
    }

    answer.symmetrized();
}


void
TransportElement :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                MatResponseMode rMode, GaussPoint *gp,
                                                TimeStep *tStep)
{
    static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicMatrix(answer, rMode, gp, tStep);
}

void
TransportElement :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FloatArray unknowns;
    this->computeVectorOf(VM_Total, tStep, unknowns);

    TransportMaterial *mat = static_cast< TransportMaterial* >( this->giveMaterial() );
    FloatArray flux, grad, field;
    FloatMatrix B, N;

    answer.clear();
    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        this->computeNmatrixAt(N, lcoords);
        this->computeBmatrixAt(B, lcoords);

        field.beProductOf(N, unknowns);
        grad.beProductOf(B, unknowns);

        mat->giveFluxVector(flux, gp, grad, field, tStep);

        double dV = this->computeVolumeAround(gp);
        answer.plusProduct(B, flux, -dV);

        ///@todo Can/should this part not be built into the flux itself? Probably not? / Mikael
        if ( mat->hasInternalSource() ) {
            // add internal source produced by material (if any)
            FloatArray val;
            mat->computeInternalSourceVector(val, gp, tStep, VM_Total);
            answer.plusProduct(N, val, -dV);
        }
    }

    ///@todo With sets, we can make this much nicer than to add it here, though it works for now. / Mikael
    // Neumann b.c.s
    FloatMatrix bc_tangent;
    this->computeBCMtrxAt(bc_tangent, tStep, VM_Total);
    if ( bc_tangent.isNotEmpty() ) {
        answer.plusProduct(bc_tangent, unknowns, 1.0);
    }
}


void
TransportElement :: computeInertiaForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FloatMatrix cap;
    FloatArray vel;
    this->computeVectorOf(VM_Velocity, tStep, vel);
    this->computeCapacityMatrix(cap, tStep);
    answer.beProductOf(cap, vel);
}


void
TransportElement :: computeLumpedCapacityVector(FloatArray &answer, TimeStep *tStep)
{
    FloatMatrix cap;
    this->computeCapacityMatrix(cap, tStep);

    int size = cap.giveNumberOfRows();
    answer.resize(size);
    answer.zero();
    for ( int i = 1; i <= size; i++ ) {
        for ( int j = 1; j <= size; j++ ) {
            answer.at(i) += cap.at(i, j);
        }
    }
}


void
TransportElement :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.clear();

    ///@todo We could technically support other types of volume loads. Even a bulk load is a type of transmission b.c.
#if 1
    if ( type != ExternalForcesVector ) {
        return;
    }
#else
    if ( !( load->giveType() == TransmissionBC && type == ExternalForcesVector ) &&
        !( load->giveType() == ConvectionBC && type == InternalForcesVector ) ) {
        return;
    }
#endif

    FloatArray gcoords, val, n, unknowns;
    FloatMatrix N;
    IntArray dofid;
    int unknownsPerNode = this->emode == HeatMass1TransferEM ? 2 : 1;

    this->giveElementDofIDMask(dofid);

    ///@todo Deal with coupled fields (I think they should be another class of problems completely).
    FEInterpolation *interp = this->giveInterpolation();
    //std :: unique_ptr< IntegrationRule > iRule( interp->giveIntegrationRule( load->giveApproxOrder() + 1 + interp->giveInterpolationOrder() ) );
    ///@todo FIXME backwards compatibility, the old tests used insufficient integration points for axisymm elements.
    std :: unique_ptr< IntegrationRule > iRule( interp->giveIntegrationRule( load->giveApproxOrder() ) );

    if ( load->giveType() == ConvectionBC || load->giveType() == RadiationBC ) {
        this->computeVectorOf(dofid, VM_Total, tStep, unknowns);
    }

    for ( GaussPoint *gp: *iRule ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        interp->evalN( n, lcoords, FEIElementGeometryWrapper(this) );
        N.beNMatrixOf(n, unknownsPerNode);

        double dV = this->computeVolumeAround(gp);

        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValueAt(val, tStep, lcoords, mode);
        } else {
            interp->local2global( gcoords, lcoords, FEIElementGeometryWrapper(this) );
            load->computeValueAt(val, tStep, gcoords, mode);
        }

        answer.plusProduct(N, val, dV);
    }
}


void
TransportElement :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.clear();

    if ( !( load->giveType() == TransmissionBC && type == ExternalForcesVector ) &&
        !( load->giveType() == ConvectionBC && type == InternalForcesVector ) &&
        !( load->giveType() == RadiationBC && type == InternalForcesVector ) ) {
        return;
    }

    FloatArray gcoords, val, n, unknowns, field;
    FloatMatrix N;
    IntArray dofid;
    int unknownsPerNode = this->emode == HeatMass1TransferEM ? 2 : 1;

    this->giveElementDofIDMask(dofid);

    ///@todo Deal with coupled fields (I think they should be another class of problems completely).
    FEInterpolation *interp = this->giveInterpolation();
    std :: unique_ptr< IntegrationRule > iRule( interp->giveBoundaryIntegrationRule(load->giveApproxOrder() + 1 + interp->giveInterpolationOrder(), boundary) );

    if ( load->giveType() == ConvectionBC || load->giveType() == RadiationBC ) {
        IntArray bNodes;
        interp->boundaryGiveNodes(bNodes, boundary);
        this->computeBoundaryVectorOf(bNodes, dofid, VM_Total, tStep, unknowns);
    }

    for ( GaussPoint *gp: *iRule ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        interp->boundaryEvalN( n, boundary, lcoords, FEIElementGeometryWrapper(this) );
        N.beNMatrixOf(n, unknownsPerNode);

        interp->boundaryLocal2Global( gcoords, boundary, lcoords, FEIElementGeometryWrapper(this) );
        double detJ = interp->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(this) );
        double dA = this->giveThicknessAt(gcoords) * gp->giveWeight() * detJ;

        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValueAt(val, tStep, lcoords, mode);
        } else {
            load->computeValueAt(val, tStep, gcoords, mode);
        }

        if ( load->giveType() == TransmissionBC ) {
            val.negated();
        } else if ( load->giveType() == ConvectionBC ) {
            field.beProductOf(N, unknowns);
            val.subtract(field);
            val.times( -1.0 * load->giveProperty('a', tStep) );
        } else if ( load->giveType() == RadiationBC ) {
            field.beProductOf(N, unknowns);
            val.subtract(field);
            val.times( -1.0 * getRadiativeHeatTranferCoef(load, tStep) );
        }

        answer.plusProduct(N, val, dA);
    }
}


void
TransportElement :: computeTangentFromBoundaryLoad(FloatMatrix &answer, BoundaryLoad *load, int boundary, MatResponseMode rmode, TimeStep *tStep)
{
    answer.clear();

    if ( load->giveType() != ConvectionBC && load->giveType() != RadiationBC ) {
        return;
    }

    FloatArray gcoords, n;
    FloatMatrix N;
    int unknownsPerNode = this->emode == HeatMass1TransferEM ? 2 : 1;

    ///@todo Deal with coupled fields (I think they should be another class of problems completely).
    FEInterpolation *interp = this->giveInterpolation();
    std :: unique_ptr< IntegrationRule > iRule( interp->giveBoundaryIntegrationRule(load->giveApproxOrder() + 1 + interp->giveInterpolationOrder(), boundary) );

    for ( auto &gp : *iRule ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        interp->boundaryEvalN( n, boundary, lcoords, FEIElementGeometryWrapper(this) );
        interp->boundaryLocal2Global( gcoords, boundary, lcoords, FEIElementGeometryWrapper(this) );
        double detJ = interp->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(this) );
        double dA = this->giveThicknessAt(gcoords) * gp->giveWeight() * detJ;
        N.beNMatrixOf(n, unknownsPerNode);

        if ( load->giveType() == ConvectionBC ){
            answer.plusProductSymmUpper(N, N, load->giveProperty('a', tStep) * dA);
        } else if ( load->giveType() == RadiationBC ){
            answer.plusProductSymmUpper(N, N, getRadiativeHeatTranferCoef(load, tStep) * dA);
        }
    }
    answer.symmetrized();
}


void
TransportElement :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.clear();

    if ( !( load->giveType() == TransmissionBC && type == ExternalForcesVector ) &&
        !( load->giveType() == ConvectionBC && type == InternalForcesVector ) &&
        !( load->giveType() == RadiationBC && type == InternalForcesVector ) ) {
        return;
    }

    FloatArray gcoords, val, n, unknowns, field;
    FloatMatrix N;
    IntArray dofid;
    int unknownsPerNode = 1;
    if ( this->emode == HeatMass1TransferEM ) {
        unknownsPerNode = 2;
    }

    this->giveElementDofIDMask(dofid);

    ///@todo Deal with coupled fields (I think they should be another class of problems completely).
    FEInterpolation *interp = this->giveInterpolation();
    std :: unique_ptr< IntegrationRule > iRule( interp->giveBoundaryEdgeIntegrationRule(load->giveApproxOrder() + 1 + interp->giveInterpolationOrder(), boundary) );

    if ( load->giveType() == ConvectionBC || load->giveType() == RadiationBC ) {
        IntArray bNodes;
        interp->boundaryGiveNodes(bNodes, boundary);
        this->computeBoundaryVectorOf(bNodes, dofid, VM_Total, tStep, unknowns);
    }

    for ( GaussPoint *gp: *iRule ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        interp->boundaryEdgeEvalN( n, boundary, lcoords, FEIElementGeometryWrapper(this) );
        N.beNMatrixOf(n, unknownsPerNode);

        double detJ = interp->boundaryEdgeGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(this) );
        interp->boundaryEdgeLocal2Global( gcoords, boundary, lcoords, FEIElementGeometryWrapper(this) );
        double dL = this->giveThicknessAt(gcoords) * gp->giveWeight() * detJ;

        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValueAt(val, tStep, lcoords, mode);
        } else {
            load->computeValueAt(val, tStep, gcoords, mode);
        }

        if ( load->giveType() == TransmissionBC ) {
            val.negated();
        } else if ( load->giveType() == ConvectionBC ) {
            field.beProductOf(N, unknowns);
            val.subtract(field);
            val.times( -1.0 * load->giveProperty('a', tStep) );
        } else if ( load->giveType() == RadiationBC  ) {
            //actual Temperature in C in field
            field.beProductOf(N, unknowns);
            val.subtract(field);
            val.times( -1.0 * getRadiativeHeatTranferCoef(load, tStep) );
        }
        answer.plusProduct(N, val, dL);
    }
}


void
TransportElement :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    this->computeBCVectorAt(answer, tStep, mode);
}

void
TransportElement :: computeBCVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    answer.resize( this->computeNumberOfDofs() );
    answer.zero();

    if ( emode == HeatTransferEM || emode == Mass1TransferEM ) {
        this->computeBCSubVectorAt(answer, tStep, mode, 1);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatArray subAnswer;

        for ( int i = 1; i <= 2; i++ ) {
            this->computeBCSubVectorAt(subAnswer, tStep, mode, i);
            this->assembleLocalContribution(answer, subAnswer, 2, i);
        }
    } else {
        OOFEM_ERROR("Unknown ElementMode");
    }
}


void
TransportElement :: computeBCMtrxAt(FloatMatrix &answer, TimeStep *tStep, ValueModeType mode)
{
    int ndofs = this->computeNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    if ( emode == HeatTransferEM || emode == Mass1TransferEM ) {
        this->computeBCSubMtrxAt(answer, tStep, mode, 1);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatMatrix subAnswer;

        for ( int i = 1; i <= 2; i++ ) {
            this->computeBCSubMtrxAt(subAnswer, tStep, mode, i);
            if ( subAnswer.isNotEmpty() ) {
                this->assembleLocalContribution(answer, subAnswer, 2, i, i);
            }
        }
    } else {
        OOFEM_ERROR("Unknown ElementMode");
    }
}


void
TransportElement :: computeBCSubVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, int indx)
{
    FloatArray vec;

    answer.resize( this->giveNumberOfDofManagers() );
    answer.zero();

    // loop over boundary load array
    int nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        Load *load = domain->giveLoad(n);
        if ( load->isDofExcluded(indx) || !load->isImposed(tStep) ) {
            continue;
        }

        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, load, id, tStep, mode, indx);
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceBCSubVectorAt(vec, load, id, tStep, mode, indx);
        } else {
            OOFEM_ERROR("unsupported bc type encountered");
        }

        answer.add(vec);
    }

    for ( int i = 1; i <= bodyLoadArray.giveSize(); ++i ) {
        Load *load = domain->giveLoad(bodyLoadArray.at(i));
        if ( !load->isImposed(tStep) ) {
            continue;
        }
        this->computeBodyBCSubVectorAt(answer, load, tStep, mode, indx);
    }
}


void
TransportElement :: computeBodyBCSubVectorAt(FloatArray &answer, Load *load,
                                             TimeStep *tStep, ValueModeType mode, int indx)
{
    FloatArray val, globalIPcoords, n;
    answer.resize( this->giveNumberOfDofManagers() );

    std :: unique_ptr< IntegrationRule > iRule( this->giveInterpolation()->giveIntegrationRule(load->giveApproxOrder()) );
    for ( GaussPoint *gp : *iRule ) {
        double dV = this->computeVolumeAround(gp);
        this->computeNAt( n, gp->giveNaturalCoordinates() );
        this->computeGlobalCoordinates( globalIPcoords, gp->giveNaturalCoordinates() );
        load->computeValueAt(val, tStep, globalIPcoords, mode);
        answer.add(val.at(indx) * dV, n);
    }
}


void
TransportElement :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge,
                                             TimeStep *tStep, ValueModeType mode, int indx)
{
    answer.resize( this->giveNumberOfDofManagers() );
    answer.zero();

    if ( ( load->giveType() == TransmissionBC ) || ( load->giveType() == ConvectionBC ) || load->giveType() == RadiationBC ) {
        BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);

        int approxOrder = edgeLoad->giveApproxOrder() + this->giveApproxOrder(indx);
        int numberOfEdgeIPs = ( int ) ceil( ( approxOrder + 1. ) / 2. );
        GaussIntegrationRule iRule(1, this, 1, 1);
        iRule.SetUpPointsOnLine(numberOfEdgeIPs, _Unknown);
        FloatArray reducedAnswer, val, n;
        IntArray mask;
        double dV, coeff = 1.0;

        for ( GaussPoint *gp: iRule ) {
            if( edgeLoad->propertyMultExpr.isDefined()) {//dependence on state variable
                ///@todo Deal with coupled fields
                if ( emode!=HeatTransferEM && emode!=Mass1TransferEM ){
                    OOFEM_ERROR("Not implemented for >=2 coupled fields");
                }
                FloatArray unknowns;
                IntArray dofid;
                this->computeEgdeNAt( n, iEdge, gp->giveNaturalCoordinates() );
                this->giveElementDofIDMask(dofid);
                this->giveEdgeDofMapping(mask, iEdge);
                this->computeBoundaryVectorOf(mask, dofid, VM_Total, tStep, unknowns);
                double value = n.dotProduct(unknowns);//unknown in IP
                edgeLoad->setVariableState('x', value);
            }
            if ( load->giveType() == TransmissionBC ) {
                coeff = -1.0;
            } else if ( load->giveType() == ConvectionBC ) {
                coeff = edgeLoad->giveProperty('a', tStep);
            } else if ( load->giveType() == RadiationBC ){
                coeff = getRadiativeHeatTranferCoef(edgeLoad, tStep);
            }

            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            this->computeEgdeNAt(n, iEdge, lcoords);
            dV = this->computeEdgeVolumeAround(gp, iEdge);

            if ( edgeLoad->giveFormulationType() == Load :: FT_Entity ) {
                edgeLoad->computeValueAt(val, tStep, lcoords, mode);
            } else {
                FloatArray globalIPcoords;
                this->computeEdgeIpGlobalCoords(globalIPcoords, lcoords, iEdge);
                edgeLoad->computeValueAt(val, tStep, globalIPcoords, mode);
            }

            reducedAnswer.add(val.at(indx) * coeff * dV, n);
        }

        this->giveEdgeDofMapping(mask, iEdge);
        answer.assemble(reducedAnswer, mask);
    } else {
        OOFEM_ERROR("unsupported bc type encountered");
    }
}


void
TransportElement :: computeSurfaceBCSubVectorAt(FloatArray &answer, Load *load,
                                                int iSurf, TimeStep *tStep, ValueModeType mode, int indx)
{
    if ( !this->testElementExtension(Element_SurfaceLoadSupport) ) {
        OOFEM_ERROR("no surface load support");
    }

    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( surfLoad ) {
        FloatArray reducedAnswer, val, globalIPcoords, n;
        IntArray mask;

        answer.resize( this->giveNumberOfDofManagers() );
        answer.zero();

        double coeff=0.;
        int approxOrder = surfLoad->giveApproxOrder() + this->giveApproxOrder(indx);

        std :: unique_ptr< IntegrationRule > iRule( this->GetSurfaceIntegrationRule(approxOrder) );
        for ( GaussPoint *gp: *iRule ) {
            if( surfLoad->propertyMultExpr.isDefined()) {//dependence on state variable
                    if ( emode!=HeatTransferEM && emode!=Mass1TransferEM ){
                        ///@todo Deal with coupled fields
                        OOFEM_ERROR("Not implemented for >=2 coupled fields");
                    }
                    FloatArray unknowns;
                    IntArray dofid;
                    this->computeSurfaceNAt( n, iSurf, gp->giveNaturalCoordinates() );
                    this->giveElementDofIDMask(dofid);
                    this->giveSurfaceDofMapping(mask, iSurf);
                    this->computeBoundaryVectorOf(mask, dofid, VM_Total, tStep, unknowns);
                    double value = n.dotProduct(unknowns);//unknown in IP
                    surfLoad->setVariableState('x', value);
                }
            if ( load->giveType() == TransmissionBC ) {
                coeff = -1.0;
            } else if ( load->giveType() == ConvectionBC ) {
                coeff = surfLoad->giveProperty('a', tStep);
            } else if ( load->giveType() == RadiationBC ) {
                coeff = getRadiativeHeatTranferCoef(surfLoad, tStep);
            } else {
                OOFEM_ERROR("Unknown load type");
            }

            this->computeSurfaceNAt( n, iSurf, gp->giveNaturalCoordinates() );
            double dV = this->computeSurfaceVolumeAround(gp, iSurf);

            if ( surfLoad->giveFormulationType() == Load :: FT_Entity ) {
                surfLoad->computeValueAt(val, tStep, gp->giveNaturalCoordinates(), mode);
            } else {
                this->computeSurfIpGlobalCoords(globalIPcoords, gp->giveNaturalCoordinates(), iSurf);
                surfLoad->computeValueAt(val, tStep, globalIPcoords, mode);
            }

            reducedAnswer.add(val.at(indx) * coeff * dV, n);
        }

        this->giveSurfaceDofMapping(mask, iSurf);
        answer.assemble(reducedAnswer, mask);
    } else {
        OOFEM_ERROR("unsupported bc type encountered");
    }
}


void
TransportElement :: computeBCSubMtrxAt(FloatMatrix &answer, TimeStep *tStep, ValueModeType mode, int indx)
{
    int defined = 0;

    answer.resize( this->giveNumberOfDofManagers(), this->giveNumberOfDofManagers() );
    answer.zero();

    // loop over boundary load array
    int nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        int k = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        GeneralBoundaryCondition *load = domain->giveLoad(k);
        if ( load->giveType() == ConvectionBC || load->giveType() == RadiationBC ) {
            bcGeomType ltype = load->giveBCGeoType();
            if ( ltype == EdgeLoadBGT ) {
                BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
                if ( edgeLoad->isDofExcluded(indx) || !edgeLoad->isImposed(tStep) ) {
                    continue;
                }

                defined = 1;
                int approxOrder = 2 * this->giveApproxOrder(indx);
                int numberOfEdgeIPs = ( int ) ceil( ( approxOrder + 1. ) / 2. );
                GaussIntegrationRule iRule(1, this, 1, 1);
                iRule.SetUpPointsOnLine(numberOfEdgeIPs, _Unknown);
                FloatArray n;
                IntArray mask;
                FloatMatrix subAnswer;

                for ( GaussPoint *gp: iRule ) {
                    this->computeEgdeNAt( n, id, gp->giveNaturalCoordinates() );
                    double dV = this->computeEdgeVolumeAround(gp, id);
                    if( edgeLoad->propertyMultExpr.isDefined()) {//dependence on state variable
                        if ( emode!=HeatTransferEM && emode!=Mass1TransferEM ){
                            ///@todo Deal with coupled fields
                            OOFEM_ERROR("Not implemented for >=2 coupled fields");
                        }
                        FloatArray unknowns;
                        IntArray dofid;
                        this->giveElementDofIDMask(dofid);
                        this->giveEdgeDofMapping(mask, id);
                        this->computeBoundaryVectorOf(mask, dofid, VM_Total, tStep, unknowns);
                        double value = n.dotProduct(unknowns);//unknown in IP
                        edgeLoad->setVariableState('x', value);
                    }
                    if ( load->giveType() == ConvectionBC ) {
                        subAnswer.plusDyadSymmUpper( n, dV * edgeLoad->giveProperty('a', tStep) );
                    } else if ( load->giveType() == RadiationBC ) {
                        subAnswer.plusDyadSymmUpper( n, dV * getRadiativeHeatTranferCoef(edgeLoad, tStep) );
                    }
                }

                subAnswer.symmetrized();
                this->giveEdgeDofMapping(mask, id);
                answer.assemble(subAnswer, mask);
            } else if ( ltype == SurfaceLoadBGT ) {
                FloatArray n;
                IntArray mask;
                FloatMatrix subAnswer;

                BoundaryLoad *surfLoad = static_cast< BoundaryLoad * >(load);
                if ( surfLoad->isDofExcluded(indx) || !surfLoad->isImposed(tStep) ) {
                    continue;
                }

                defined = 1;
                int approxOrder = 2 * this->giveApproxOrder(indx);
                std :: unique_ptr< IntegrationRule > iRule( this->GetSurfaceIntegrationRule(approxOrder) );

                for ( GaussPoint *gp: *iRule ) {
                    this->computeSurfaceNAt( n, id, gp->giveNaturalCoordinates() );
                    double dV = this->computeSurfaceVolumeAround(gp, id);
                    if( surfLoad->propertyMultExpr.isDefined()) {//dependence on state variable
                        if ( emode!=HeatTransferEM && emode!=Mass1TransferEM ){
                            ///@todo Deal with coupled fields
                            OOFEM_ERROR("Not implemented for >=2 coupled fields");
                        }
                        FloatArray unknowns;
                        IntArray dofid;
                        this->giveElementDofIDMask(dofid);
                        this->giveSurfaceDofMapping(mask, id);
                        this->computeBoundaryVectorOf(mask, dofid, VM_Total, tStep, unknowns);
                        double value = n.dotProduct(unknowns);//unknown in IP
                        surfLoad->setVariableState('x', value);
                    }
                    if ( load->giveType() == ConvectionBC ) {
                        subAnswer.plusDyadSymmUpper( n, dV * surfLoad->giveProperty('a', tStep) );
                    } else if ( load->giveType() == RadiationBC ) {
                        subAnswer.plusDyadSymmUpper( n, dV * getRadiativeHeatTranferCoef(surfLoad, tStep) );
                    }
                }

                subAnswer.symmetrized();
                this->giveSurfaceDofMapping(mask, id);
                answer.assemble(subAnswer, mask);
            } else {
                OOFEM_ERROR("unsupported bc type encountered");
            }
        }
    }

    // end loop over applied bc

    if ( !defined ) {
        answer.clear();
    }
}


void
TransportElement :: assembleLocalContribution(FloatMatrix &answer, FloatMatrix &src,
                                              int ndofs, int rdof, int cdof)
{
    int nnodes = this->giveNumberOfDofManagers();

    for ( int i = 1; i <= nnodes; i++ ) {
        int ti = ( i - 1 ) * ndofs + rdof;
        for ( int j = 1; j <= nnodes; j++ ) {
            int tj = ( j - 1 ) * ndofs + cdof;
            answer.at(ti, tj) += src.at(i, j);
        }
    }
}


void
TransportElement :: assembleLocalContribution(FloatArray &answer, FloatArray &src,
                                              int ndofs, int rdof)
{
    int nnodes = this->giveNumberOfDofManagers();

    for ( int i = 1; i <= nnodes; i++ ) {
        int ti = ( i - 1 ) * ndofs + rdof;
        answer.at(ti) += src.at(i);
    }
}

double
TransportElement :: getRadiativeHeatTranferCoef(BoundaryLoad *bLoad, TimeStep *tStep)
{
    double answer = 0;
    ///@todo Why aren't this code using the standard approach of calling computeComponentArrayAt(...) to get the time function scaling and all?
    const FloatArray &components = bLoad->giveComponentArray();
    
    answer = components.at(1);//T_infty
    answer += 273.15;
    answer = answer*answer*answer;
    answer *= 4 * bLoad->giveProperty('e', tStep) * stefanBoltzmann;
    return answer;
}


void
TransportElement :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray u;
    FloatMatrix n;
    this->computeNmatrixAt(n, lcoords);
    this->computeVectorOf(mode, tStep, u);
    answer.beProductOf(n, u);
}

void
TransportElement :: computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray r, br;
    FloatMatrix b, d;
    IntArray dofid;

    this->giveElementDofIDMask(dofid);
    this->computeVectorOf(dofid, VM_Total, tStep, r);
    this->computeGradientMatrixAt(b, gp->giveNaturalCoordinates());

    if ( emode == HeatTransferEM ||  emode == Mass1TransferEM ) {
        this->computeConstitutiveMatrixAt(d, Conductivity_hh, gp, tStep);
        br.beProductOf(b, r);
        answer.beProductOf(d, br);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatArray r_h, r_w;
        FloatMatrix b_tot, d_hh, d_hw, d_wh, d_ww;
        r_h.resize(r.giveSize() / 2);
        r_h.zero();
        r_w.resize(r.giveSize() / 2);
        r_w.zero();

        this->computeConstitutiveMatrixAt(d_hh, Conductivity_hh, gp, tStep);
        this->computeConstitutiveMatrixAt(d_hw, Conductivity_hw, gp, tStep);
        this->computeConstitutiveMatrixAt(d_wh, Conductivity_wh, gp, tStep);
        this->computeConstitutiveMatrixAt(d_ww, Conductivity_ww, gp, tStep);

        b_tot.resize( 2 * b.giveNumberOfRows(), 2 * b.giveNumberOfColumns() );
        b_tot.setSubMatrix(b, 1, 1);
        b_tot.setSubMatrix(b, b.giveNumberOfRows() + 1, b.giveNumberOfColumns() + 1);
        br.beProductOf(b_tot, r);

        d.resize( 2 * d_hh.giveNumberOfRows(), 2 * d_hh.giveNumberOfColumns() );
        d.setSubMatrix(d_hh, 1, 1);
        d.setSubMatrix(d_hw, 1, d_hh.giveNumberOfColumns() + 1);
        d.setSubMatrix(d_wh, d_hh.giveNumberOfRows() + 1, 1);
        d.setSubMatrix(d_ww, d_hh.giveNumberOfRows() + 1, d_wh.giveNumberOfColumns() + 1);

        answer.beProductOf(d, br);
    } else {
        OOFEM_ERROR("Unknown element mode");
    }

    answer.negated();

    if ( !isActivated(tStep) ) {
        answer.zero();
    }
}


void
TransportElement :: updateInternalState(TimeStep *tStep)
// Updates the receiver at the end of a solution step
{
    FloatArray stateVector, r;
#if 0
    FloatArray gradient, flux;
    FloatMatrix B;
#endif
    FloatMatrix n;
    IntArray dofid;
    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    this->giveElementDofIDMask(dofid);
    this->computeVectorOf(dofid, VM_Total, tStep, r);
    // force updating ip values
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {

            ///@todo Why is the state vector the unknown solution at the gauss point? / Mikael
            this->computeNmatrixAt( n, gp->giveNaturalCoordinates() );
            stateVector.beProductOf(n, r);
            mat->updateInternalState(stateVector, gp, tStep);

            ///@todo We need to sort out multiple materials for coupled (heat+mass) problems
#if 0
            this->computeGradientMatrixAt( B, gp->giveNaturalCoordinates() );
            gradient.beProductOf(B, r);
            mat->giveFluxVector(flux, gp, gradient, tStep);
#endif
        }
    }
}

int
TransportElement :: EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                          FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                          TimeStep *tStep)
{
    int indx;
    FloatArray elemvector, lc;
    FloatMatrix n;
    IntArray elemdofs;
    // determine element dof ids
    this->giveElementDofIDMask(elemdofs);
    // first evaluate element unknown vector
    this->computeVectorOf(pf, elemdofs, mode, tStep, elemvector);
    // determine corresponding local coordinates
    if ( this->computeLocalCoordinates(lc, coords) ) {
        // compute interpolation matrix
        this->computeNmatrixAt(n, lc);
        // compute answer
        answer.resize( dofId.giveSize() );
        answer.zero();
        for ( int i = 1; i <= dofId.giveSize(); i++ ) {
            if ( ( indx = elemdofs.findFirstIndexOf( dofId.at(i) ) ) ) {
                double sum = 0.0;
                for ( int j = 1; j <= elemvector.giveSize(); j++ ) {
                    sum += n.at(indx, j) * elemvector.at(j);
                }

                answer.at(i) = sum;
            } else {
                //_error("EIPrimaryFieldI_evaluateFieldVectorAt: unknown dof id encountered");
                answer.at(i) = 0.0;
            }
        }

        return 0; // ok
    } else {
        OOFEM_ERROR("target point not in receiver volume");
        return 1; // failed
    }
}


TransportCrossSection *
TransportElement :: giveTransportCrossSection()
{
    return static_cast< TransportCrossSection* >( this->giveCrossSection() );
}


Material *
TransportElement :: giveMaterial()
{
    // This method is overloaded in order to have some backwards compatibility
    if ( material ) {
        return domain->giveMaterial(this->material);
    }
    return static_cast< TransportCrossSection* >( this->giveCrossSection() )->giveMaterial();
}


#ifdef __OOFEG
int
TransportElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                            int node, TimeStep *tStep)
{
    Node *n = this->giveNode(node);
    if ( type == IST_Temperature ) {
        auto dofindx = n->findDofWithDofId(T_f);
        if ( dofindx != n->end() ) {
            answer.resize(1);
            answer.at(1) = (*dofindx)->giveUnknown(VM_Total, tStep);
            return 1;
        } else {
            return 0;
        }
    } else if ( type == IST_MassConcentration_1 ) {
        auto dofindx = n->findDofWithDofId(C_1);
        if ( dofindx != n->end() ) {
            answer.resize(1);
            answer.at(1) = (*dofindx)->giveUnknown(VM_Total, tStep);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, tStep);
    }
}

#endif
} // end namespace oofem
