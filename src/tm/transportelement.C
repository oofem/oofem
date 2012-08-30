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

#include "transportelement.h"
#include "domain.h"
#include "transportmaterial.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "verbose.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
TransportElement :: TransportElement(int n, Domain *aDomain, ElementMode em) :
    Element(n, aDomain)
{
    emode = em;
}


TransportElement :: ~TransportElement()
{ }


void
TransportElement :: giveElementDofIDMask(EquationID, IntArray &answer) const
{
    if ( emode == HeatTransferEM ) {
        answer.setValues(1, T_f);
    } else if ( emode == HeatMass1TransferEM ) {
        answer.setValues(2, T_f, C_1);
    } else {
        _error("Unknown ElementMode");
    }
}


void
TransportElement :: giveCharacteristicMatrix(FloatMatrix &answer,
                                             CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == ConductivityMatrix ) {
        this->computeConductivityMatrix(answer, Conductivity, tStep);
    } else if ( mtrx == CapacityMatrix ) {
        this->computeCapacityMatrix(answer, tStep);
    } else if ( mtrx == LHSBCMatrix ) {
        this->computeBCMtrxAt(answer, tStep, VM_Total);
    } else if ( mtrx == IntSourceLHSMatrix ) {
        this->computeIntSourceLHSMatrix(answer, tStep);
    } else {
        _error2( "giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}


void
TransportElement :: giveDofManDofIDMask(int inode, EquationID eid, IntArray &answer) const
{
    if ( eid == EID_ConservationEquation ) {
        if ( emode == HeatTransferEM ) {
            answer.setValues(1, T_f);
        } else if ( emode == HeatMass1TransferEM ) {
            answer.setValues(2, T_f, C_1);
        } else {
            _error("Unknown ElementMode");
        }
    } else {
        answer.resize(0);
    }
}


void
TransportElement ::  giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                              TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == ElementBCTransportVector ) {
        this->computeBCVectorAt(answer, tStep, mode);
    } else if ( mtrx == ElementInternalSourceVector ) {
        this->computeInternalSourceRhsVectorAt(answer, tStep, mode);
    } else {
        _error2( "giveCharacteristicVector: Unknown Type of characteristic mtrx (%s)",
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
    if ( !dynamic_cast<TransportMaterial*>(giveMaterial()) ) {
        _warning("checkConsistency : material without support for transport problems");
        result = 0;
    }

    /*
     * if (!this->giveCrossSection()->testCrossSectionExtension(CS_TransportCapability)) {
     * this->warning("checkConsistency : cross-section without support for transport problems", 1);
     * result =0;
     * }
     */
    return result;
}

void
TransportElement :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    Element :: printOutputAt(file, stepN);
}

void
TransportElement :: computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize( computeNumberOfDofs(EID_ConservationEquation), computeNumberOfDofs(EID_ConservationEquation) );
    answer.zero();

    if ( emode == HeatTransferEM ) {
        this->computeCapacitySubMatrix(answer, Capacity, 0, tStep);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatMatrix subAnswer;
        MatResponseMode rmode [ 2 ] = {
            Capacity_hh, Capacity_ww
        };
        //double coeff = 1.0; //this->giveMaterial()->give('d');

        for (int i = 1; i <= 2; i++ ) {
            this->computeCapacitySubMatrix(subAnswer, rmode [ i - 1 ], 0, tStep);
            this->assembleLocalContribution(answer, subAnswer, 2, i, i);
        }
    } else {
        _error("Unknown ElementMode");
    }
}

void
TransportElement :: computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    answer.resize( computeNumberOfDofs(EID_ConservationEquation), computeNumberOfDofs(EID_ConservationEquation) );
    answer.zero();
    if ( emode == HeatTransferEM ) {
        this->computeConductivitySubMatrix(answer, 2, 0, Conductivity_hh, tStep);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatMatrix subAnswer;
        MatResponseMode rmode [ 2 ] [ 2 ] = { { Conductivity_hh, Conductivity_hw }, { Conductivity_wh, Conductivity_ww } };

        for (int i = 1; i <= 2; i++ ) {
            for (int j = 1; j <= 2; j++ ) {
                this->computeConductivitySubMatrix(subAnswer, 2, 0, rmode [ i - 1 ] [ j - 1 ], tStep);
                this->assembleLocalContribution(answer, subAnswer, 2, i, j);
            }
        }
    } else {
        _error("Unknown ElementMode");
    }
}


void 
TransportElement :: computeNAt(FloatArray &answer, const FloatArray &lcoord)
{
    this->giveInterpolation()->evalN(answer, lcoord, FEIElementGeometryWrapper(this));
}


void 
TransportElement :: computeNmatrixAt(FloatMatrix &answer, const FloatArray &lcoords)
{
    FloatArray n;
    this->computeNAt(n, lcoords);
    int q = n.giveSize();    
    if ( this->emode == HeatTransferEM ) {
        answer.resize(1, q);
        for (int i = 1; i <= q; i++ ) {
            answer.at(1, i) = n.at(i);
        }
    } else {
        answer.resize(2, 2*q);
        for (int i = 0; i < 2; i++ ) {
            for (int j = 0; j < q; j++ ) {
                answer(i, j * 2 + i) = n(j);
            }
        }
    }
}


void
TransportElement :: computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dnx;
    this->giveInterpolation()->evaldNdx(dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    ///@todo We should change the transposition in evaldNdx;
    answer.beTranspositionOf(dnx);
}


void
TransportElement :: computeEgdeNAt(FloatArray &answer, const FloatArray &lcoords)
{
    FEInterpolation *interp = this->giveInterpolation();
    if (dynamic_cast<FEInterpolation2d*>(interp)) {
        dynamic_cast<FEInterpolation2d*>(interp)->edgeEvalN(answer, lcoords, FEIElementGeometryWrapper(this));
    } else if (dynamic_cast<FEInterpolation3d*>(interp)) {
        dynamic_cast<FEInterpolation3d*>(interp)->edgeEvalN(answer, lcoords, FEIElementGeometryWrapper(this));
    }    
}


void
TransportElement :: giveEdgeDofMapping(IntArray &answer, int iEdge)
{
    FEInterpolation *interp = this->giveInterpolation();
    if (dynamic_cast<FEInterpolation2d*>(interp)) {
        dynamic_cast<FEInterpolation2d*>(interp)->computeLocalEdgeMapping(answer, iEdge);
    } else if (dynamic_cast<FEInterpolation3d*>(interp)) {
        dynamic_cast<FEInterpolation3d*>(interp)->computeLocalEdgeMapping(answer, iEdge);
    }
}


void
TransportElement :: computeEdgeIpGlobalCoords(FloatArray &answer, const FloatArray &lcoords, int iEdge)
{
    FEInterpolation *interp = this->giveInterpolation();
    if (dynamic_cast<FEInterpolation2d*>(interp)) {
        dynamic_cast<FEInterpolation3d*>(interp)->edgeLocal2global(answer, iEdge, lcoords, FEIElementGeometryWrapper(this));
    } else if (dynamic_cast<FEInterpolation3d*>(interp)) {
        dynamic_cast<FEInterpolation3d*>(interp)->edgeLocal2global(answer, iEdge, lcoords, FEIElementGeometryWrapper(this));
    }
}

void
TransportElement :: computeSurfaceNAt(FloatArray &answer, const FloatArray &lcoord)
{
    FEInterpolation *interp = this->giveInterpolation();
    if (dynamic_cast<FEInterpolation3d*>(interp)) {
        dynamic_cast<FEInterpolation3d*>(interp)->surfaceEvalN(answer, lcoord, FEIElementGeometryWrapper(this));
    }
}

void
TransportElement :: giveSurfaceDofMapping(IntArray &answer, int iSurf)
{
    FEInterpolation *interp = this->giveInterpolation();
    if (dynamic_cast<FEInterpolation3d*>(interp)) {
        dynamic_cast<FEInterpolation3d*>(interp)->computeLocalSurfaceMapping(answer, iSurf);
    }
}

void
TransportElement :: computeSurfIpGlobalCoords(FloatArray &answer, const FloatArray &lcoord, int iSurf)
{
    FEInterpolation *interp = this->giveInterpolation();
    if (dynamic_cast<FEInterpolation3d*>(interp)) {
        dynamic_cast<FEInterpolation3d*>(interp)->edgeLocal2global(answer, iSurf, lcoord, FEIElementGeometryWrapper(this));
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
    double dV, c;
    FloatArray n;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ iri ];

    answer.beEmptyMtrx();
    for (int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeNAt( n, *gp->giveCoordinates() );
        // ask for capacity coefficient. In basic units [J/K/m3]
        c = ( ( TransportMaterial * ) this->giveMaterial() )->giveCharacteristicValue(rmode, gp, tStep);
        dV = this->computeVolumeAround(gp);
        answer.plusDyadSymmUpper(n, n, dV * c);
    }

    answer.symmetrized();
}

void
TransportElement :: computeConductivitySubMatrix(FloatMatrix &answer, int nsd, int iri, MatResponseMode rmode, TimeStep *tStep)
{
    double dV;
    FloatMatrix b, d, db;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ iri ];

    answer.resize( this->giveNumberOfDofManagers(), this->giveNumberOfDofManagers() );
    answer.zero();
    for (int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeConstitutiveMatrixAt(d, rmode, gp, tStep);
        this->computeGradientMatrixAt(b, gp);
        dV = this->computeVolumeAround(gp);

        db.beProductOf(d, b);
        answer.plusProductSymmUpper(b, db, dV);
        //answer.plusProductUnsym(b,db,dV) ;
    }

    answer.symmetrized();
}

void
TransportElement :: computeInternalSourceRhsSubVectorAt(FloatArray &answer, TimeStep *atTime, ValueModeType mode, int indx)
{
    // Computes numerically the generator Rhs vector of the receiver due to the generator
    //  at stepN.
    // // load is first transformed to local cs.
    // // load vector is then transformed to coordinate system in each node.
    // // (should be global coordinate system, but there may be defined
    // //  different coordinate system in each node)
    int i, igp, k, nLoads;
    double dV;
    bcGeomType ltype;
    Load *load;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    TransportMaterial *mat = ( ( TransportMaterial * ) this->giveMaterial() );
    GaussPoint *gp;


    FloatArray val, helpLoadVector, globalIPcoords, n;
    answer.resize(0);

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        k     = bodyLoadArray.at(i);
        load  = ( Load * ) domain->giveLoad(k);
        ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT ) {
            for ( igp = 0; igp < iRule->getNumberOfIntegrationPoints(); igp++ ) {
                gp  = iRule->getIntegrationPoint(igp);
                this->computeNAt( n, *gp->giveCoordinates() );
                dV  = this->computeVolumeAround(gp);
                this->computeGlobalCoordinates( globalIPcoords, * gp->giveCoordinates() );
                load->computeValueAt(val, atTime, globalIPcoords, mode);

                helpLoadVector.add(val.at(indx) * dV, n);
            }

            answer.add(helpLoadVector);
        }
    }

    // add internal source produced by material (if any)
    if ( mat->hasInternalSource() ) {
        for ( igp = 0; igp < iRule->getNumberOfIntegrationPoints(); igp++ ) {
            gp  = iRule->getIntegrationPoint(igp);
            this->computeNAt( n, *gp->giveCoordinates() );
            dV  = this->computeVolumeAround(gp);
            mat->computeInternalSourceVector(val, gp, atTime, mode);
            
            helpLoadVector.add(val.at(indx) * dV, n);
        }

        answer.add(helpLoadVector);
    }
}


void
TransportElement :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *atTime, ValueModeType mode)
{
    if ( emode == HeatTransferEM ) {
        this->computeInternalSourceRhsSubVectorAt(answer, atTime, mode, 1);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatArray subAnswer;

        for (int i = 1; i <= 2; i++ ) {
            this->computeInternalSourceRhsSubVectorAt(subAnswer, atTime, mode, i);
            if ( subAnswer.isNotEmpty() ) {
                if ( answer.isEmpty() ) {
                    answer.resize(2*subAnswer.giveSize());
                    answer.zero();
                }

                this->assembleLocalContribution(answer, subAnswer, 2, i);
            }
        }
    } else {
        _error("Unknown ElementMode");
    }
}


void
TransportElement :: computeIntSourceLHSMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    TransportMaterial *mat = ( TransportMaterial * ) this->giveMaterial();
    if ( mat->hasInternalSource() ) {
        answer.resize( computeNumberOfDofs(EID_ConservationEquation), computeNumberOfDofs(EID_ConservationEquation) );
        answer.zero();

        if ( emode == HeatTransferEM ) {
            this->computeIntSourceLHSSubMatrix(answer, IntSource, 0, tStep);
        } else if ( emode == HeatMass1TransferEM ) {
            FloatMatrix subAnswer;
            int i;
            MatResponseMode rmode [ 2 ] = {
                IntSource_hh, IntSource_ww
            };
            //double coeff = 1.0; //this->giveMaterial()->give('d');

            for ( i = 1; i <= 2; i++ ) {
                this->computeIntSourceLHSSubMatrix(subAnswer, rmode [ i - 1 ], 0, tStep);
                this->assembleLocalContribution(answer, subAnswer, 2, i, i);
            }
        } else {
            _error("Unknown ElementMode");
        }
    } else {
        answer.beEmptyMtrx();
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
    double dV, c;
    FloatArray n;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ iri ];

    answer.beEmptyMtrx();
    for (int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeNAt( n, *gp->giveCoordinates() );
        // ask for coefficient from material
        c = ( ( TransportMaterial * ) this->giveMaterial() )->giveCharacteristicValue(rmode, gp, tStep);
        dV = this->computeVolumeAround(gp);
        answer.plusDyadSymmUpper(n, n, dV * c);
    }

    answer.symmetrized();
}


void
TransportElement :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                MatResponseMode rMode, GaussPoint *gp,
                                                TimeStep *tStep)
{
    ( ( TransportMaterial * ) this->giveMaterial() )->giveCharacteristicMatrix(answer, FullForm, rMode, gp, tStep);
}


void
TransportElement :: computeBCVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    answer.resize( computeNumberOfDofs(EID_ConservationEquation) );
    answer.zero();

    if ( emode == HeatTransferEM ) {
        this->computeBCSubVectorAt(answer, tStep, mode, 1);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatArray subAnswer;

        for (int i = 1; i <= 2; i++ ) {
            this->computeBCSubVectorAt(subAnswer, tStep, mode, i);
            this->assembleLocalContribution(answer, subAnswer, 2, i);
        }
    } else {
        _error("Unknown ElementMode");
    }
}

void
TransportElement :: computeBCMtrxAt(FloatMatrix &answer, TimeStep *tStep, ValueModeType mode)
{
    int ndofs = computeNumberOfDofs(EID_ConservationEquation);
    answer.resize(ndofs, ndofs);
    answer.zero();

    if ( emode == HeatTransferEM ) {
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
        _error("Unknown ElementMode");
    }
}


void
TransportElement :: computeBCSubVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, int indx)
{
    int n, id;
    GeneralBoundaryCondition *load;
    bcGeomType ltype;
    FloatArray vec;

    answer.resize( this->giveNumberOfDofManagers() );
    answer.zero();

    // loop over boundary load array
    int nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        id    = boundaryLoadArray.at(i * 2);
        load  = ( GeneralBoundaryCondition * ) domain->giveLoad(n);
        ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, ( Load * ) load, id, tStep, mode, indx);
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceBCSubVectorAt(vec, ( Load * ) load, id, tStep, mode, indx);
        } else {
            _error("computeBCSubVectorAt : unsupported bc type encountered");
        }

        answer.add(vec);
    }

    // end loop over applied bc
}

void
TransportElement :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge,
                                             TimeStep *tStep, ValueModeType mode, int indx)
{
    int approxOrder, numberOfEdgeIPs;

    answer.resize( this->giveNumberOfDofManagers() );
    answer.zero();

    if ( ( load->giveType() == TransmissionBC ) || ( load->giveType() == ConvectionBC ) ) {
        BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
        if ( edgeLoad->isDofExcluded(indx) ) {
            return;
        }


        approxOrder = edgeLoad->giveApproxOrder() + this->giveApproxOrder(indx);
        numberOfEdgeIPs = ( int ) ceil( ( approxOrder + 1. ) / 2. );
        GaussIntegrationRule iRule(1, this, 1, 1);
        iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);
        GaussPoint *gp;
        FloatArray reducedAnswer, val, ntf, n;
        IntArray mask;
        double dV, coeff = 1.0;

        if ( load->giveType() == TransmissionBC ) {
            coeff = -1.0;
        } else {
            coeff = edgeLoad->giveProperty('a');
        }

        for (int i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();
            this->computeEgdeNAt(n, *lcoords);
            dV  = this->computeEdgeVolumeAround(gp, iEdge);

            if ( edgeLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                edgeLoad->computeValueAt(val, tStep, *lcoords, mode);
            } else {
                FloatArray globalIPcoords;
                this->computeEdgeIpGlobalCoords(globalIPcoords, *lcoords, iEdge);
                edgeLoad->computeValueAt(val, tStep, globalIPcoords, mode);
            }

            reducedAnswer.add(val.at(indx) * coeff * dV, n);
        }

        this->giveEdgeDofMapping(mask, iEdge);
        answer.assemble(reducedAnswer, mask);
    } else {
        _error("computeBCSubVectorAt : unsupported bc type encountered");
    }
}

void
TransportElement :: computeSurfaceBCSubVectorAt(FloatArray &answer, Load *load,
                                                int iSurf, TimeStep *tStep, ValueModeType mode, int indx)
{
    int approxOrder;
    double dV, coeff = 1.0;

    if ( !this->testElementExtension(Element_SurfaceLoadSupport) ) {
        _error("computeSurfaceBCSubVectorAt : no surface load support");
    }

    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( surfLoad ) {
        IntegrationRule *iRule;
        GaussPoint *gp;
        FloatArray reducedAnswer, val, globalIPcoords, n;
        IntArray mask;

        answer.resize( this->giveNumberOfDofManagers() );
        answer.zero();

        if ( surfLoad->isDofExcluded(indx) ) {
            return;
        }

        if ( load->giveType() == TransmissionBC ) {
            coeff = -1.0;
        } else {
            coeff = surfLoad->giveProperty('a');
        }

        approxOrder = surfLoad->giveApproxOrder() + this->giveApproxOrder(indx);

        iRule = this->GetSurfaceIntegrationRule(approxOrder);
        for (int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            this->computeSurfaceNAt(n, *gp->giveCoordinates());
            dV  = this->computeSurfaceVolumeAround(gp, iSurf);

            if ( surfLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                surfLoad->computeValueAt(val, tStep, *gp->giveCoordinates(), mode);
            } else {
                this->computeSurfIpGlobalCoords(globalIPcoords, *gp->giveCoordinates(), iSurf);
                surfLoad->computeValueAt(val, tStep, globalIPcoords, mode);
            }

            reducedAnswer.add(val.at(indx) * coeff * dV, n);
        }

        this->giveSurfaceDofMapping(mask, iSurf);
        answer.assemble(reducedAnswer, mask);

        delete iRule;
    } else {
        _error("computeSurfaceBCSubVectorAt : unsupported bc type encountered");
    }
}




void
TransportElement :: computeBCSubMtrxAt(FloatMatrix &answer, TimeStep *tStep, ValueModeType mode, int indx)
{
    int k, id, defined = 0;
    GeneralBoundaryCondition *load;
    double dV;
    bcGeomType ltype;

    answer.resize( this->giveNumberOfDofManagers(), this->giveNumberOfDofManagers() );
    answer.zero();

    // loop over boundary load array
    int nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
    for (int i = 1; i <= nLoads; i++ ) {
        k     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        id    = boundaryLoadArray.at(i * 2);
        load  = ( Load * ) domain->giveLoad(k);
        if ( load->giveType() == ConvectionBC ) {
            ltype = load->giveBCGeoType();
            if ( ltype == EdgeLoadBGT ) {
                BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
                if ( edgeLoad->isDofExcluded(indx) ) {
                    continue;
                }

                defined = 1;
                int approxOrder = 2 * this->giveApproxOrder(indx);
                int numberOfEdgeIPs = ( int ) ceil( ( approxOrder + 1. ) / 2. );
                GaussIntegrationRule iRule(1, this, 1, 1);
                iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);
                GaussPoint *gp;
                FloatArray val, n;
                IntArray mask;
                FloatMatrix subAnswer;

                for (int igp = 0; igp < iRule.getNumberOfIntegrationPoints(); igp++ ) {
                    gp  = iRule.getIntegrationPoint(igp);
                    this->computeEgdeNAt(n, *gp->giveCoordinates());
                    dV  = this->computeEdgeVolumeAround(gp, id);
                    subAnswer.plusDyadSymmUpper(n, n, dV * edgeLoad->giveProperty('a') );
                }

                subAnswer.symmetrized();
                this->giveEdgeDofMapping(mask, id);
                answer.assemble(subAnswer, mask);
            } else if ( ltype == SurfaceLoadBGT ) {
                IntegrationRule *iRule;
                GaussPoint *gp;
                FloatArray val, n;
                IntArray mask;
                FloatMatrix subAnswer;

                BoundaryLoad *surfLoad = static_cast< BoundaryLoad * >(load);
                if ( surfLoad->isDofExcluded(indx) ) {
                    continue;
                }

                defined = 1;
                int approxOrder = 2 * this->giveApproxOrder(indx);
                iRule = this->GetSurfaceIntegrationRule(approxOrder);

                for (int igp = 0; igp < iRule->getNumberOfIntegrationPoints(); igp++ ) {
                    gp  = iRule->getIntegrationPoint(igp);
                    this->computeSurfaceNAt(n, *gp->giveCoordinates());
                    dV  = this->computeSurfaceVolumeAround(gp, id);
                    subAnswer.plusDyadSymmUpper( n, n, dV * surfLoad->giveProperty('a') );
                }

                delete iRule;
                subAnswer.symmetrized();
                this->giveSurfaceDofMapping(mask, id);
                answer.assemble(subAnswer, mask);
            } else {
                _error("computeBCSubMtrxAt : unsupported bc type encountered");
            }
        }
    }

    // end loop over applied bc

    if ( !defined ) {
        answer.resize(0, 0);
    }
}


void
TransportElement :: assembleLocalContribution(FloatMatrix &answer, FloatMatrix &src,
                                              int ndofs, int rdof, int cdof)
{
    int ti, tj;
    int nnodes = this->giveNumberOfDofManagers();

    for (int i = 1; i <= nnodes; i++ ) {
        ti = ( i - 1 ) * ndofs + rdof;
        for (int j = 1; j <= nnodes; j++ ) {
            tj = ( j - 1 ) * ndofs + cdof;
            answer.at(ti, tj) += src.at(i, j);
        }
    }
}


void
TransportElement :: assembleLocalContribution(FloatArray &answer, FloatArray &src,
                                              int ndofs, int rdof)
{
    int ti;
    int nnodes = this->giveNumberOfDofManagers();

    for (int i = 1; i <= nnodes; i++ ) {
        ti = ( i - 1 ) * ndofs + rdof;
        answer.at(ti) += src.at(i);
    }
}

void
TransportElement :: computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatArray r, br;
    FloatMatrix b, d;

    this->computeVectorOf(EID_ConservationEquation, VM_Total, stepN, r);
    this->computeGradientMatrixAt(b, gp);

    if ( emode == HeatTransferEM ) {
        this->computeConstitutiveMatrixAt(d, Conductivity_hh, gp, stepN);
        br.beProductOf(b, r);
        answer.beProductOf(d, br);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatArray r_h, r_w;
        FloatMatrix b_tot, d_hh, d_hw, d_wh, d_ww;
        r_h.resize(r.giveSize() / 2);
        r_h.zero();
        r_w.resize(r.giveSize() / 2);
        r_w.zero();

        this->computeConstitutiveMatrixAt(d_hh, Conductivity_hh, gp, stepN);
        this->computeConstitutiveMatrixAt(d_hw, Conductivity_hw, gp, stepN);
        this->computeConstitutiveMatrixAt(d_wh, Conductivity_wh, gp, stepN);
        this->computeConstitutiveMatrixAt(d_ww, Conductivity_ww, gp, stepN);
        d.resize( 2 * d_hh.giveNumberOfRows(), 2 * d_hh.giveNumberOfColumns() );
        d.zero();

        b_tot.resize( 2 * b.giveNumberOfRows(), 2 * b.giveNumberOfColumns() );
        b_tot.zero();
        b_tot.addSubMatrix(b, 1, 1);
        b_tot.addSubMatrix(b, b.giveNumberOfRows() + 1, b.giveNumberOfColumns() + 1);
        br.beProductOf(b_tot, r);

        d.addSubMatrix(d_hh, 1, 1);
        d.addSubMatrix(d_hw, 1, d_hh.giveNumberOfColumns() + 1);
        d.addSubMatrix(d_wh, d_hh.giveNumberOfRows() + 1, 1);
        d.addSubMatrix(d_ww, d_hh.giveNumberOfRows() + 1, d_wh.giveNumberOfColumns() + 1);

        answer.beProductOf(d, br);
    } else {
        OOFEM_ERROR1("Unknown element mode");
    }

    answer.negated();
}


void
TransportElement :: updateInternalState(TimeStep *stepN)
// Updates the receiver at the end of a solution step
{
    int i, j;
    IntegrationRule *iRule;
    FloatArray stateVector, r, br;
    FloatMatrix n, b, d;
    TransportMaterial *mat = ( ( TransportMaterial * ) this->giveMaterial() );
    GaussPoint *gp;

    // force updating ip values
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            this->computeNmatrixAt( n, *gp->giveCoordinates() );
            this->computeVectorOf(EID_ConservationEquation, VM_Total, stepN, r);
            stateVector.beProductOf(n, r);
            mat->updateInternalState(stateVector, gp, stepN);
        }
    }
}

int
TransportElement :: EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                          FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                          TimeStep *atTime)
{
    int indx;
    FloatArray elemvector, f, lc;
    FloatMatrix n;
    IntArray elemdofs;
    // determine element dof ids
    this->giveElementDofIDMask(pf.giveEquationID(), elemdofs);
    // first evaluate element unknown vector
    this->computeVectorOf(pf, mode, atTime, elemvector);
    // determine corresponding local coordinates
    if ( this->computeLocalCoordinates(lc, coords) ) {
        // compute interpolation matrix
        this->computeNmatrixAt(n, lc);
        // compute answer
        answer.resize( dofId.giveSize() );
        for (int i = 1; i <= dofId.giveSize(); i++ ) {
            if ( ( indx = elemdofs.findFirstIndexOf( dofId.at(i) ) ) ) {
                double sum = 0.0;
                for (int j = 1; j <= elemvector.giveSize(); j++ ) {
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
        _error("EIPrimaryFieldI_evaluateFieldVectorAt: target point not in receiver volume");
        return 1; // failed
    }
}


#ifdef __OOFEG
int
TransportElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                            int node, TimeStep *atTime)
{
    Node *n = this->giveNode(node);
    if ( type == IST_Temperature ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(T_f) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
            return 1;
        } else {
            return 0;
        }
    } else if ( type == IST_MassConcentration_1 ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(C_1) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

#endif

int
TransportElement :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_Temperature ) || ( type == IST_MassConcentration_1 ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return Element :: giveIntVarCompFullIndx(answer, type);
    }
}
} // end namespace oofem
