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

#include "structuralelement.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "nonlocalmaterialext.h"
#include "load.h"
#include "boundaryload.h"
#include "pointload.h"
#include "structtemperatureload.h"
#include "structeigenstrainload.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "mathfem.h"
#include "materialmapperinterface.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif


namespace oofem {
StructuralElement :: StructuralElement(int n, Domain *aDomain) :
    Element(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    activityLtf = 0;
    initialDisplacements = NULL;
}


StructuralElement :: ~StructuralElement()
// Destructor.
{
    if ( initialDisplacements ) {
        delete initialDisplacements;
    }
}


void
StructuralElement :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                 MatResponseMode rMode, GaussPoint *gp,
                                                 TimeStep *tStep)
// Returns the  material matrix {E} of the receiver.
// type of matrix is determined by this->giveMaterialMode()
// rMode parameter determines type of stiffness matrix to be requested
// (tangent, secant, ...)
{
    ( ( StructuralCrossSection * ) this->giveCrossSection() )
    ->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
}


void
StructuralElement :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body
// loads, at stepN.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    int i;
    double dens, dV;
    GaussPoint *gp;
    FloatArray force, ntf;
    FloatMatrix n, nt, T;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, stepN, mode);
    //force.times( this->giveMaterial()->give('d') );

    answer.resize(0);

    if ( force.giveSize() ) {
        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            this->computeNmatrixAt(gp, n);
            dV  = this->computeVolumeAround(gp);
            dens = this->giveMaterial()->give('d', gp);
            nt.beTranspositionOf(n);
            ntf.beProductOf(nt, force);
            ntf.times(dV * dens);
            answer.add(ntf);
        }
    } else {
        return;
    }
}

void
StructuralElement :: computePointLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    FloatArray force, coords, lcoords;
    FloatMatrix T, n;

    PointLoad *pointLoad = dynamic_cast< PointLoad * >(load);
    pointLoad->giveCoordinates(coords);
    pointLoad->computeValueAt(force, tStep, coords, mode);
    if ( this->computeLocalCoordinates(lcoords, coords) ) {
        GaussPoint __gp(NULL, 0, (new FloatArray(lcoords)), 1.0, _Unknown);
        this->computeNmatrixAt(& __gp, n);
        answer.beTProductOf(n, force);
    } else {
        _warning("computePointLoadVectorAt: point load outside element");
    }


    // transform force
    if ( pointLoad->giveCoordSystMode() == PointLoad :: PL_GlobalMode ) {
        // transform from global to element local c.s
        if ( this->computeLoadGToLRotationMtrx(T) ) {
            answer.rotatedWith(T, 'n');
        }
    }
}

void
StructuralElement :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load,
                                             int iEdge, TimeStep *tStep, ValueModeType mode)
{
    // computes numericaly the edge load vector of the receiver for given load
    // Each element edge must have unique number assigned to identify it.
    // Integration is done in local edge space (i.e. one dimensional integration is
    // performed on line). This general implementation requires that element must
    // provide following functions:
    // - ComputeEgdeNMatrixAt - returns interpolation matrix of the edge in the
    //   local edge space.
    // - computeEdgeVolumeAround - returns volumeAround in local edge space
    // - GiveEdgeDofMapping - returns integer array specifying local edge dof mapping to
    //   element dofs.
    //
    int i, approxOrder, numberOfGaussPoints;
    double dV;
    FloatMatrix T;
    FloatArray globalIPcoords;

    if ( !this->testElementExtension(Element_EdgeLoadSupport) ) {
        _error("computeEdgeLoadVectorAt : no edge load support");
    }

    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( edgeLoad ) {
        approxOrder = edgeLoad->giveApproxOrder() + this->giveApproxOrder();
        numberOfGaussPoints = ( int ) ceil( ( approxOrder + 1. ) / 2. );
        GaussIntegrationRule iRule(1, this, 1, 1);
        iRule.setUpIntegrationPoints(_Line, numberOfGaussPoints, _Unknown);
        GaussPoint *gp;
        FloatArray reducedAnswer, force, ntf;
        IntArray mask;
        FloatMatrix n, nt;

        for ( i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule.getIntegrationPoint(i);
            this->computeEgdeNMatrixAt(n, gp);
            dV  = this->computeEdgeVolumeAround(gp, iEdge);
            nt.beTranspositionOf(n);

            if ( edgeLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                edgeLoad->computeValueAt(force, tStep, * ( gp->giveCoordinates() ), mode);
            } else {
                this->computeEdgeIpGlobalCoords(globalIPcoords, gp, iEdge);
                edgeLoad->computeValueAt(force, tStep, globalIPcoords, mode);
            }

            // transform force
            if ( edgeLoad->giveCoordSystMode() == BoundaryLoad :: BL_GlobalMode ) {
                // transform from global to element local c.s
                if ( this->computeLoadGToLRotationMtrx(T) ) {
                    force.rotatedWith(T, 'n');
                }
            } else {
                // transform from local edge to element local c.s
                if ( this->computeLoadLEToLRotationMatrix(T, iEdge, gp) ) {
                    force.rotatedWith(T, 'n');
                }
            }

            ntf.beProductOf(nt, force);
            ntf.times(dV);
            reducedAnswer.add(ntf);
        }


        this->giveEdgeDofMapping(mask, iEdge);
        answer.resize( this->computeNumberOfDofs(EID_MomentumBalance) );
        answer.zero();
        answer.assemble(reducedAnswer, mask);

        return;
    } else {
        _error("computeEdgeLoadVectorAt: incompatible load");
        return;
    }
}





void
StructuralElement :: computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                                int iSurf, TimeStep *tStep, ValueModeType mode)
{
    // computes numericaly the surface load vector of the receiver for given load
    // Each element surface must have unique number assigned to identify it.
    // Integration is done in local surface space (i.e. two dimensional integration is
    // performed on triangle or square). This general implementation requires that element must
    // provide following functions:
    // - GetSurfaceIntegrationRule - returns integration rule for surface for given number of
    //   integration rules.
    // - ComputeSurfaceNMatrixAt - returns interpolation matrix of the surface in the
    //   local edge space.
    // - computeSurfaceVolumeAround - returns volumeAround in local edge space
    // - GiveSurfaceDofMapping - returns integer array specifying local edge dof mapping to
    //   element dofs.
    //

    int i, approxOrder;
    double dV;
    FloatMatrix T;

    if ( !this->testElementExtension(Element_SurfaceLoadSupport) ) {
        _error("computeSurfaceLoadVectorAt : no surface load support");
    }

    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( surfLoad ) {
        IntegrationRule *iRule;
        GaussPoint *gp;
        FloatArray reducedAnswer, force, ntf;
        IntArray mask;
        FloatMatrix n, nt;
        FloatArray globalIPcoords;

        approxOrder = surfLoad->giveApproxOrder() + this->giveApproxOrder();

        iRule = this->GetSurfaceIntegrationRule(approxOrder);
        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            this->computeSurfaceNMatrixAt(n, gp);
            dV  = this->computeSurfaceVolumeAround(gp, iSurf);
            nt.beTranspositionOf(n);

            if ( surfLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                surfLoad->computeValueAt(force, tStep, * ( gp->giveCoordinates() ), mode);
            } else {
                this->computeSurfIpGlobalCoords(globalIPcoords, gp, iSurf);
                surfLoad->computeValueAt(force, tStep, globalIPcoords, mode);
            }

            // transform force
            if ( surfLoad->giveCoordSystMode() == BoundaryLoad :: BL_GlobalMode ) {
                // transform from global to element local c.s
                if ( this->computeLoadGToLRotationMtrx(T) ) {
                    force.rotatedWith(T, 'n');
                }
            } else {
                // transform from local edge to element local c.s
                if ( this->computeLoadLSToLRotationMatrix(T, iSurf, gp) ) {
                    force.rotatedWith(T, 'n');
                }
            }

            ntf.beProductOf(nt, force);
            ntf.times(dV);
            reducedAnswer.add(ntf);
        }

        delete iRule;
        this->giveSurfaceDofMapping(mask, iSurf);
        answer.resize( this->computeNumberOfDofs(EID_MomentumBalance) );
        answer.zero();
        answer.assemble(reducedAnswer, mask);

        return;
    } else {
        _error("computeSurfaceLoadVectorAt: incompatible load");
        return;
    }
}


void
StructuralElement :: computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
//
// Computes numerically temperature and eigenstrain load vector
// Assumes that temperature is constant over the whole element
{
    int i;
    // TemperatureLoad   *load;
    double dV;
    GaussPoint *gp;
    FloatArray et, de, bde;
    FloatMatrix b, bt, d;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    StructuralCrossSection *cs = ( StructuralCrossSection * ) this->giveCrossSection();
    //   if (this -> giveBodyLoadArray() -> isEmpty())         // no loads
    //      return NULL ;

    //   else {
    // perform assembling of load vector over
    // complete volume
    answer.resize(0);
    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        cs->computeStressIndependentStrainVector(et, gp, tStep, mode);
        if ( et.giveSize() ) {
            this->computeBmatrixAt(gp, b);
            bt.beTranspositionOf(b);
            this->computeConstitutiveMatrixAt(d, TangentStiffness, gp, tStep);
            dV  = this->computeVolumeAround(gp);
            de.beProductOf(d, et);
            bde.beProductOf(bt, de);
            bde.times(dV);
            answer.add(bde);
        }
    }
}


void
StructuralElement :: computePrescribedStrainLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
//
// Computes numerically temperature load vector
// Assumes that temperature is constant over the
// whole element
//
{
    this->computePrescribedStrainLocalLoadVectorAt(answer, tStep, mode);
}


void
StructuralElement :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
// Computes numerically the consistent (full) mass matrix of the receiver.
{
    int i, nip, ndofs = computeNumberOfDofs(EID_MomentumBalance);
    double density, dV;
    FloatMatrix n;
    GaussPoint *gp;
    GaussIntegrationRule iRule(1, this, 1, 1);
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    if ( ( nip = this->giveNumberOfIPForMassMtrxIntegration() ) == 0 ) {
        _error("computeConsistentMassMatrix no integration points available");
    }

    iRule.setUpIntegrationPoints( this->giveIntegrationDomain(),
                                 nip, this->giveMaterialMode() );

    this->giveMassMtrxIntegrationgMask(mask);

    //density = this->giveMaterial()->give('d');
    mass = 0.;

    for ( i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
        gp      = iRule.getIntegrationPoint(i);
        this->computeNmatrixAt(gp, n);
        density = this->giveMaterial()->give('d', gp);
        dV      = this->computeVolumeAround(gp);
        mass   += density * dV;

        if ( mask.isEmpty() ) {
            answer.plusProductSymmUpper(n, n, density * dV);
        } else {
            int i, j, k;
            double summ;

            for ( i = 1; i <= ndofs; i++ ) {
                for ( j = i; j <= ndofs; j++ ) {
                    summ = 0.;
                    for ( k = 1; k <= n.giveNumberOfRows(); k++ ) {
                        if ( mask.at(k) == 0 ) {
                            continue;
                        }

                        summ += n.at(k, i) * n.at(k, j);
                    }

                    answer.at(i, j) += summ * density * dV;
                }
            }
        }
    }

    answer.symmetrized();
}


void
StructuralElement :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// Why is this function taken separately ?
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further subtract part corresponding to non-nodal loading.
{
    int i, n, id, nLoads;
    bcGeomType ltype;
    GeneralBoundaryCondition *load;
    FloatArray helpLoadVector;

    answer.resize(0);

    // loop over body load array first
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        n     = bodyLoadArray.at(i);
        load  = ( GeneralBoundaryCondition * ) domain->giveLoad(n);
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            this->computeBodyLoadVectorAt(helpLoadVector, ( Load * ) load, stepN, mode);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            if ( load->giveBCValType() != TemperatureBVT && load->giveBCValType() != EigenstrainBVT ) {
                // temperature and eigenstrain is handled separately at computeLoadVectorAt subroutine
                OOFEM_ERROR("StructuralElement :: computeForceLoadVector - unsupported load type class");
            }
        }
    }

    // loop over boundary load array
    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( i = 1; i <= nLoads; i++ ) {
        n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        id    = boundaryLoadArray.at(i * 2);
        load  = ( Load * ) domain->giveLoad(n);
        ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVectorAt(helpLoadVector, ( Load * ) load, id, stepN, mode);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVectorAt(helpLoadVector, ( Load * ) load, id, stepN, mode);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == PointLoadBGT ) {
            // id not used
            this->computePointLoadVectorAt(helpLoadVector, ( Load * ) load, stepN, mode);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            OOFEM_ERROR("StructuralElement :: computeForceLoadVector -unsupported load type class");
        }
    }
}


void
StructuralElement :: computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// Why is this function taken separately ?
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    this->computeLocalForceLoadVector(answer, stepN, mode);
}


void
StructuralElement :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    double mass;

    this->computeConsistentMassMatrix(answer, tStep, mass);
}

void
StructuralElement :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    double mass;

    IntArray nodeDofIDMask, dimFlag(3);
    IntArray nodalArray;
    int i, j, indx = 0, k, ldofs, dim;
    double summ;

    if ( !this->isActivated(tStep) ) {
        int ndofs = computeNumberOfDofs(EID_MomentumBalance);
        answer.resize(ndofs, ndofs);
        answer.zero();
        return;
    }

    this->computeConsistentMassMatrix(answer, tStep, mass);
    ldofs = answer.giveNumberOfRows();

    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, EID_MomentumBalance, nodeDofIDMask);
        //this->giveDofManager(i)->giveLocationArray(nodeDofIDMask, nodalArray);
        for ( j = 1; j <= nodeDofIDMask.giveSize(); j++ ) {
            indx++;
            // zero all off-diagonal terms
            for ( k = 1; k <= ldofs; k++ ) {
                if ( k != indx ) {
                    answer.at(indx, k) = 0.;
                    answer.at(k, indx) = 0.;
                }
            }

            if ( ( nodeDofIDMask.at(j) != D_u ) && ( nodeDofIDMask.at(j) != D_v ) && ( nodeDofIDMask.at(j) != D_w ) ) {
                // zero corresponding diagonal member too <= no displacement dof
                answer.at(indx, indx) = 0.;
            } else if ( nodeDofIDMask.at(j) == D_u ) {
                dimFlag.at(1) = 1;
            } else if ( nodeDofIDMask.at(j) == D_v ) {
                dimFlag.at(2) = 1;
            } else if ( nodeDofIDMask.at(j) == D_w ) {
                dimFlag.at(3) = 1;
            }
        }
    }

    if ( indx != ldofs ) {
        _error("computeMassMatrix : internal consistency check failed");
    }

    dim = dimFlag.at(1) + dimFlag.at(2) + dimFlag.at(3);
    for ( summ = 0., k = 1; k <= ldofs; k++ ) {
        summ += answer.at(k, k);
    }

    answer.times(dim * mass / summ);
}


void
StructuralElement :: computeResultingIPTemperatureAt(FloatArray &answer, TimeStep *stepN, GaussPoint *gp, ValueModeType mode)
// Computes at stepN the resulting force due to all temperature loads that act
// on the receiver. This force is used by the element for computing its
// body load vector.
{
    int i, n, nLoads;
    StructuralTemperatureLoad *load;
    FloatArray gCoords, temperature;

    if ( this->computeGlobalCoordinates( gCoords, * ( gp->giveCoordinates() ) ) == 0 ) {
        _error("computeResultingIPTemperatureAt: computeGlobalCoordinates failed");
    }

    answer.resize(0);
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        n     = bodyLoadArray.at(i);
        load  = ( StructuralTemperatureLoad * ) domain->giveLoad(n);
        if ( load->giveBCValType() == TemperatureBVT ) {
            load->computeValueAt(temperature, stepN, gCoords, mode);
            answer.add(temperature);
        }
    }
}

void
StructuralElement :: computeResultingIPEigenstrainAt(FloatArray &answer, TimeStep *stepN, GaussPoint *gp, ValueModeType mode)
// Computes at stepN all eigenstrains that act on the receiver. Eigenstrains are used by the element for computing its body load vector.
{
    int i, n, nLoads;
    StructuralEigenstrainLoad *load;
    FloatArray gCoords, eigenstrain;

    if ( this->computeGlobalCoordinates( gCoords, * ( gp->giveCoordinates() ) ) == 0 ) {
        _error("computeResultingIPTemperatureAt: computeGlobalCoordinates failed");
    }

    answer.resize(0);
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        n     = bodyLoadArray.at(i);
        load  = ( StructuralEigenstrainLoad * ) domain->giveLoad(n);
        if ( load->giveBCValType() == EigenstrainBVT ) {
            load->computeValueAt(eigenstrain, stepN, gCoords, mode);
            answer.add(eigenstrain);
        }
    }
}



void
StructuralElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                            TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    int i, j, k, iStartIndx, iEndIndx, jStartIndx, jEndIndx;
    double dV;
    FloatMatrix d, bi, bj, dbj, dij;
    GaussPoint *gp;
    IntegrationRule *iRule;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, this->material);

    answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    if ( numberOfIntegrationRules > 1 ) {
        for ( i = 0; i < numberOfIntegrationRules; i++ ) {
            iStartIndx = integrationRulesArray [ i ]->getStartIndexOfLocalStrainWhereApply();
            iEndIndx   = integrationRulesArray [ i ]->getEndIndexOfLocalStrainWhereApply();
            for ( j = 0; j < numberOfIntegrationRules; j++ ) {
                jStartIndx = integrationRulesArray [ j ]->getStartIndexOfLocalStrainWhereApply();
                jEndIndx   = integrationRulesArray [ j ]->getEndIndexOfLocalStrainWhereApply();
                if ( i == j ) {
                    iRule = integrationRulesArray [ i ];
                } else if ( integrationRulesArray [ i ]->getNumberOfIntegrationPoints() < integrationRulesArray [ j ]->getNumberOfIntegrationPoints() ) {
                    iRule = integrationRulesArray [ i ];
                } else {
                    iRule = integrationRulesArray [ j ];
                }

                for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
                    gp = iRule->getIntegrationPoint(k);
                    this->computeBmatrixAt(gp, bi, iStartIndx, iEndIndx);
                    this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
                    dij.beSubMatrixOf(d, iStartIndx, iEndIndx, jStartIndx, jEndIndx);
                    if ( i != j ) {
                        this->computeBmatrixAt(gp, bj, jStartIndx, jEndIndx);
                    } else {
                        bj = bi;
                    }

                    dV  = this->computeVolumeAround(gp);
                    dbj.beProductOf(dij, bj);
                    if ( matStiffSymmFlag ) {
                        answer.plusProductSymmUpper(bi, dbj, dV);
                    } else {
                        answer.plusProductUnsym(bi, dbj, dV);
                    }
                }
            }
        }
    } else { // numberOfIntegrationRules == 1
        iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            this->computeBmatrixAt(gp, bj);
            this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
            dV = this->computeVolumeAround(gp);
            dbj.beProductOf(d, bj);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(bj, dbj, dV);
            } else {
                answer.plusProductUnsym(bj, dbj, dV);
            }
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}

void StructuralElement :: computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer,
                                                                      MatResponseMode rMode, TimeStep *tStep)
{
    int ir, j;
    FloatMatrix temp, bj, d, dbj;
    IntegrationRule *iRule;
    GaussPoint *gp;
    int ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric( rMode, this->giveMaterial()->giveNumber() );
    IntArray irlocnum;
    double dV;

    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // loop over individual integration rules
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = integrationRulesArray [ ir ];
        // loop over individual integration points
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            this->computeBmatrixAt(gp, bj);
            this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

            dV = this->computeVolumeAround(gp);
            dbj.beProductOf(d, bj);
            if ( matStiffSymmFlag ) {
                m->plusProductSymmUpper(bj, dbj, dV);
            } else {
                m->plusProductUnsym(bj, dbj, dV);
            }
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
            answer.assemble(* m, irlocnum);
            m->resize(0, 0);
        }
    } // end loop over irules

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
StructuralElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;

    if ( !this->isActivated(stepN) ) {
        answer.resize( this->giveCrossSection()->giveIPValueSize(IST_StrainTensor, gp) );
        answer.zero();
        return;
    }

    this->computeBmatrixAt(gp, b);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(*initialDisplacements);
    }
    answer.beProductOf(b, u);
}


void
StructuralElement :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step stepN. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{
    /*
     * StructuralCrossSection* cs = (StructuralCrossSection*) this->giveCrossSection();
     * FloatArray totalEpsilon;
     * // FloatArray *help;
     *
     *
     * this->computeStrainVector(totalEpsilon, gp,stepN) ;
     * cs->giveRealStresses (answer, ReducedForm, gp,totalEpsilon,stepN);
     *
     * return ;
     */
    FloatArray Epsilon;
    StructuralCrossSection *cs = ( StructuralCrossSection * ) this->giveCrossSection();

    this->computeStrainVector(Epsilon, gp, stepN);
    cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
}


void
StructuralElement :: giveInternalForcesVector(FloatArray &answer,
                                              TimeStep *tStep, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    GaussPoint *gp;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix b, bt, R, GNT;
    FloatArray bs, TotalStressVector;
    double dV;

    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeBmatrixAt(gp, b);
        bt.beTranspositionOf(b);
        // TotalStressVector = gp->giveStressVector() ;
        if ( useUpdatedGpRecord == 1 ) {
            TotalStressVector = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )
                                ->giveStressVector();
        } else {
            this->computeStressVector(TotalStressVector, gp, tStep);
        }

        //
        // updates gp stress and strain record  acording to current
        // increment of displacement
        //
        if ( TotalStressVector.giveSize() == 0 ) {
            break;
        }

        //
        // now every gauss point has real stress vector
        //
        // compute nodal representation of internal forces using f = B^T*Sigma dV
        //
        dV  = this->computeVolumeAround(gp);
        bs.beProductOf(bt, TotalStressVector);
        bs.times(dV);
        answer.add(bs);
    }

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}



void
StructuralElement :: giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer,
                                                                   TimeStep *tStep, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    GaussPoint *gp;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule;

    FloatMatrix b, bt, R, GNT;
    int ir;
    FloatArray temp, bs, TotalStressVector;
    IntArray irlocnum;
    double dV;

    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    FloatArray *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // loop over individual integration rules
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = integrationRulesArray [ ir ];

        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);
            this->computeBmatrixAt(gp, b);
            bt.beTranspositionOf(b);
            // TotalStressVector = gp->giveStressVector() ;
            if ( useUpdatedGpRecord == 1 ) {
                TotalStressVector = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )
                                    ->giveStressVector();
            } else {
                this->computeStressVector(TotalStressVector, gp, tStep);
            }

            //
            // updates gp stress and strain record  acording to current
            // increment of displacement
            //
            if ( TotalStressVector.giveSize() == 0 ) {
                break;
            }

            //
            // now every gauss point has real stress vector
            //
            // compute nodal representation of internal forces using f = B^T*Sigma dV
            //
            dV  = this->computeVolumeAround(gp);
            bs.beProductOf(bt, TotalStressVector);
            bs.times(dV);

            m->add(bs);

            // localize irule contribution into element matrix
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
                answer.assemble(* m, irlocnum);
                m->resize(0, 0);
            }
        }
    }

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


void
StructuralElement :: giveCharacteristicMatrix(FloatMatrix &answer,
                                               CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else if ( mtrx == SecantStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, SecantStiffness, tStep);
    } else if ( mtrx == ElasticStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, ElasticStiffness, tStep);
    } else if ( mtrx == MassMatrix ) {
        this->computeMassMatrix(answer, tStep);
    } else if ( mtrx == LumpedMassMatrix ) {
        this->computeLumpedMassMatrix(answer, tStep);
    } else if ( mtrx == InitialStressMatrix ) {
        this->computeInitialStressMatrix(answer, tStep);
    } else {
        _error2( "giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}



void
StructuralElement :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                               TimeStep *tStep)
//
// returns characteristics vector of receiver according to mtrx
//
{
    if ( mtrx == ExternalForcesVector ) {
        this->computeForceLoadVector(answer, tStep, mode);
    } else if ( ( mtrx == InternalForcesVector ) && ( mode == VM_Total ) ) {
        this->giveInternalForcesVector(answer, tStep);
    } else if ( ( mtrx == LastEquilibratedInternalForcesVector ) && ( mode == VM_Total ) ) {
        /* here tstep is not relevant, we set useUpdatedGpRecord = 1
         * and this will cause to integrate internal forces using existing (nontemp, equlibrated) stresses in
         * statuses. Mainly used to compute reaction forces */
        this->giveInternalForcesVector(answer, tStep, 1);
    } else {
        _error2( "giveCharacteristicVector: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}

void
StructuralElement :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);

    // record initial displacement if element not active
    if ( activityLtf && !isActivated(tStep) ) {
        if ( !initialDisplacements ) {
            initialDisplacements = new FloatArray();
        }

        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, * initialDisplacements);
    }
}


void
StructuralElement :: updateInternalState(TimeStep *stepN)
// Updates the receiver at end of step.
{
    int i, j;
    IntegrationRule *iRule;
    FloatArray stress;

    // force updating strains & stresses
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            computeStressVector(stress, iRule->getIntegrationPoint(j), stepN);
        }
    }
}

void
StructuralElement :: updateBeforeNonlocalAverage(TimeStep *atTime)
// Updates the local material quantities before nonlocal averaging
{
    /*
     * Nonlocal material, related to receiver is also passed, because it is not possible to
     * ask receiver  for its material model  pointer using giveMaterial() service (returns Material type)
     * and then to cast this Material pointer into NonlocalMaterial pointer type,
     * because it is possible to cast only the pointer of derived class to pointer to base class.
     */
    int i, j;
    IntegrationRule *iRule;
    FloatArray epsilon;

#ifdef __PARALLEL_MODE
    if ( this->giveParallelMode() == Element_remote ) {
        return;
    }

#endif

    // not possible - produces wrong result
    //StructuralNonlocalMaterial* material = DYNAMIC_CAST(StructuralNonlocalMaterial, this->giveMaterial());
    StructuralNonlocalMaterialExtensionInterface *materialExt =
        ( StructuralNonlocalMaterialExtensionInterface * ) this->giveMaterial()->
        giveInterface(NonlocalMaterialExtensionInterfaceType);

    if ( !materialExt ) {
        return;             //_error ("updateBeforeNonlocalAverage: material with no StructuralNonlocalMaterial support");
    }

    // force updating local quantities
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            this->computeStrainVector(epsilon, iRule->getIntegrationPoint(j), atTime);
            // provide material local strain increment - as is provided to computeRealStresVector
            // allows to update internal vars to be averaged to new state
            materialExt->updateBeforeNonlocAverage(epsilon, iRule->getIntegrationPoint(j), atTime);
        }
    }
}


int
StructuralElement :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    if ( !dynamic_cast<StructuralMaterial*>(this->giveMaterial()) ) {
        _warning("checkConsistency : material without structural support");
        result = 0;
    }

    if ( !this->giveCrossSection()->testCrossSectionExtension(CS_StructuralCapability) ) {
        _warning("checkConsistency : cross-section without structural support");
        result = 0;
    }

    return result;
}

void
StructuralElement :: condense(FloatMatrix *stiff, FloatMatrix *mass, FloatArray *load, IntArray *what)
{
    /*
     * function for condensation of stiffness matrix and if requested of load vector and
     * initial stress or mass matrices (if mass and load arguments are nonzero (not NULL) then
     * their condensation is done).
     * Based on Rayleigh-Ritz method
     */
    int i, ii, j, k;
    int nkon = what->giveSize();
    int size = stiff->giveNumberOfRows();
    int ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    long double coeff, dii, lii = 0;
    FloatArray *gaussCoeff = NULL;
    // stored gauss coefficient
    // for mass condensation

    // check
    if ( !stiff->isSquare() ) {
        _error("condense: stiffness size mismatch");
    }

    if ( mass ) {
        if ( !( mass->isSquare() && mass->giveNumberOfRows() == size ) ) {
            _error("condense: mass size mismatch");
        }
    }

    if ( load ) {
        if ( !( load->giveSize() == size ) ) {
            _error("condense: load size mismatch");
        }
    }

    // create gauss coeff array if mass condensation requested
    if ( mass ) {
        gaussCoeff = new FloatArray(size);
    }

    for ( i = 1; i <= nkon; i++ ) {
        ii  = what->at(i);
        if ( ( ii > ndofs ) || ( ii <= 0 ) ) {
            _error("condense: wrong dof number");
        }

        dii = stiff->at(ii, ii);
        if ( load ) {
            lii = load->at(ii);
        }

        // stiffness matrix condensation
        for ( j = 1; j <= size; j++ ) {
            coeff = -stiff->at(j, ii) / dii;
            if ( ii != j ) {
                for ( k = 1; k <= size; k++ ) {
                    stiff->at(j, k) += stiff->at(ii, k) * coeff;
                }
            }

            if ( load ) {
                load->at(j) += coeff * lii;
            }

            if ( mass ) {
                gaussCoeff->at(j) = coeff;
            }
        }

        for ( k = 1; k <= size; k++ ) {
            stiff->at(ii, k) = 0.;
            stiff->at(k, ii) = 0.;
        }

        // mass or initial stress matrix condensation
        // uses gauss coefficients stored in gaussCoeff.
        if ( mass ) {
            for ( j = 1; j <= size; j++ ) {
                for ( k = 1; k <= size; k++ ) {
                    if ( ( ii != j ) && ( ii != k ) ) {
                        mass->at(j, k) += mass->at(j, ii) * gaussCoeff->at(k) +
                                          mass->at(ii, k) * gaussCoeff->at(j) +
                                          mass->at(ii, ii) * gaussCoeff->at(j) * gaussCoeff->at(k);
                    }
                }
            }

            for ( k = 1; k <= size; k++ ) {
                mass->at(ii, k) = 0.;
                mass->at(k, ii) = 0.;
            }
        }
    }

    if ( mass ) {
        delete gaussCoeff;
    }
}


int
StructuralElement :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    return Element :: giveIPValue(answer, aGaussPoint, type, atTime);
}


void
StructuralElement :: giveNonlocalLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s)
{
    NonlocalMaterialStiffnessInterface *interface;
    // test for material model interface
    interface = ( NonlocalMaterialStiffnessInterface * )
                this->giveMaterial()->giveInterface(NonlocalMaterialStiffnessInterfaceType);
    if ( interface == NULL ) {
        locationArray.resize(0);
        return;
    } else {
        IntArray elemLocArry;
        // create lit of remote elements, contributing to receiver
        dynaList< localIntegrationRecord > *integrationDomainList;
        IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        dynaList< localIntegrationRecord > :: iterator pos;

        locationArray.resize(0);
        // loop over element IP
        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            integrationDomainList = interface->
                                    NonlocalMaterialStiffnessInterface_giveIntegrationDomainList( iRule->getIntegrationPoint(i) );
            // loop over IP influencing IPs, extract corresponding element numbers and their code numbers
            for ( pos = integrationDomainList->begin(); pos != integrationDomainList->end(); ++pos ) {
                ( * pos ).nearGp->giveElement()->giveLocationArray(elemLocArry, EID_MomentumBalance, s);
                /*
                 * Currently no care given to multiple occurences of code number in locationArray.
                 */
                locationArray.followedBy(elemLocArry, 20);
            }
        } // end loop over IPs

    }
}


void
StructuralElement :: addNonlocalStiffnessContributions(SparseMtrx &dest, const UnknownNumberingScheme &s, TimeStep *atTime)
{
    ///@todo Take into account cross section model (slaves)
    NonlocalMaterialStiffnessInterface *interface;

    if ( !this->isActivated(atTime) ) {
        return;
    }

    // test for material model interface
    interface = ( NonlocalMaterialStiffnessInterface * )
                this->giveMaterial()->giveInterface(NonlocalMaterialStiffnessInterfaceType);
    if ( interface == NULL ) {
        return;
    } else {
        IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        // loop over element IP
        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            interface->NonlocalMaterialStiffnessInterface_addIPContribution(dest, s, iRule->getIntegrationPoint(i), atTime);
        }
    }
}


int
StructuralElement :: adaptiveUpdate(TimeStep *tStep)
{
    int i, j, result = 1;
    IntegrationRule *iRule;
    FloatArray strain;

    MaterialModelMapperInterface *interface = ( MaterialModelMapperInterface * )
                                              this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType);

    if ( !interface ) {
        return 0;
    }

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            this->computeStrainVector(strain, iRule->getIntegrationPoint(j), tStep);
            result &= interface->MMI_update(iRule->getIntegrationPoint(j), tStep, & strain);
        }
    }

    return result;
}

IRResultType
StructuralElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                          // Required by IR_GIVE_FIELD macro

    result = Element :: initializeFrom(ir);

    activityLtf = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, activityLtf, IFT_StructuralElement_activityltf, "activityltf"); // Macro
    return result;
}



#ifdef __OOFEG

//
int
StructuralElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                             int node, TimeStep *atTime)
{
    if ( type == IST_DisplacementVector ) {
        Node *n = this->giveNode(node);
        answer.resize(3);
        answer.at(1) = n->giveUpdatedCoordinate(1, atTime, EID_MomentumBalance) - n->giveCoordinate(1);
        answer.at(2) = n->giveUpdatedCoordinate(2, atTime, EID_MomentumBalance) - n->giveCoordinate(2);
        answer.at(3) = n->giveUpdatedCoordinate(3, atTime, EID_MomentumBalance) - n->giveCoordinate(3);
        return 1;
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}



void
StructuralElement :: showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime)
{
    if ( ( mtrx == StiffnessMatrix ) || ( mtrx == TangentStiffnessMatrix ) ||
        ( mtrx == SecantStiffnessMatrix ) || ( mtrx == ElasticStiffnessMatrix ) ) {
        int i, j, n;
        IntArray loc;
        this->giveLocationArray( loc, EID_MomentumBalance, EModelDefaultEquationNumbering() );

        WCRec p [ 4 ];
        GraphicObj *go;

        EASValsSetLineWidth(OOFEG_SPARSE_PROFILE_WIDTH);
        EASValsSetColor( gc.getStandardSparseProfileColor() );
        EASValsSetLayer(OOFEG_SPARSE_PROFILE_LAYER);
        EASValsSetFillStyle(FILL_SOLID);

        n = loc.giveSize();
        for ( i = 1; i <= n; i++ ) {
            if ( loc.at(i) == 0 ) {
                continue;
            }

            for ( j = i; j <= n; j++ ) {
                if ( loc.at(j) == 0 ) {
                    continue;
                }

                if ( gc.getSparseProfileMode() == 0 ) {
                    p [ 0 ].x = ( FPNum ) loc.at(i) - 0.5;
                    p [ 0 ].y = ( FPNum ) loc.at(j) - 0.5;
                    p [ 0 ].z = 0.;
                    p [ 1 ].x = ( FPNum ) loc.at(i) + 0.5;
                    p [ 1 ].y = ( FPNum ) loc.at(j) - 0.5;
                    p [ 1 ].z = 0.;
                    p [ 2 ].x = ( FPNum ) loc.at(i) + 0.5;
                    p [ 2 ].y = ( FPNum ) loc.at(j) + 0.5;
                    p [ 2 ].z = 0.;
                    p [ 3 ].x = ( FPNum ) loc.at(i) - 0.5;
                    p [ 3 ].y = ( FPNum ) loc.at(j) + 0.5;
                    p [ 3 ].z = 0.;
                    go =  CreateQuad3D(p);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | LAYER_MASK, go);
                    EMAddGraphicsToModel(ESIModel(), go);

                    p [ 0 ].x = ( FPNum ) loc.at(j) - 0.5;
                    p [ 0 ].y = ( FPNum ) loc.at(i) - 0.5;
                    p [ 0 ].z = 0.;
                    p [ 1 ].x = ( FPNum ) loc.at(j) + 0.5;
                    p [ 1 ].y = ( FPNum ) loc.at(i) - 0.5;
                    p [ 1 ].z = 0.;
                    p [ 2 ].x = ( FPNum ) loc.at(j) + 0.5;
                    p [ 2 ].y = ( FPNum ) loc.at(i) + 0.5;
                    p [ 2 ].z = 0.;
                    p [ 3 ].x = ( FPNum ) loc.at(j) - 0.5;
                    p [ 3 ].y = ( FPNum ) loc.at(i) + 0.5;
                    p [ 3 ].z = 0.;
                    go =  CreateQuad3D(p);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | LAYER_MASK, go);
                    EMAddGraphicsToModel(ESIModel(), go);
                } else {
                    p [ 0 ].x = ( FPNum ) loc.at(i);
                    p [ 0 ].y = ( FPNum ) loc.at(j);
                    p [ 0 ].z = 0.;

                    EASValsSetMType(SQUARE_MARKER);
                    go = CreateMarker3D(p);
                    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | VECMTYPE_MASK, go);
                    EMAddGraphicsToModel(ESIModel(), go);

                    p [ 0 ].x = ( FPNum ) loc.at(j);
                    p [ 0 ].y = ( FPNum ) loc.at(i);
                    p [ 0 ].z = 0.;

                    EASValsSetMType(SQUARE_MARKER);
                    go = CreateMarker3D(p);
                    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | VECMTYPE_MASK, go);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    }
}

void
StructuralElement :: showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime)
{
    NonlocalMaterialStiffnessInterface *interface;
    int i;
    if ( ( ( mtrx == StiffnessMatrix ) || ( mtrx == TangentStiffnessMatrix ) ) ) {
        interface = ( NonlocalMaterialStiffnessInterface * )
                    this->giveMaterial()->giveInterface(NonlocalMaterialStiffnessInterfaceType);
        if ( interface == NULL ) {
            return;
        }

        IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        // loop over element IP
        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            interface->NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(iRule->getIntegrationPoint(i), gc, atTime);
        }
    }
}

#endif
} // end namespace oofem
