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

#include "structuralelementevaluator2.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
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
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "verbose.h"
#include "elementgeometry.h"
#include "feinterpol.h"
#include "elementside.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
 #include <stdlib.h>
 #include <stdio.h>
 #include <math.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

#include "materialmapperinterface.h"


namespace oofem {

StructuralElementEvaluator2 :: StructuralElementEvaluator2()
// Constructor. Creates an element evaluator
{
    initialDisplacements = NULL;
	displacementInterpolationNumber = 1;
	
}
	
	
StructuralElementEvaluator2 :: StructuralElementEvaluator2(int displInterpNumber)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    initialDisplacements = NULL;
	displacementInterpolationNumber = displInterpNumber;
}


StructuralElementEvaluator2 :: ~StructuralElementEvaluator2()
// Destructor.
{
    if ( initialDisplacements ) {
        delete initialDisplacements;
    }
}


void
StructuralElementEvaluator2 :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                 MatResponseMode rMode, GaussPoint *gp,
												 TimeStep *tStep, ElementGeometry *elemGeometry)
// Returns the  material matrix {E} of the receiver.
// type of matrix is determined by this->giveMaterialMode()
// rMode parameter determines type of stiffness matrix to be requested
// (tangent, secant, ...)
{
    this->giveStructuralCrossSection(elemGeometry)->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
}


void StructuralElementEvaluator2 :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep, ElementGeometry* elemGeometry)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }
    // Just a wrapper for the deadweight body load computations:
    PointLoad *p = dynamic_cast< PointLoad * >( load );
    if ( p ) {
        FloatArray gcoords, lcoords;
        p->giveCoordinates(gcoords);
        if ( elemGeometry->computeLocalCoordinates(lcoords, gcoords) ) {
            this->computePointLoadVectorAt(answer, load, tStep, mode, elemGeometry);
        }
    } else {
        ///@todo This assumption of dead-weight loads needs to be lifted. We can have other body loads, such as
        this->computeBodyLoadVectorAt(answer, load, tStep, mode, elemGeometry);
    }
}

void StructuralElementEvaluator2 :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, ElementGeometry* elemGeometry)
{
	
    answer.resize(0);
    if ( type != ExternalForcesVector ) {
        return;
    }

	FEInterpolation *fei = elemGeometry->giveInterpolation();
    if ( !fei ) {
        OOFEM_SIMPLE_ERROR("StructuralElement :: computeBoundaryLoadVector - No interpolator available\n");
    }

    FloatArray n_vec;
    FloatMatrix n, T;
    FloatArray force, globalIPcoords;
    int nsd = fei->giveNsd();

    IntegrationRule *iRule = fei->giveBoundaryIntegrationRule(load->giveApproxOrder(), boundary);

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray &lcoords = * gp->giveCoordinates();

        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValueAt(force, tStep, lcoords, mode);
        } else {
            fei->boundaryLocal2Global( globalIPcoords, boundary, lcoords, FEIElementGeometryWrapper(elemGeometry) );
            load->computeValueAt(force, tStep, globalIPcoords, mode);
        }

        ///@todo Make sure this part is correct.
        // We always want the global values in the end, so we might as well compute them here directly:
        // transform force
        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            // then just keep it in global c.s
        } else {
            ///@todo Support this...
            // transform from local boundary to element local c.s
            /*if ( this->computeLoadLSToLRotationMatrix(T, boundary, gp) ) {
             *  force.rotatedWith(T, 'n');
             * }*/
            // then to global c.s
            if ( this->computeLoadGToLRotationMtrx(T) ) {
                force.rotatedWith(T, 't');
            }
        }

        // Construct n-matrix
        fei->boundaryEvalN( n_vec, boundary, lcoords, FEIElementGeometryWrapper(elemGeometry) );
        n.beNMatrixOf(n_vec, nsd);

        ///@todo Some way to ask for the thickness at a global coordinate maybe?
        double thickness = 1.0; // Should be the circumference for axisymm-elements.
        double dV = thickness * gp->giveWeight() * fei->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(elemGeometry) );
        answer.plusProduct(n, force, dV);
    }

    delete iRule;
}


void
StructuralElementEvaluator2 :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemGeometry)
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV;
    GaussPoint *gp;
    FloatArray force, ntf;
    FloatMatrix n, T;
	IntegrationRule *iRule = elemGeometry->giveDefaultIntegrationRulePtr();

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
       // _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);
    // transform from global to element local c.s
    if ( this->computeLoadGToLRotationMtrx(T) ) {
        force.rotatedWith(T, 'n');
    }

    answer.resize(0);

    if ( force.giveSize() ) {
        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
			this->computeNmatrixAt(* ( gp->giveLocalCoordinates() ), n, elemGeometry);
            dV  = elemGeometry->computeVolumeAround(gp);
            dens = elemGeometry->giveCrossSection()->give('d', gp);
            ntf.beTProductOf(n, force);
            answer.add(dV * dens, ntf);
        }
    } else {
        return;
    }
}

void
StructuralElementEvaluator2 :: computePointLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemGeometry)
{
    FloatArray force, coords, lcoords;
    FloatMatrix T, n;

    PointLoad *pointLoad = dynamic_cast< PointLoad * >( load );
    pointLoad->giveCoordinates(coords);
    pointLoad->computeValueAt(force, tStep, coords, mode);
    if ( elemGeometry->computeLocalCoordinates(lcoords, coords) ) {
        GaussPoint __gp(NULL, 0, ( new FloatArray(lcoords) ), 1.0, _Unknown);
        this->computeNmatrixAt(* ( __gp.giveLocalCoordinates() ), n, elemGeometry);
        answer.beTProductOf(n, force);
    } else {
       // _warning("computePointLoadVectorAt: point load outside element");
    }

    // transform force
    if ( pointLoad->giveCoordSystMode() == Load :: CST_Global ) {
        // transform from global to element local c.s
        if ( this->computeLoadGToLRotationMtrx(T) ) {
            answer.rotatedWith(T, 'n');
        }
    }
}

void
StructuralElementEvaluator2 :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load,
                                             int iEdge, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemGeometry)
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
    int approxOrder, numberOfGaussPoints;
    double dV;
    FloatMatrix T;
    FloatArray globalIPcoords;

	
    if ( !this->testElementExtension(Element_EdgeLoadSupport) ) {
        //_error("computeEdgeLoadVectorAt : no edge load support");
    }

    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >( load );
    if ( edgeLoad ) {
        approxOrder = edgeLoad->giveApproxOrder() + this->giveApproxOrder();
        numberOfGaussPoints = ( int ) ceil( ( approxOrder + 1. ) / 2. );
        GaussIntegrationRule iRule(1, elemGeometry, 1, 1);
        iRule.SetUpPointsOnLine(numberOfGaussPoints, _Unknown);
        GaussPoint *gp;
        FloatArray reducedAnswer, force, ntf;
        IntArray mask;
        FloatMatrix n;

        for ( int i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule.getIntegrationPoint(i);
            this->computeEgdeNMatrixAt(n, iEdge, gp,elemGeometry);
            dV  = this->computeEdgeVolumeAround(gp, iEdge,elemGeometry);

            if ( edgeLoad->giveFormulationType() == Load :: FT_Entity ) {
                edgeLoad->computeValueAt(force, tStep, * ( gp->giveCoordinates() ), mode);
            } else {
                this->computeEdgeIpGlobalCoords(globalIPcoords, gp, iEdge,elemGeometry);
                edgeLoad->computeValueAt(force, tStep, globalIPcoords, mode);
            }

            // transform force
            if ( edgeLoad->giveCoordSystMode() == Load :: CST_Global ) {
                // transform from global to element local c.s
                if ( this->computeLoadGToLRotationMtrx(T) ) {
                    force.rotatedWith(T, 'n');
                }
            } else {
                // transform from local edge to element local c.s
                if ( this->computeLoadLEToLRotationMatrix(T, iEdge, gp,elemGeometry) ) {
                    force.rotatedWith(T, 'n');
                }
            }

            ntf.beTProductOf(n, force);
            reducedAnswer.add(dV, ntf);
        }


        this->giveEdgeDofMapping(mask, iEdge, elemGeometry);
        answer.resize( this->computeNumberOfDofs() );
        answer.zero();
        answer.assemble(reducedAnswer, mask);

        return;
    } else {
        //_error("computeEdgeLoadVectorAt: incompatible load");
        return;
    }
}


void
StructuralElementEvaluator2 :: computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                                int iSurf, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemGeometry)
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

    int approxOrder;
    double dV;
    FloatMatrix T;

    if ( !this->testElementExtension(Element_SurfaceLoadSupport) ) {
        //_error("computeSurfaceLoadVectorAt : no surface load support");
    }

    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >( load );
    if ( surfLoad ) {
        IntegrationRule *iRule;
        GaussPoint *gp;
        FloatArray reducedAnswer, force, ntf;
        IntArray mask;
        FloatMatrix n;
        FloatArray globalIPcoords;

        approxOrder = surfLoad->giveApproxOrder() + this->giveApproxOrder();

        iRule = this->GetSurfaceIntegrationRule(approxOrder,elemGeometry);
        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            this->computeSurfaceNMatrixAt(n, iSurf, gp,elemGeometry);
            dV  = this->computeSurfaceVolumeAround(gp, iSurf,elemGeometry);

            if ( surfLoad->giveFormulationType() == Load :: FT_Entity ) {
                surfLoad->computeValueAt(force, tStep, * ( gp->giveCoordinates() ), mode);
            } else {
                this->computeSurfIpGlobalCoords(globalIPcoords, gp, iSurf,elemGeometry);
                surfLoad->computeValueAt(force, tStep, globalIPcoords, mode);
            }

            // transform force
            if ( surfLoad->giveCoordSystMode() == Load :: CST_Global ) {
                // transform from global to element local c.s
                if ( this->computeLoadGToLRotationMtrx(T) ) {
                    force.rotatedWith(T, 'n');
                }
            } else {
                // transform from local edge to element local c.s
                if ( this->computeLoadLSToLRotationMatrix(T, iSurf, gp,elemGeometry) ) {
                    force.rotatedWith(T, 'n');
                }
            }

            ntf.beTProductOf(n, force);
            reducedAnswer.add(dV, ntf);
        }

        delete iRule;
        this->giveSurfaceDofMapping(mask, iSurf,elemGeometry);
        answer.resize( this->computeNumberOfDofs() );
        answer.zero();
        answer.assemble(reducedAnswer, mask);

        return;
    } else {
        //_error("computeSurfaceLoadVectorAt: incompatible load");
        return;
    }
}


void
StructuralElementEvaluator2 :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass,  ElementGeometry* elemGeometry, const double *ipDensity)
// Computes numerically the consistent (full) mass matrix of the receiver.
{
    int nip, ndofs = this->computeNumberOfDofs();
    double density, dV;
    FloatMatrix n;
    GaussIntegrationRule iRule(1, elemGeometry, 1, 1);
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !elemGeometry->isActivated(tStep) ) {
        return;
    }

    if ( ( nip = this->giveNumberOfIPForMassMtrxIntegration() ) == 0 ) {
        //_error("computeConsistentMassMatrix no integration points available");
    }

    iRule.setUpIntegrationPoints( elemGeometry->giveIntegrationDomain(),
                                  nip, elemGeometry->giveMaterialMode() );

    this->giveMassMtrxIntegrationgMask(mask);

    mass = 0.;

    for ( int i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule.getIntegrationPoint(i);
		this->computeNmatrixAt(* ( gp->giveLocalCoordinates() ), n, elemGeometry);
        density = elemGeometry->giveCrossSection()->give('d', gp);

        if ( ipDensity != NULL ) {
            // Override density if desired
            density = * ipDensity;
        }

        dV = elemGeometry->computeVolumeAround(gp);
        mass += density * dV;

        if ( mask.isEmpty() ) {
            answer.plusProductSymmUpper(n, n, density * dV);
        } else {
            for ( int i = 1; i <= ndofs; i++ ) {
                for ( int j = i; j <= ndofs; j++ ) {
                    double summ = 0.;
                    for ( int k = 1; k <= n.giveNumberOfRows(); k++ ) {
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
StructuralElementEvaluator2 :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemGeometry)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// Why is this function taken separately ?
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further subtract part corresponding to non-nodal loading.
{
    FloatArray helpLoadVector(1);
    answer.resize(0);
	Domain* domain = elemGeometry->giveDomain();
	IntArray *bodyLoadArray = this->giveBodyLoadArray();

    // loop over body load array first
    int nBodyLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nBodyLoads; i++ ) {
        int id = bodyLoadArray->at(i);
        Load *load = domain->giveLoad(id);
        bcGeomType ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            this->computeBodyLoadVectorAt(helpLoadVector, load, tStep, mode,elemGeometry);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            if ( load->giveBCValType() != TemperatureBVT && load->giveBCValType() != EigenstrainBVT ) {
                // temperature and eigenstrain is handled separately at computeLoadVectorAt subroutine
                OOFEM_SIMPLE_ERROR("StructuralElement :: computeLocalForceLoadVector - body load %d is of unsupported type (%d)", id, ltype);
            }
        }
    }

    // loop over boundary load array
    int nBoundaryLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nBoundaryLoads; i++ ) {
		IntArray *boundaryLoadArray = this->giveBoundaryLoadArray();
        int n = boundaryLoadArray->at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray->at(i * 2);
        Load *load = domain->giveLoad(n);
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
			this->computeEdgeLoadVectorAt(helpLoadVector, load, id, tStep, mode, elemGeometry);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVectorAt(helpLoadVector, load, id, tStep, mode, elemGeometry);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == PointLoadBGT ) {
            // id not used
            this->computePointLoadVectorAt(helpLoadVector, load, tStep, mode, elemGeometry);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
	  OOFEM_SIMPLE_ERROR("StructuralElement :: computeLocalForceLoadVector - boundary load %d is of unsupported type (%d)", id, ltype);
        }
    }
}


void
StructuralElementEvaluator2 :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemGeometry)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// Why is this function taken separately ?
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    this->computeLocalForceLoadVector(answer, tStep, mode, elemGeometry);
}


void
StructuralElementEvaluator2 :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep, ElementGeometry* elemGeometry)
// Returns the lumped mass matrix of the receiver.
{
    double mass;

    this->computeConsistentMassMatrix(answer, tStep, mass, elemGeometry);
}

void
StructuralElementEvaluator2 :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep, ElementGeometry* elemGeometry)
// Returns the lumped mass matrix of the receiver.
{
    double mass;

    IntArray nodeDofIDMask, dimFlag(3);
    IntArray nodalArray;
    int indx = 0, k, ldofs, dim;
    double summ;

	if ( !elemGeometry->isActivated(tStep) ) {
        int ndofs = this->computeNumberOfDofs();
        answer.resize(ndofs, ndofs);
        answer.zero();
        return;
    }

    this->computeConsistentMassMatrix(answer, tStep, mass, elemGeometry);
    ldofs = answer.giveNumberOfRows();
	
    for ( int i = 1; i <= elemGeometry->giveNumberOfDofManagers(); i++ ) {
		elemGeometry->giveDofManDofIDMask(i, EID_MomentumBalance, nodeDofIDMask);
        //this->giveDofManager(i)->giveLocationArray(nodeDofIDMask, nodalArray);
        for ( int j = 1; j <= nodeDofIDMask.giveSize(); j++ ) {
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
        //_error("computeMassMatrix : internal consistency check failed");
    }

    dim = dimFlag.at(1) + dimFlag.at(2) + dimFlag.at(3);
    summ = 0.;
    for ( int k = 1; k <= ldofs; k++ ) {
        summ += answer.at(k, k);
    }

    answer.times(dim * mass / summ);
}


void
StructuralElementEvaluator2 :: computeResultingIPTemperatureAt(FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode, ElementGeometry* elemGeometry)
// Computes at tStep the resulting force due to all temperature loads that act
// on the receiver. This force is used by the element for computing its
// body load vector.
{
    int n, nLoads;
    Load *load;
    FloatArray gCoords, temperature;
	Domain* domain = elemGeometry->giveDomain();


    if ( elemGeometry->computeGlobalCoordinates( gCoords, * ( gp->giveCoordinates() ) ) == 0 ) {
       // _error("computeResultingIPTemperatureAt: computeGlobalCoordinates failed");
    }

    answer.resize(0);
    nLoads = this->giveBodyLoadArray()->giveSize();
	IntArray* bodyLoadArray = this->giveBodyLoadArray();
    for ( int i = 1; i <= nLoads; i++ ) {
        n = bodyLoadArray->at(i);
        load = domain->giveLoad(n);
        if ( load->giveBCValType() == TemperatureBVT ) {
            static_cast< StructuralTemperatureLoad * >( load )->computeValueAt(temperature, tStep, gCoords, mode);
            answer.add(temperature);
        }
    }
}

void
StructuralElementEvaluator2 :: computeResultingIPEigenstrainAt(FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode, ElementGeometry* elemGeometry)
// Computes at tStep all eigenstrains that act on the receiver. Eigenstrains are used by the element for computing its body load vector.
{
    int n, nLoads;
    Load *load;
    FloatArray gCoords, eigenstrain;
	Domain* domain = elemGeometry->giveDomain();


    if ( elemGeometry->computeGlobalCoordinates( gCoords, * ( gp->giveCoordinates() ) ) == 0 ) {
        //_error("computeResultingIPTemperatureAt: computeGlobalCoordinates failed");
    }

    answer.resize(0);
    nLoads = this->giveBodyLoadArray()->giveSize();
	IntArray* bodyLoadArray = this->giveBodyLoadArray();
    for ( int i = 1; i <= nLoads; i++ ) {
        n = bodyLoadArray->at(i);
        load = domain->giveLoad(n);
        if ( load->giveBCValType() == EigenstrainBVT ) {
            static_cast< StructuralEigenstrainLoad * >( load )->computeValueAt(eigenstrain, tStep, gCoords, mode);
            answer.add(eigenstrain);
        }
    }
}



void
StructuralElementEvaluator2 :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                            TimeStep *tStep, ElementGeometry* elemGeometry)
// Computes numerically the stiffness matrix of the receiver.
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection(elemGeometry);
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.resize(0, 0);

    if ( !elemGeometry->isActivated(tStep) ) {
        return;
    }

    // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
    if ( elemGeometry->giveNumberOfIntegrationRules() == 1 ) {
        FloatMatrix B, D, DB;
        IntegrationRule *iRule = elemGeometry->giveDefaultIntegrationRulePtr() ;
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(j);

            // Engineering (small strain) stiffness
            if ( nlGeometry == 0 ) {
				this->computeBmatrixAt(gp, B, elemGeometry);
                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep, elemGeometry);
			// geometrical nonlinearity
            } else if ( nlGeometry == 1 ) {
                if ( elemGeometry->giveDomain()->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                    this->computeBmatrixAt(gp, B, elemGeometry);
                    /// @todo We probably need overloaded function (like above) here as well.
                    cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                } else { // Material stiffness dP/dF
                    this->computeBHmatrixAt(gp, B, elemGeometry);
                    /// @todo We probably need overloaded function (like above) here as well.
                    cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                }
            }

            double dV = elemGeometry->computeVolumeAround(gp);
            DB.beProductOf(D, B);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(B, DB, dV);
            } else {
                answer.plusProductUnsym(B, DB, dV);
            }
        }
		// geometrical nonlinearity - updated lagrangian formulation
        if ( elemGeometry->giveDomain()->giveEngngModel()->giveFormulation() == AL ) {
            FloatMatrix initialStressMatrix;
            this->computeInitialStressMatrix(initialStressMatrix, tStep, elemGeometry);
            answer.add(initialStressMatrix);
        }
    } else { /// @todo Verify that it works with large deformations
        int iStartIndx, iEndIndx, jStartIndx, jEndIndx;
        FloatMatrix Bi, Bj, D, Dij, DBj;
        for ( int i = 0; i < elemGeometry->giveNumberOfIntegrationRules(); i++ ) {
			IntegrationRule* iRi;
			iRi = elemGeometry->giveIntegrationRule( i );
            iStartIndx = iRi->getStartIndexOfLocalStrainWhereApply();
            iEndIndx   = iRi->getEndIndexOfLocalStrainWhereApply();
			for ( int j = 0; j < elemGeometry->giveNumberOfIntegrationRules(); j++ ) {
                IntegrationRule *iRule, *iRj;
				iRj = elemGeometry->giveIntegrationRule( j );
                jStartIndx = iRj->getStartIndexOfLocalStrainWhereApply();
                jEndIndx   = iRj->getEndIndexOfLocalStrainWhereApply();
                if ( i == j ) {
                    iRule = iRi;
                } else if ( iRi->giveNumberOfIntegrationPoints() < iRj->giveNumberOfIntegrationPoints() ) {
                    iRule = iRi;
                } else {
                    iRule = iRj;
                }

                for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                    GaussPoint *gp = iRule->getIntegrationPoint(k);

                    // Engineering (small strain) stiffness dSig/dEps
                    if ( nlGeometry == 0 ) {
                        this->computeBmatrixAt(gp, Bi, elemGeometry, iStartIndx, iEndIndx);
						this->computeConstitutiveMatrixAt(D, rMode, gp, tStep, elemGeometry);
                    } else if ( nlGeometry == 1 ) {
						if ( elemGeometry->giveDomain()->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                            this->computeBmatrixAt(gp, Bi, elemGeometry);
                            cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                        } else { // Material stiffness dP/dF
                            this->computeBHmatrixAt(gp, Bi, elemGeometry);
                            cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                        }
                    }


                    if ( i != j ) {
                        if ( nlGeometry == 0 ) {
                            this->computeBmatrixAt(gp, Bj, elemGeometry, jStartIndx, jEndIndx);
                        } else if ( nlGeometry == 1 ) {
                            if ( elemGeometry->giveDomain()->giveEngngModel()->giveFormulation() == AL ) {
                                this->computeBmatrixAt(gp, Bj, elemGeometry);
                            } else {
                                this->computeBHmatrixAt(gp, Bj, elemGeometry);
                            }
                        }
                    } else {
                        Bj  = Bi;
                    }

                    Dij.beSubMatrixOf(D, iStartIndx, iEndIndx, jStartIndx, jEndIndx);
                    double dV = elemGeometry->computeVolumeAround(gp);
                    DBj.beProductOf(Dij, Bj);
                    if ( matStiffSymmFlag ) {
                        answer.plusProductSymmUpper(Bi, DBj, dV);
                    } else {
                        answer.plusProductUnsym(Bi, DBj, dV);
                    }
                }
            }
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}

void StructuralElementEvaluator2 :: computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer,
                                                                      MatResponseMode rMode, TimeStep *tStep, ElementGeometry* elemGeometry)
{
	StructuralCrossSection *cs = this->giveStructuralCrossSection(elemGeometry);
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.resize(0, 0);
    if ( !elemGeometry->isActivated(tStep) ) {
        return;
    }

    FloatMatrix temp;
    FloatMatrix *m = & answer;
	if (elemGeometry->giveInterpolation() && elemGeometry->giveInterpolation()->hasSubPatchFormulation()) {
        m = & temp;
    }

    // Compute matrix from material stiffness
    FloatMatrix B, D, DB;
    IntArray irlocnum;
    for ( int ir = 0; ir < elemGeometry->giveNumberOfIntegrationRules(); ir++ ) {
        IntegrationRule *iRule = elemGeometry->giveIntegrationRule( ir );
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(j);

            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B, elemGeometry);
                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep, elemGeometry);
			// geometrical nonlinearity
            } else if ( nlGeometry == 1 ) {
				if ( elemGeometry->giveDomain()->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                    this->computeBmatrixAt(gp, B, elemGeometry);
                    cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                } else { // Material stiffness dP/dF
                    this->computeBHmatrixAt(gp, B, elemGeometry);
                    cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                }
            }

            double dV = elemGeometry->computeVolumeAround(gp);
            DB.beProductOf(D, B);
            if ( matStiffSymmFlag ) {
                m->plusProductSymmUpper(B, DB, dV);
            } else {
                m->plusProductUnsym(B, DB, dV);
            }
        }

        // localize irule contribution into element matrix
        if ( elemGeometry->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
            answer.assemble(* m, irlocnum);
            m->resize(0, 0);
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
StructuralElementEvaluator2 :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemGeometry)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;


    if ( !elemGeometry->isActivated(tStep) ) {
        answer.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        answer.zero();
        return;
    }

    this->computeBmatrixAt(gp, b, elemGeometry);
	this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u, elemGeometry);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }
    answer.beProductOf(b, u);
}


void
	StructuralElementEvaluator2 :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemGeometry)
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step tStep. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{
    this->giveStructuralCrossSection(elemGeometry)->giveRealStresses(answer, gp, strain, tStep);
}


void
StructuralElementEvaluator2 :: giveInternalForcesVector(FloatArray &answer,
                                              TimeStep *tStep, ElementGeometry* elemGeometry, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// or this->computeFirstPKStress for geometrically nonlinear analysis
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    FloatMatrix B;
    FloatArray vStress, vStrain, u;
	
    // This function can be quite costly to do inside the loops when one has many slave dofs.
	this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u, elemGeometry);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

	IntegrationRule *iRule = elemGeometry->giveDefaultIntegrationRulePtr();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

        // Engineering (small strain) stress
        if ( nlGeometry == 0 ) {
            this->computeBmatrixAt(gp, B, elemGeometry);
            if ( useUpdatedGpRecord == 1 ) {
                vStress = matStat->giveStressVector();
            } else {
                ///@todo Is this really what we should do for inactive elements?
                if ( !elemGeometry->isActivated(tStep) ) {
                    vStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
                    vStrain.zero();
                }
                vStrain.beProductOf(B, u);
                this->computeStressVector(vStress, vStrain, gp, tStep,elemGeometry);
            }
		// nonlinear geometry
        } else if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
            if ( elemGeometry->giveDomain()->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->giveCVector();
                } else {
                    this->computeCauchyStressVector(vStress, gp, tStep, elemGeometry);
                }

                this->computeBmatrixAt(gp, B, elemGeometry);
            } else { // First Piola-Kirchhoff stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->givePVector();
                } else {
                    this->computeFirstPKStressVector(vStress, gp, tStep, elemGeometry);
                    ///@todo This is actaully inefficient since it constructs B and twice and collects the nodal unknowns over and over.
                }

                this->computeBHmatrixAt(gp, B, elemGeometry);
            }
        }

        if ( vStress.giveSize() == 0 ) { /// @todo is this check really necessary?
            break;
        }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = elemGeometry->computeVolumeAround(gp);

		// nonlinear geometry
        if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
			if(vStress.giveSize() == 9) {
				FloatArray stressTemp;
				StructuralMaterial::giveReducedVectorForm(stressTemp, vStress, gp->giveMaterialMode());
				answer.plusProduct(B, stressTemp, dV);
			}
			else {
				answer.plusProduct(B, vStress, dV);
			}
        }
        else {
			if(vStress.giveSize() == 6) {
				// It may happen that e.g. plane strain is computed
				// using the default 3D implementation. If so,
				// the stress needs to be reduced.
				// (Note that no reduction will take place if
				//  the simulation is actually 3D.)
				FloatArray stressTemp;
				StructuralMaterial::giveReducedSymVectorForm(stressTemp, vStress, gp->giveMaterialMode());
				answer.plusProduct(B, stressTemp, dV);
			}
			else {
				answer.plusProduct(B, vStress, dV);
			}
        }

    }

    // If inactive: update fields but do not give any contribution to the internal forces
    if ( !elemGeometry->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}



void
StructuralElementEvaluator2 :: giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer,
                                                                   TimeStep *tStep, ElementGeometry* elemGeometry, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    /*
     * Returns nodal representation of real internal forces computed from first Piola-Kirchoff stress
     * if useGpRecord == 1 then stresses stored in the gp are used, otherwise stresses are computed
     * this must be done if you want internal forces after element->updateYourself() has been called
     * for the same time step.
     * The integration procedure uses an integrationRulesArray for numerical integration.
     * Each integration rule is considered to represent a separate sub-cell/element. Typically this would be used when
     * integration of the element domain needs special treatment, e.g. when using the XFEM.
     */

    FloatMatrix B;
    FloatArray vStress, vStrain;

    IntArray irlocnum;
    FloatArray *m = & answer, temp;
	if (elemGeometry->giveInterpolation() && elemGeometry->giveInterpolation()->hasSubPatchFormulation()) {
        m = & temp;
    }

    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);


    // loop over individual integration rules
    for ( int ir = 0; ir < elemGeometry->giveNumberOfIntegrationRules(); ir++ ) {
		IntegrationRule *iRule = elemGeometry->giveIntegrationRule( ir );

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B, elemGeometry);
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->giveStressVector();
                } else {
                    this->computeStrainVector(vStrain, gp, tStep, elemGeometry);
                    this->computeStressVector(vStress, vStrain, gp, tStep,elemGeometry);
                }
            } else if ( nlGeometry == 1 ) {
				if ( elemGeometry->giveDomain()->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = matStat->giveCVector();
                    } else {
                        this->computeCauchyStressVector(vStress, gp, tStep,elemGeometry);
                    }

                    this->computeBmatrixAt(gp, B, elemGeometry);
                } else { // First Piola-Kirchhoff stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = matStat->givePVector();
                    } else {
                        this->computeFirstPKStressVector(vStress, gp, tStep, elemGeometry);
                    }

                    this->computeBHmatrixAt(gp, B, elemGeometry);
                }
            }

            if ( vStress.giveSize() == 0 ) { //@todo is this really necessary?
                break;
            }

            // compute nodal representation of internal forces at nodes as f = B^T*stress dV
            double dV = elemGeometry->computeVolumeAround(gp);
            m->plusProduct(B, vStress, dV);

            // localize irule contribution into element matrix
            if ( elemGeometry->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
                answer.assemble(* m, irlocnum);
                m->resize(0);
            }
        }
    }

    // if inactive: update fields but do not give any contribution to the structure
    if ( !elemGeometry->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


void
StructuralElementEvaluator2 :: giveCharacteristicMatrix(FloatMatrix &answer,
                                              CharType mtrx, TimeStep *tStep, ElementGeometry* elemGeometry)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep, elemGeometry);
    } else if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep, elemGeometry);
    } else if ( mtrx == SecantStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, SecantStiffness, tStep, elemGeometry);
    } else if ( mtrx == ElasticStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, ElasticStiffness, tStep, elemGeometry);
    } else if ( mtrx == MassMatrix ) {
        this->computeMassMatrix(answer, tStep, elemGeometry);
    } else if ( mtrx == LumpedMassMatrix ) {
        this->computeLumpedMassMatrix(answer, tStep, elemGeometry);
    } else if ( mtrx == InitialStressMatrix ) {
        this->computeInitialStressMatrix(answer, tStep, elemGeometry);
    } else {
       // _error2( "giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}



void
StructuralElementEvaluator2 :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
TimeStep *tStep, ElementGeometry* elemGeometry)
//
// returns characteristics vector of receiver according to mtrx
//
{
    if ( mtrx == ExternalForcesVector ) {
        this->computeForceLoadVector(answer, tStep, mode, elemGeometry);
    } else if ( ( mtrx == InternalForcesVector ) && ( mode == VM_Total ) ) {
        this->giveInternalForcesVector(answer, tStep, elemGeometry);
    } else if ( ( mtrx == LastEquilibratedInternalForcesVector ) && ( mode == VM_Total ) ) {
        /* here tstep is not relevant, we set useUpdatedGpRecord = 1
         * and this will cause to integrate internal forces using existing (nontemp, equlibrated) stresses in
         * statuses. Mainly used to compute reaction forces */
        this->giveInternalForcesVector(answer, tStep, elemGeometry, 1);
    } else {
        //_error2( "giveCharacteristicVector: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}

// ** methods related to geometrical nonlinearity
void
StructuralElementEvaluator2 :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemGeometry)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
	
	this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u, elemGeometry); // solution vector
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt(gp, B, elemGeometry);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
        answer.at(3) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at(1) += 1.0;
    } else {
        OOFEM_SIMPLE_ERROR( "computeDeformationGradientVector : MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}


void
	StructuralElementEvaluator2 :: computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemGeometry)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection(elemGeometry);

    FloatArray vF;
    this->computeDeformationGradientVector(vF, gp, tStep, elemGeometry);
    cs->giveFirstPKStresses(answer, gp, vF, tStep);
}

void
StructuralElementEvaluator2 :: computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemGeometry)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection(elemGeometry);

    FloatArray vF;
    this->computeDeformationGradientVector(vF, gp, tStep, elemGeometry);
    cs->giveCauchyStresses(answer, gp, vF, tStep);
}

// end of methods related to geometrical nonlinearity ** 



void
StructuralElementEvaluator2 :: updateYourself(TimeStep *tStep, ElementGeometry* elemGeometry)
{
//    elemGeometry->updateYourself(tStep);

    // record initial displacement if element not active
    if ( !elemGeometry->isActivated(tStep) ) {
        if ( !initialDisplacements ) {
            initialDisplacements = new FloatArray();
        }

		this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, *initialDisplacements, elemGeometry);
    }
	
}


void
StructuralElementEvaluator2 :: updateInternalState(TimeStep *tStep, ElementGeometry* elemGeometry)
// Updates the receiver at end of step.
{
    FloatArray stress, strain;
	
    // force updating strains & stresses
    for ( int i = 0; i < elemGeometry->giveNumberOfIntegrationRules(); i++ ) {
        IntegrationRule *iRule = elemGeometry->giveIntegrationRule( i );
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(j);
            this->computeStrainVector(strain, gp, tStep,elemGeometry);
            this->computeStressVector(stress, strain, gp, tStep,elemGeometry);
        }
    }
}

void
StructuralElementEvaluator2 :: updateBeforeNonlocalAverage(TimeStep *tStep, ElementGeometry* elemGeometry)
// Updates the local material quantities before nonlocal averaging
{
    /*
     * Nonlocal material, related to receiver is also passed, because it is not possible to
     * ask receiver  for its material model  pointer using giveMaterial() service (returns Material type)
     * and then to cast this Material pointer into NonlocalMaterial pointer type,
     * because it is possible to cast only the pointer of derived class to pointer to base class.
     */
    IntegrationRule *iRule;
    FloatArray epsilon;
	

#ifdef __PARALLEL_MODE
    if ( this->giveParallelMode() == Element_remote ) {
        return;
    }

#endif

    // force updating local quantities
    for ( int i = 0; i < elemGeometry->giveNumberOfIntegrationRules(); i++ ) {
		iRule = elemGeometry->giveIntegrationRule( i );
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(j);
			this->computeStrainVector(epsilon, ip, tStep, elemGeometry);
            // provide material local strain increment - as is provided to computeRealStresVector
            // allows to update internal vars to be averaged to new state

            // not possible - produces wrong result
            StructuralNonlocalMaterialExtensionInterface *materialExt;
            materialExt =  static_cast< StructuralNonlocalMaterialExtensionInterface * >( this->giveStructuralCrossSection(elemGeometry)->
                                                                                          giveMaterialInterface(NonlocalMaterialExtensionInterfaceType, ip) );

            if ( !materialExt ) {
                return;             //_error ("updateBeforeNonlocalAverage: material with no StructuralNonlocalMaterial support");
            }
            materialExt->updateBeforeNonlocAverage(epsilon, iRule->getIntegrationPoint(j), tStep);
        }
    }
}


int
StructuralElementEvaluator2 :: checkConsistency(ElementGeometry* elemGeometry)
//
// check internal consistency
// mainly tests, whether crossSection data
// are safe for conversion to "Structural" version
//
{
	
    int result = 1;
    if ( !elemGeometry->giveCrossSection()->testCrossSectionExtension(CS_StructuralCapability) ) {
		//@todo correct this
		//_warning2( "checkConsistency : cross-section %s without structural support", this->giveCrossSection()->giveClassName() );
        result = 0;
    }


	if ( this->nlGeometry != 0  &&  this->nlGeometry != 1 ) {
        OOFEM_SIMPLE_ERROR("StructuralElementEvaluator :: checkConsistency - nlGeometry must be either 0 or 1 (%d not supported)", this->nlGeometry);
        return 0;
    } else {
        return 1;
    }


    return result;
	
}

void
	StructuralElementEvaluator2 :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer, ElementGeometry* elemGeometry)
{
	const int numNodes = elemGeometry->giveNumberOfDofManagers();
    FloatArray N(numNodes);

    const int dim = elemGeometry->giveSpatialDimension();

    answer.resize(dim, dim * numNodes);
    answer.zero();
    elemGeometry->giveInterpolation()->evalN( N, iLocCoord, FEIElementGeometryWrapper(elemGeometry) );

    answer.beNMatrixOf(N, dim);
	
}

void
StructuralElementEvaluator2 :: condense(FloatMatrix *stiff, FloatMatrix *mass, FloatArray *load, IntArray *what)
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
    double coeff, dii, lii = 0;
    FloatArray *gaussCoeff = NULL;
    // stored gauss coefficient
    // for mass condensation

    // check
    if ( !stiff->isSquare() ) {
        //_error("condense: stiffness size mismatch");
    }

    if ( mass ) {
        if ( !( mass->isSquare() && mass->giveNumberOfRows() == size ) ) {
           // _error("condense: mass size mismatch");
        }
    }

    if ( load ) {
        if ( !( load->giveSize() == size ) ) {
           // _error("condense: load size mismatch");
        }
    }

    // create gauss coeff array if mass condensation requested
    if ( mass ) {
        gaussCoeff = new FloatArray(size);
    }

    for ( i = 1; i <= nkon; i++ ) {
        ii  = what->at(i);
        if ( ( ii > size ) || ( ii <= 0 ) ) {
            //_error("condense: wrong dof number");
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
	StructuralElementEvaluator2 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep, ElementGeometry* elemGeometry)
{
    if ( type == IST_DisplacementVector ) {
        FloatArray u;
        FloatMatrix N;
		this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u, elemGeometry);
        this->computeNmatrixAt(* ( gp->giveLocalCoordinates() ), N, elemGeometry);
        answer.beProductOf(N, u);
        return 1;
    }
    return elemGeometry->giveIPValue(answer, gp, type, tStep);
	
}


void
	StructuralElementEvaluator2 :: giveNonlocalLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s, ElementGeometry* elemGeometry)
{
   
	NonlocalMaterialStiffnessInterface *interface;	
    IntArray elemLocArry;
    // create lit of remote elements, contributing to receiver
    std :: list< localIntegrationRecord > *integrationDomainList;
	IntegrationRule *iRule = elemGeometry->giveDefaultIntegrationRulePtr();
    std :: list< localIntegrationRecord > :: iterator pos;

    locationArray.resize(0);
    // loop over element IP
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        interface =  static_cast< NonlocalMaterialStiffnessInterface * >( this->giveStructuralCrossSection(elemGeometry)->
                                                                          giveMaterialInterface(NonlocalMaterialStiffnessInterfaceType, ip) );


        if ( interface == NULL ) {
            locationArray.resize(0);
            return;
        }

        integrationDomainList = interface->
                                NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(ip);
        // loop over IP influencing IPs, extract corresponding element numbers and their code numbers
        for ( pos = integrationDomainList->begin(); pos != integrationDomainList->end(); ++pos ) {
			pos->nearGp->giveElementGeometry()->giveDomain()->giveElementEvaluator(elemGeometry->giveNumber())->giveLocationArray(elemLocArry, EID_MomentumBalance, s, pos->nearGp->giveElementGeometry());
			
			/*
             * Currently no care given to multiple occurences of code number in locationArray.
             */
            locationArray.followedBy(elemLocArry, 20);
        }
    } // end loop over IPs
	
}


void
	StructuralElementEvaluator2 :: addNonlocalStiffnessContributions(SparseMtrx &dest, const UnknownNumberingScheme &s, TimeStep *tStep, ElementGeometry* elemGeometry)
{
    ///@todo Take into account cross section model (slaves)
    NonlocalMaterialStiffnessInterface *interface;
    if ( !elemGeometry->isActivated(tStep) ) {
        return;
    }

	IntegrationRule *iRule = elemGeometry->giveDefaultIntegrationRulePtr();
    // loop over element IP
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        interface = static_cast< NonlocalMaterialStiffnessInterface * >( this->giveStructuralCrossSection(elemGeometry)->
                                                                         giveMaterialInterface(NonlocalMaterialStiffnessInterfaceType, ip) );
        if ( interface == NULL ) {
            return;
        }

        interface->NonlocalMaterialStiffnessInterface_addIPContribution(dest, s, iRule->getIntegrationPoint(i), tStep);
    }
	
}


int
	StructuralElementEvaluator2 :: adaptiveUpdate(TimeStep *tStep, ElementGeometry* elemGeometry)
{
	
    int result = 1;
    IntegrationRule *iRule;
    FloatArray strain;
    MaterialModelMapperInterface *interface;
    for ( int i = 0; i < elemGeometry->giveNumberOfIntegrationRules(); i++ ) {
        iRule = elemGeometry->giveIntegrationRule( i );
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(i);
            interface = static_cast< MaterialModelMapperInterface * >( this->giveStructuralCrossSection(elemGeometry)->
                                                                       giveMaterialInterface(MaterialModelMapperInterfaceType, ip) );

            if ( interface == NULL ) {
                return 0;
            }
            this->computeStrainVector(strain, ip, tStep, elemGeometry);
            result &= interface->MMI_update(ip, tStep, & strain);
        }
    }

    return result;
	
}

IRResultType
	StructuralElementEvaluator2 :: initializeFrom(InputRecord *ir)
{
    /*ElementGeometry :: initializeFrom(ir);
	nlGeometry = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, _IFT_StructuralElementEvaluator2_nlgeoflag);

    return IRRT_OK;
	*/
	return IRRT_OK;
}

void StructuralElementEvaluator2 :: giveInputRecord(DynamicInputRecord &input)
{
    //ElementGeometry :: giveInputRecord(input);
	//    input.setField(nlGeometry, _IFT_NLStructuralElement_nlgeoflag);

    /// TODO: Should initialDisplacements be stored? /ES
}


StructuralCrossSection *StructuralElementEvaluator2 :: giveStructuralCrossSection(ElementGeometry* elemGeometry)
{
    return static_cast< StructuralCrossSection * >( elemGeometry->giveCrossSection() );
}

void StructuralElementEvaluator2 :: createMaterialStatus()
{
  /*  StructuralCrossSection *cs = giveStructuralCrossSection();
    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        IntegrationRule *iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint &gp = * ( iRule->getIntegrationPoint(j) );
            cs->createMaterialStatus(gp);
        }
    }
	*/
}


#ifdef __OOFEG

//
int
StructuralElementEvaluator2 :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode, int node, TimeStep *tStep)
{
    if ( type == IST_DisplacementVector ) {
        Node *n = this->giveNode(node);
        answer.resize(3);
        answer.at(1) = n->giveUpdatedCoordinate(1, tStep) - n->giveCoordinate(1);
        answer.at(2) = n->giveUpdatedCoordinate(2, tStep) - n->giveCoordinate(2);
        answer.at(3) = n->giveUpdatedCoordinate(3, tStep) - n->giveCoordinate(3);
        return 1;
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, tStep);
    }
}



void
StructuralElementEvaluator2 :: showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep)
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
StructuralElementEvaluator2 :: showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep)
{
    NonlocalMaterialStiffnessInterface *interface;
    int i;
    if ( ( ( mtrx == StiffnessMatrix ) || ( mtrx == TangentStiffnessMatrix ) ) ) {
        //interface = static_cast< NonlocalMaterialStiffnessInterface * >
        //            ( this->giveMaterial()->giveInterface(NonlocalMaterialStiffnessInterfaceType) );
        //if ( interface == NULL ) {
        //    return;
        //}

        IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        // loop over element IP
        for ( i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(i);
            interface = static_cast< NonlocalMaterialStiffnessInterface * >( this->giveStructuralCrossSection()->
                                                                             giveMaterialInterface(NonlocalMaterialStiffnessInterfaceType, ip) );

            if ( interface == NULL ) {
                return;
            }
            interface->NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(iRule->getIntegrationPoint(i), gc, tStep);
        }
    }
}

#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void OneDimensionalStructuralElementEvaluator :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry,
																  int lI, int uI)
{
	int nDofMan;
	nDofMan  = elemGeometry->giveNumberOfDofManagers();
	answer.resize(1, nDofMan);
    answer.zero();
    FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx( answer, * gp->giveCoordinates(), FEIElementGeometryWrapper(elemGeometry) );
}

void
OneDimensionalStructuralElementEvaluator :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry)
//
// Returns the one dimensional displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
	this->computeBmatrixAt(gp, answer, elemGeometry);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
void PlaneStressElementEvaluator2 :: computeBmatrixAt( GaussPoint *gp,FloatMatrix &answer, ElementGeometry* elemGeometry, int lI, int uI) 
{   
	int nDofMan;
    FloatMatrix d;
    nDofMan  = elemGeometry->giveNumberOfDofManagers();
	FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));

    answer.resize(3, nDofMan);
    answer.zero();

    for ( int i = 1; i <= nDofMan/2.; i++ )
      {
	answer.at(1, i * 2 - 1) = d.at(i, 1);
	answer.at(2, i * 2 - 0) = d.at(i, 2);
	
	answer.at(3, 2 * i - 1) = d.at(i, 2);
	answer.at(3, 2 * i - 0) = d.at(i, 1);
      }
  }

void PlaneStressElementEvaluator2 :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry)
// Returns the displacement gradient matrix {BH} of the receiver,
// evaluated at aGaussPoint.
// BH matrix  -  9 rows : du/dx, dv/dy, du/dy, dv/dx
// @todo not checked if correct
{
	int nDofMan;
    FloatMatrix d;
    nDofMan  = elemGeometry->giveNumberOfDofManagers();
    FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));
       
	answer.resize(4, nDofMan);
    answer.zero();
    // 3rd row is zero -> dw/dz = 0
    for ( int i = 1; i <= nDofMan/2; i++ ) {
        answer.at(1, 2 * i - 1) = d.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = d.at(i, 2);     // dv/dy -2
        answer.at(3, 2 * i - 1) = d.at(i, 2);     // du/dy -6
        answer.at(4, 2 * i - 0) = d.at(i, 1);     // dv/dx -9
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PlaneStrainElementEvaluator ::  computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry,
													  int lI, int uI)
{
	int nDofMan;
    FloatMatrix d;
    nDofMan  = elemGeometry->giveNumberOfDofManagers();
	FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));

    answer.resize(4, nDofMan);
    answer.zero();

    for ( int i = 1; i <= nDofMan/2.; i++ ) {
        answer.at(1, 2 * i - 1) = d.at(i, 1);
        answer.at(2, 2 * i - 0) = d.at(i, 2);
        answer.at(4, 2 * i - 1) = d.at(i, 2);
        answer.at(4, 2 * i - 0) = d.at(i, 1);
    }
}

void PlaneStrainElementEvaluator :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry)
// Returns the displacement gradient matrix {BH} of the receiver,
// evaluated at aGaussPoint.
// BH matrix  -  9 rows : du/dx, dv/dy, du/dy, dv/dx
// @todo not checked if correct
{
	int nDofMan;
    FloatMatrix d;
    nDofMan  = elemGeometry->giveNumberOfDofManagers();
	FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));
       
	answer.resize(5, nDofMan);
    answer.zero();
    // 3rd row is zero -> dw/dz = 0
    for ( int i = 1; i <= nDofMan/2; i++ ) {
        answer.at(1, 2 * i - 1) = d.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = d.at(i, 2);     // dv/dy -2
        answer.at(4, 2 * i - 1) = d.at(i, 2);     // du/dy -6
        answer.at(5, 2 * i - 0) = d.at(i, 1);     // dv/dx -9
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AxisymmetricStructuralElementEvaluator :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry,
                                  int li, int ui)
{
	
    int nDofMan = elemGeometry->giveNumberOfDofManagers();
    double r = 0., x;
    FloatArray n(nDofMan/2);
    FloatMatrix d;

   
    FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
    interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));
    interp->evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(elemGeometry) );
    
    

    answer.resize(6, nDofMan);
    answer.zero();


	    for ( int i = 1; i <= nDofMan; i++ ) {
			x  = elemGeometry->giveNode(i)->giveCoordinate(1);
            r += x * n.at(i);
        }

        for ( int i = 1; i <= nDofMan/2; i++ ) {
            answer.at(1, 2 * i - 1) = d.at(i, 1);
            answer.at(2, 2 * i - 0) = d.at(i, 2);
			answer.at(3,2*i-1)      = n.at(i) / r;
			answer.at(6, 2 * i - 1) = d.at(i, 2);
            answer.at(6, 2 * i - 0) = d.at(i, 1);
		}
}

void AxisymmetricStructuralElementEvaluator :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry)
{

// Returns the [9x8] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz = u/r, 0, 0, du/dy,  0, 0, dv/dx
// @todo not checked if correct
{
	int nDofMan;
	FloatArray n;
    FloatMatrix d;
    nDofMan  = elemGeometry->giveNumberOfDofManagers();
	FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));
    
    answer.resize(9, nDofMan);
    answer.zero();

    double r = 0., x;
    for ( int i = 1; i <= nDofMan; i++ ) {
		x  = elemGeometry->giveNode(i)->giveCoordinate(1);
        r += x * n.at(i);
    }


    // mode is _3dMat !!!!!! answer.at(4,*), answer.at(5,*), answer.at(7,*), and answer.at(8,*) is zero
    for ( int i = 1; i <= nDofMan/2; i++ ) {
        answer.at(1, 3 * i - 2) = d.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = d.at(i, 2);     // dv/dy
        answer.at(6, 3 * i - 2) = d.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = d.at(i, 1);     // dv/dx
		// dw/dz
		answer.at(3,2*i-1) = n.at(i) / r;
    }
    
}


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ThreeDimensionalStructuralElementEvaluator :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry,
                                  int lI, int uI )
{
	int nDofMan;
    FloatMatrix d;
    nDofMan  = elemGeometry->giveNumberOfDofManagers();
	FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));

    answer.resize(6, nDofMan);
    answer.zero();

    for ( int i = 1; i <= nDofMan/3.; i++ ) {
        answer.at(1, 3 * i - 2) = d.at(i, 1);
        answer.at(2, 3 * i - 1) = d.at(i, 2);
        answer.at(3, 3 * i - 0) = d.at(i, 3);

        answer.at(5, 3 * i - 2) = answer.at(4, 3 * i - 1) = d.at(i, 3);
        answer.at(6, 3 * i - 2) = answer.at(4, 3 * i - 0) = d.at(i, 2);
        answer.at(6, 3 * i - 1) = answer.at(5, 3 * i - 0) = d.at(i, 1);
	}
}

void ThreeDimensionalStructuralElementEvaluator :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemGeometry)
{
// Returns the displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz, dv/dz, du/dz, du/dy, dw/dy, dw/dx, dv/dx
    int nDofMan;
    FloatMatrix d;
    nDofMan  = elemGeometry->giveNumberOfDofManagers();
	FEInterpolation *interp = elemGeometry->giveInterpolation(displacementInterpolationNumber);
	interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( elemGeometry));

    answer.resize(9, nDofMan);
    answer.zero();

    for ( int i = 1; i <= nDofMan/3.; i++ ) {
        answer.at(1, 3 * i - 2) = d.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = d.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = d.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = d.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = d.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = d.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = d.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = d.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = d.at(i, 1);     // dv/dx
    }

}






} // end namespace oofem
