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

#include "Elements/structuralelement.h"
#include "CrossSections/structuralcrosssection.h"
#include "Materials/structuralmaterial.h"
#include "Materials/structuralms.h"
#include "Loads/structtemperatureload.h"
#include "Materials/structuralnonlocalmaterialext.h"
#include "Loads/structeigenstrainload.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "nonlocalmaterialext.h"
#include "load.h"
#include "boundaryload.h"
#include "pointload.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "nonlocmatstiffinterface.h"
#include "mathfem.h"
#include "materialmapperinterface.h"
#include "unknownnumberingscheme.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif


namespace oofem {
StructuralElement :: StructuralElement(int n, Domain *aDomain) :
    Element(n, aDomain)
{
}


StructuralElement :: ~StructuralElement()
{
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
    this->giveStructuralCrossSection()->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
}


void StructuralElement :: computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }
    // Just a wrapper for the deadweight body load computations:
    PointLoad *p = dynamic_cast< PointLoad * >(load);
    if ( p ) {
        FloatArray lcoords;
        if ( this->computeLocalCoordinates(lcoords, p->giveCoordinates()) ) {
            this->computePointLoadVectorAt(answer, load, tStep, mode);
        }
    } else {
        ///@todo This assumption of dead-weight loads needs to be lifted. We can have other body loads, such as
        this->computeBodyLoadVectorAt(answer, load, tStep, mode);
    }
}

  void StructuralElement :: computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    answer.clear();
    if ( type != ExternalForcesVector ) {
        return;
    }

    FEInterpolation *fei = this->giveInterpolation();
    if ( !fei ) {
        OOFEM_ERROR("No interpolator available");
    }

    FloatArray n_vec;
    FloatMatrix n, T;
    FloatArray force, globalIPcoords;
    //int nsd = fei->giveNsd();

    std :: unique_ptr< IntegrationRule >iRule( fei->giveBoundarySurfaceIntegrationRule(load->giveApproxOrder(), boundary) );

    for ( GaussPoint *gp: *iRule ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValueAt(force, tStep, lcoords, mode);
        } else {
            fei->boundaryLocal2Global( globalIPcoords, boundary, lcoords, FEIElementGeometryWrapper(this) );
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
        this->computeSurfaceNMatrix(n, boundary, lcoords); // to allow adapttation on element level

        ///@todo Some way to ask for the thickness at a global coordinate maybe?
        double thickness = 1.0; // Should be the circumference for axisymm-elements.
        double dV = thickness * gp->giveWeight() * fei->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(this) );
        answer.plusProduct(n, force, dV);
    }
}


void
StructuralElement::computeSurfaceNMatrix (FloatMatrix &answer, int boundaryID, const FloatArray& lcoords)
{
  FloatArray n_vec;
  this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, FEIElementGeometryWrapper(this) );
  answer.beNMatrixOf(n_vec, this->giveInterpolation()->giveNsd());
}

  
void StructuralElement :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    answer.clear();
    if ( type != ExternalForcesVector ) {
        return;
    }

    FEInterpolation *fei = this->giveInterpolation();
    if ( !fei ) {
        OOFEM_ERROR("No interpolator available");
    }

    FloatMatrix n, T;
    FloatArray force, globalIPcoords;

    std :: unique_ptr< IntegrationRule >iRule( fei->giveBoundaryEdgeIntegrationRule(load->giveApproxOrder(), boundary) );

    for ( GaussPoint *gp: *iRule ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValueAt(force, tStep, lcoords, mode);
        } else {
            fei->boundaryEdgeLocal2Global( globalIPcoords, boundary, lcoords, FEIElementGeometryWrapper(this) );
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
            if ( this->computeLoadLEToLRotationMatrix(T, boundary, gp) ) {
               force.rotatedWith(T, 'n');
            }
            // then to global c.s
            if ( this->computeLoadGToLRotationMtrx(T) ) {
                force.rotatedWith(T, 't');
            }
        }

        // Construct n-matrix
        //fei->boundaryEdgeEvalN( n_vec, boundary, lcoords, FEIElementGeometryWrapper(this) );
	//n.beNMatrixOf(n_vec, nsd);
	this->computeEdgeNMatrix(n, boundary, lcoords); // to allow adapttation on element level

        double dV = gp->giveWeight() * fei->boundaryEdgeGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(this) );
        answer.plusProduct(n, force, dV);
    }
}

void
StructuralElement::computeEdgeNMatrix (FloatMatrix &answer, int boundaryID, const FloatArray& lcoords)
{
  FloatArray n_vec;
  this->giveInterpolation()->boundaryEdgeEvalN (n_vec, boundaryID, lcoords, FEIElementGeometryWrapper(this) );
  answer.beNMatrixOf(n_vec, this->giveInterpolation()->giveNsd());
}

				       

  
void
StructuralElement :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV;
    FloatArray force, ntf;
    FloatMatrix n, T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);
    // transform from global to element local c.s
    if ( this->computeLoadGToLRotationMtrx(T) ) {
        force.rotatedWith(T, 'n');
    }

    answer.clear();

    if ( force.giveSize() ) {
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            this->computeNmatrixAt(gp->giveSubPatchCoordinates(), n);
            dV  = this->computeVolumeAround(gp);
            dens = this->giveCrossSection()->give('d', gp);
            ntf.beTProductOf(n, force);
            answer.add(dV * dens, ntf);
        }
    } else {
        return;
    }
}

void
StructuralElement :: computePointLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode, bool global)
{
    FloatArray force, lcoords;
    FloatMatrix T, n;

    PointLoad *pointLoad = dynamic_cast< PointLoad * >(load);
    FloatArray coords = pointLoad->giveCoordinates();
    pointLoad->computeValueAt(force, tStep, coords, mode);
    if ( this->computeLocalCoordinates(lcoords, pointLoad->giveCoordinates()) ) {
        this->computeNmatrixAt(lcoords, n);
        answer.beTProductOf(n, force);
    } else {
        OOFEM_WARNING("point load outside element");
    }

    // transform force
    if ( pointLoad->giveCoordSystMode() == Load :: CST_Global ) {
        // transform from global to element local c.s
        if ( this->computeLoadGToLRotationMtrx(T) ) {
            answer.rotatedWith(T, 'n');
        }
    }
}

int 
StructuralElement :: giveNumberOfIPForMassMtrxIntegration()
{
    IntegrationRule *iRule = this->giveIntegrationRule(0);
    // returns necessary number of ip to integrate the mass matrix 
    // \int_V N^T*N dV => (order of the approximation)*2 (constant density assumed)
    ///TODO this is without the jacobian and density
    int order = this->giveInterpolation()->giveInterpolationOrder();
    return iRule->getRequiredNumberOfIntegrationPoints(this->giveInterpolation()->giveIntegrationDomain(), 2*order);
}

void
StructuralElement :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
// Computes numerically the consistent (full) mass matrix of the receiver.
{
    int nip, ndofs = computeNumberOfDofs();
    double density, dV;
    FloatMatrix n;
    GaussIntegrationRule iRule(1, this, 1, 1);
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    if ( ( nip = this->giveNumberOfIPForMassMtrxIntegration() ) == 0 ) {
        OOFEM_ERROR("no integration points available");
    }

    iRule.setUpIntegrationPoints( this->giveIntegrationDomain(),
                                 nip, this->giveMaterialMode() );

    this->giveMassMtrxIntegrationgMask(mask);

    mass = 0.;

    for ( GaussPoint *gp: iRule ) {
        this->computeNmatrixAt(gp->giveSubPatchCoordinates(), n);
        density = this->giveCrossSection()->give('d', gp);

        if ( ipDensity != NULL ) {
            // Override density if desired
            density = * ipDensity;
        }

        dV = this->computeVolumeAround(gp);
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
    double mass = 0.;

    IntArray nodeDofIDMask, dimFlag(3);
    int indx = 0, ldofs, dim;
    double summ;

    if ( !this->isActivated(tStep) ) {
        int ndofs = computeNumberOfDofs();
        answer.resize(ndofs, ndofs);
        answer.zero();
        return;
    }

    this->computeConsistentMassMatrix(answer, tStep, mass);
    ldofs = answer.giveNumberOfRows();

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManDofIDMask(i, nodeDofIDMask);
        //this->giveDofManager(i)->giveLocationArray(nodeDofIDMask, nodalArray);
        for ( int j = 1; j <= nodeDofIDMask.giveSize(); j++ ) {
            indx++;
            // zero all off-diagonal terms
            for (int k = 1; k <= ldofs; k++ ) {
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
        OOFEM_ERROR("internal consistency check failed");
    }

    dim = dimFlag.at(1) + dimFlag.at(2) + dimFlag.at(3);
    summ = 0.;
    for ( int k = 1; k <= ldofs; k++ ) {
        summ += answer.at(k, k);
    }

    answer.times(dim * mass / summ);
}


void
StructuralElement :: computeResultingIPTemperatureAt(FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode)
// Computes at tStep the resulting force due to all temperature loads that act
// on the receiver. This force is used by the element for computing its
// body load vector.
{
    int n, nLoads;
    Load *load;
    FloatArray gCoords, temperature;

    if ( this->computeGlobalCoordinates( gCoords, gp->giveNaturalCoordinates() ) == 0 ) {
        OOFEM_ERROR("computeGlobalCoordinates failed");
    }

    answer.clear();
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        n = bodyLoadArray.at(i);
        load = domain->giveLoad(n);
        if ( load->giveBCValType() == TemperatureBVT ) {
            load->computeValueAt(temperature, tStep, gCoords, mode);
            answer.add(temperature);
        }
    }
}

void
StructuralElement :: computeResultingIPEigenstrainAt(FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode)
// Computes at tStep all eigenstrains that act on the receiver. Eigenstrains are used by the element for computing its body load vector.
{
    int n, nLoads;
    Load *load;
    FloatArray gCoords, eigenstrain;

    if ( this->computeGlobalCoordinates( gCoords, gp->giveNaturalCoordinates() ) == 0 ) {
        OOFEM_ERROR("computeGlobalCoordinates failed");
    }

    answer.clear();
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        n = bodyLoadArray.at(i);
        load = domain->giveLoad(n);
        if ( load->giveBCValType() == EigenstrainBVT ) {
            load->computeValueAt(eigenstrain, tStep, gCoords, mode);
            answer.add(eigenstrain);
        }
    }
}


void
StructuralElement :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray u;
    FloatMatrix n;
    this->computeNmatrixAt(lcoords, n);
    this->computeVectorOf(mode, tStep, u);
    answer.beProductOf(n, u);
}

void
StructuralElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                            TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    int iStartIndx, iEndIndx, jStartIndx, jEndIndx;
    double dV;
    FloatMatrix d, bi, bj, dbj, dij;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();

    if ( !this->isActivated(tStep) ) {
        return;
    }

    if ( integrationRulesArray.size() > 1 ) {
        for ( int i = 0; i < (int)integrationRulesArray.size(); i++ ) {
            iStartIndx = integrationRulesArray [ i ]->getStartIndexOfLocalStrainWhereApply();
            iEndIndx   = integrationRulesArray [ i ]->getEndIndexOfLocalStrainWhereApply();
            for ( int j = 0; j < (int)integrationRulesArray.size(); j++ ) {
                IntegrationRule *iRule;
                jStartIndx = integrationRulesArray [ j ]->getStartIndexOfLocalStrainWhereApply();
                jEndIndx   = integrationRulesArray [ j ]->getEndIndexOfLocalStrainWhereApply();
                if ( i == j ) {
                    iRule = integrationRulesArray [ i ].get();
                } else if ( integrationRulesArray [ i ]->giveNumberOfIntegrationPoints() < integrationRulesArray [ j ]->giveNumberOfIntegrationPoints() ) {
                    iRule = integrationRulesArray [ i ].get();
                } else {
                    iRule = integrationRulesArray [ j ].get();
                }

                for ( GaussPoint *gp: *iRule ) {
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
    } else {
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
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
    FloatMatrix temp, bj, d, dbj;
    int ndofs = this->computeNumberOfDofs();
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    IntArray irlocnum;

    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // loop over individual integration rules
    for ( auto &iRule: integrationRulesArray ) {
        // loop over individual integration points
        for ( GaussPoint *gp: *iRule ) {
            double dV = this->computeVolumeAround(gp);
            this->computeBmatrixAt(gp, bj);
            this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

            dbj.beProductOf(d, bj);
            if ( matStiffSymmFlag ) {
                m->plusProductSymmUpper(bj, dbj, dV);
            } else {
                m->plusProductUnsym(bj, dbj, dV);
            }
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, *iRule) ) {
            answer.assemble(* m, irlocnum);
            m->clear();
        }
    } // end loop over irules

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
StructuralElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;

    if ( !this->isActivated(tStep) ) {
        answer.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        answer.zero();
        return;
    }

    this->computeBmatrixAt(gp, b);
    this->computeVectorOf(VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }
    answer.beProductOf(b, u);
}


void
StructuralElement :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStresses(answer, gp, strain, tStep);
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
    FloatMatrix b;
    FloatArray u, stress, strain;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
        StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
        this->computeBmatrixAt(gp, b);

        if ( useUpdatedGpRecord == 1 ) {
            stress = matStat->giveStressVector();
        } else {
            if ( !this->isActivated(tStep) ) {
                strain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
                strain.zero();
            }
            strain.beProductOf(b, u);
            this->computeStressVector(stress, strain, gp, tStep);
        }

        // updates gp stress and strain record  acording to current
        // increment of displacement
        if ( stress.giveSize() == 0 ) {
            break;
        }

        // now every gauss point has real stress vector
        // compute nodal representation of internal forces using f = B^T*Sigma dV
        double dV = this->computeVolumeAround(gp);
        if ( stress.giveSize() == 6 ) {
            // It may happen that e.g. plane strain is computed
            // using the default 3D implementation. If so,
            // the stress needs to be reduced.
            // (Note that no reduction will take place if
            //  the simulation is actually 3D.)
            FloatArray stressTemp;
            StructuralMaterial :: giveReducedSymVectorForm( stressTemp, stress, gp->giveMaterialMode() );
            answer.plusProduct(b, stressTemp, dV);
        } else   {
            answer.plusProduct(b, stress, dV);
        }
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
    FloatMatrix b;
    FloatArray temp, u, stress, strain;
    IntArray irlocnum;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    FloatArray *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // loop over individual integration rules
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
            this->computeBmatrixAt(gp, b);

            if ( useUpdatedGpRecord == 1 ) {
                stress = matStat->giveStressVector();
            } else {
                if ( !this->isActivated(tStep) ) {
                    strain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
                    strain.zero();
                }
                strain.beProductOf(b, u);
                this->computeStressVector(stress, strain, gp, tStep);
            }

            //
            // updates gp stress and strain record  acording to current
            // increment of displacement
            //
            if ( stress.giveSize() == 0 ) {
                break;
            }

            //
            // now every gauss point has real stress vector
            //
            // compute nodal representation of internal forces using f = B^T*Sigma dV
            //
            double dV = this->computeVolumeAround(gp);
            m->plusProduct(b, stress, dV);

            // localize irule contribution into element matrix
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, *iRule) ) {
                answer.assemble(* m, irlocnum);
                m->clear();
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
    if ( mtrx == TangentStiffnessMatrix ) {
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
        OOFEM_ERROR("Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
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
      //this->computeForceLoadVector(answer, tStep, mode); // bp: assembled by emodel 
      answer.resize(0);
    } else if ( ( mtrx == InternalForcesVector ) && ( mode == VM_Total ) ) {
        this->giveInternalForcesVector(answer, tStep);
    } else if ( ( mtrx == LastEquilibratedInternalForcesVector ) && ( mode == VM_Total ) ) {
        /* here tstep is not relevant, we set useUpdatedGpRecord = 1
         * and this will cause to integrate internal forces using existing (nontemp, equlibrated) stresses in
         * statuses. Mainly used to compute reaction forces */
        this->giveInternalForcesVector(answer, tStep, 1);
    } else if ( mtrx == LumpedMassMatrix ) {
        FloatMatrix M;
        this->computeLumpedMassMatrix(M, tStep);
        answer.resize(M.giveNumberOfColumns());
        for ( int i = 0; i < M.giveNumberOfColumns(); ++i ) {
            answer[i] = M(i,i);
        }
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}

void
StructuralElement :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);

    // record initial displacement if element not active
    if ( activityTimeFunction && !isActivated(tStep) ) {
        if ( !initialDisplacements ) {
            initialDisplacements.reset( new FloatArray() );
        }

        this->computeVectorOf(VM_Total, tStep, * initialDisplacements);
    }
}


void
StructuralElement :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    FloatArray stress, strain;

    // force updating strains & stresses
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            this->computeStrainVector(strain, gp, tStep);
            this->computeStressVector(stress, strain, gp, tStep);
        }
    }
}

void
StructuralElement :: updateBeforeNonlocalAverage(TimeStep *tStep)
// Updates the local material quantities before nonlocal averaging
{
    /*
     * Nonlocal material, related to receiver is also passed, because it is not possible to
     * ask receiver  for its material model  pointer using giveMaterial() service (returns Material type)
     * and then to cast this Material pointer into NonlocalMaterial pointer type,
     * because it is possible to cast only the pointer of derived class to pointer to base class.
     */
    FloatArray epsilon;

    if ( this->giveParallelMode() == Element_remote ) {
        return;
    }

    // force updating local quantities
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            this->computeStrainVector(epsilon, gp, tStep);
            // provide material local strain increment - as is provided to computeRealStresVector
            // allows to update internal vars to be averaged to new state

            // not possible - produces wrong result
            StructuralNonlocalMaterialExtensionInterface *materialExt;
            materialExt =  static_cast< StructuralNonlocalMaterialExtensionInterface * >( this->giveStructuralCrossSection()->
                                                                                         giveMaterialInterface(NonlocalMaterialExtensionInterfaceType, gp) );

            if ( !materialExt ) {
                return;             //_error("updateBeforeNonlocalAverage: material with no StructuralNonlocalMaterial support");
            }
            materialExt->updateBeforeNonlocAverage(epsilon, gp, tStep);
        }
    }
}


int
StructuralElement :: checkConsistency()
//
// check internal consistency
// mainly tests, whether crossSection data
// are safe for conversion to "Structural" version
//
{
    int result = 1;
    if ( !this->giveCrossSection()->testCrossSectionExtension(CS_StructuralCapability) ) {
        OOFEM_WARNING("cross-section %s without structural support", this->giveCrossSection()->giveClassName() );
        result = 0;
    }

    return result;
}

void
StructuralElement :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    int numNodes = this->giveNumberOfDofManagers();
    FloatArray N(numNodes);

    int dim = this->giveSpatialDimension();

    answer.resize(dim, dim * numNodes);
    answer.zero();
    giveInterpolation()->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.beNMatrixOf(N, dim);
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
    double coeff, dii, lii = 0;
    FloatArray gaussCoeff;
    // stored gauss coefficient
    // for mass condensation

    // check
    if ( !stiff->isSquare() ) {
        OOFEM_ERROR("stiffness size mismatch");
    }

    if ( mass ) {
        if ( !( mass->isSquare() && mass->giveNumberOfRows() == size ) ) {
            OOFEM_ERROR("mass size mismatch");
        }
    }

    if ( load ) {
        if ( !( load->giveSize() == size ) ) {
            OOFEM_ERROR("load size mismatch");
        }
    }

    // create gauss coeff array if mass condensation requested
    if ( mass ) {
        gaussCoeff.resize(size);
    }

    for ( i = 1; i <= nkon; i++ ) {
        ii  = what->at(i);
        if ( ( ii > size ) || ( ii <= 0 ) ) {
            OOFEM_ERROR("wrong dof number");
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
                gaussCoeff.at(j) = coeff;
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
                        mass->at(j, k) += mass->at(j, ii) * gaussCoeff.at(k) +
                        mass->at(ii, k) * gaussCoeff.at(j) +
                        mass->at(ii, ii) * gaussCoeff.at(j) * gaussCoeff.at(k);
                    }
                }
            }

            for ( k = 1; k <= size; k++ ) {
                mass->at(ii, k) = 0.;
                mass->at(k, ii) = 0.;
            }
        }
    }
}


int
StructuralElement :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DisplacementVector ) {
        FloatArray u;
        FloatMatrix N;
        this->computeVectorOf(VM_Total, tStep, u);
        this->computeNmatrixAt(gp->giveSubPatchCoordinates(), N);
        answer.beProductOf(N, u);
        return 1;
    }
    return Element :: giveIPValue(answer, gp, type, tStep);
}


void
StructuralElement :: giveNonlocalLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s)
{
    NonlocalMaterialStiffnessInterface *interface;

    IntArray elemLocArry;
    // create lit of remote elements, contributing to receiver
    std :: list< localIntegrationRecord > *integrationDomainList;

    locationArray.clear();
    // loop over element IP
    for ( IntegrationPoint *ip: *this->giveDefaultIntegrationRulePtr() ) {
        interface =  static_cast< NonlocalMaterialStiffnessInterface * >( this->giveStructuralCrossSection()->
                                                                         giveMaterialInterface(NonlocalMaterialStiffnessInterfaceType, ip) );


        if ( interface == NULL ) {
            locationArray.clear();
            return;
        }

        integrationDomainList = interface->
                                NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(ip);
        // loop over IP influencing IPs, extract corresponding element numbers and their code numbers
        for ( auto &lir: *integrationDomainList ) {
            lir.nearGp->giveElement()->giveLocationArray(elemLocArry, s);
            /*
             * Currently no care given to multiple occurences of code number in locationArray.
             */
            locationArray.followedBy(elemLocArry, 20);
        }
    } // end loop over IPs
}


void
StructuralElement :: addNonlocalStiffnessContributions(SparseMtrx &dest, const UnknownNumberingScheme &s, TimeStep *tStep)
{
    ///@todo Take into account cross section model (slaves)
    NonlocalMaterialStiffnessInterface *interface;

    if ( !this->isActivated(tStep) ) {
        return;
    }

    // loop over element IP
    for ( IntegrationPoint *ip: *this->giveDefaultIntegrationRulePtr() ) {
        interface = static_cast< NonlocalMaterialStiffnessInterface * >( this->giveStructuralCrossSection()->
                                                                        giveMaterialInterface(NonlocalMaterialStiffnessInterfaceType, ip) );
        if ( interface == NULL ) {
            return;
        }

        interface->NonlocalMaterialStiffnessInterface_addIPContribution(dest, s, ip, tStep);
    }
}


int
StructuralElement :: adaptiveUpdate(TimeStep *tStep)
{
    int result = 1;
    FloatArray strain;

    MaterialModelMapperInterface *interface;
    for ( auto &iRule: integrationRulesArray ) {
        for ( IntegrationPoint *ip: *iRule ) {
            interface = static_cast< MaterialModelMapperInterface * >( this->giveStructuralCrossSection()->
                                                                      giveMaterialInterface(MaterialModelMapperInterfaceType, ip) );

            if ( interface == NULL ) {
                return 0;
            }
            this->computeStrainVector(strain, ip, tStep);
            result &= interface->MMI_update(ip, tStep, & strain);
        }
    }

    return result;
}

IRResultType
StructuralElement :: initializeFrom(InputRecord *ir)
{
    return Element :: initializeFrom(ir);
}

void StructuralElement :: giveInputRecord(DynamicInputRecord &input)
{
    Element :: giveInputRecord(input);

    /// TODO: Should initialDisplacements be stored? /ES
}


StructuralCrossSection *StructuralElement :: giveStructuralCrossSection()
{
    return static_cast< StructuralCrossSection * >( this->giveCrossSection() );
}

void StructuralElement :: createMaterialStatus()
{
    StructuralCrossSection *cs = giveStructuralCrossSection();
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            cs->createMaterialStatus(*gp);
        }
    }
}


#ifdef __OOFEG

//
int
StructuralElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode, int node, TimeStep *tStep)
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
StructuralElement :: showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep)
{
    if ( mtrx == TangentStiffnessMatrix ||
        mtrx == SecantStiffnessMatrix || mtrx == ElasticStiffnessMatrix ) {
        int i, j, n;
        IntArray loc;
        this->giveLocationArray( loc, EModelDefaultEquationNumbering() );

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
StructuralElement :: showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep)
{
    NonlocalMaterialStiffnessInterface *interface;
    if ( mtrx == TangentStiffnessMatrix ) {
        // loop over element IP
        for ( IntegrationPoint *ip: *this->giveDefaultIntegrationRulePtr() ) {
            interface = static_cast< NonlocalMaterialStiffnessInterface * >( this->giveStructuralCrossSection()->
                                                                            giveMaterialInterface(NonlocalMaterialStiffnessInterfaceType, ip) );

            if ( interface == NULL ) {
                return;
            }
            interface->NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(ip, gc, tStep);
        }
    }
}

#endif
} // end namespace oofem
