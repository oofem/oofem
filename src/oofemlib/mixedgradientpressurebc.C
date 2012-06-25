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


#include "mixedgradientpressurebc.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "engngm.h"
#include "node.h"
#include "activedof.h"
#include "masterdof.h"
#include "usrdefsub.h" // For sparse matrix creation.
#include "sparsemtrxtype.h"
#include "line2boundaryelement.h"

#include "sparsemtrx.h"
#include "sparselinsystemnm.h"

namespace oofem {

MixedGradientPressureBC :: MixedGradientPressureBC(int n, Domain *d) : ActiveBoundaryCondition(n,d)
{
    // The unknown volumetric strain
    voldman = new Node(100, d);
    voldman->appendDof(new MasterDof(1, voldman, X_1));

    int nsd = d->giveNumberOfSpatialDimensions();
    int components = (nsd+1)*nsd/2;
    // The prescribed strains.
    devdman = new Node(200, d);
    for (int i = 0; i < components; i++) {
        // Just putting in X_i id-items since they don't matter.
        // These don't actually need to be active, they are masterdofs with prescribed values, its
        // easier to just have them here rather than trying to make another Dirichlet boundary condition.
        devdman->appendDof(new ActiveDof(i+1, devdman, this->giveNumber(), (DofIDItem)(X_1+i) ));
    }
}


MixedGradientPressureBC :: ~MixedGradientPressureBC()
{
    delete voldman;
    delete devdman;
}


Dof *MixedGradientPressureBC :: giveVolDof()
{
    return voldman->giveDof(1);
}


int MixedGradientPressureBC :: giveNumberOfInternalDofManagers()
{
    return 2;
}


DofManager *MixedGradientPressureBC :: giveInternalDofManager(int i)
{
    if (i == 1) {
        return this->voldman;
    } else {
        return this->devdman;
    }
}


int MixedGradientPressureBC :: giveNumberOfMasterDofs(ActiveDof *dof)
{
    if (this->isDevDof(dof))
        return 1;
    return devdman->giveNumberOfDofs() + 1;
}


Dof *MixedGradientPressureBC :: giveMasterDof(ActiveDof *dof, int mdof)
{
    if (this->isDevDof(dof))
        return NULL;
    if (mdof == 1) {
        return voldman->giveDof(1);
    } else {
        return devdman->giveDof(mdof-1);
    }
}


void MixedGradientPressureBC :: computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();

    FloatArray dx;
    dx.beDifferenceOf(*coords, this->centerCoord);

    int nsd = dx.giveSize(); // Number of spatial dimensions

    masterContribs.resize(devdman->giveNumberOfDofs()+1);

    if ( id == D_u || id == V_u ) {
        masterContribs.at(1) = dx.at(1)/3.0;      // d_vol
        if (nsd == 2) {
            masterContribs.at(2) = dx.at(1);      // d_dev_11
            masterContribs.at(3) = 0.0;           // d_dev_22
            masterContribs.at(4) = dx.at(2)/2.0;  // gamma_12
        } else if (nsd == 3) {
            masterContribs.at(2) = dx.at(1);      // d_dev_11
            masterContribs.at(3) = 0.0;           // d_dev_22
            masterContribs.at(3) = 0.0;           // d_dev_33
            masterContribs.at(4) = 0.0;           // gamma_23
            masterContribs.at(5) = dx.at(3)/2.0;  // gamma_13
            masterContribs.at(6) = dx.at(2)/2.0;  // gamma_12
        }
    } else if ( id == D_v || id == V_v ) {
        masterContribs.at(1) = dx.at(2)/3.0;      // d_vol
        if (nsd == 2) {
            masterContribs.at(2) = 0.0;           // d_dev_11
            masterContribs.at(3) = dx.at(2);      // d_dev_22
            masterContribs.at(4) = dx.at(1)/2.0;  // gamma_12
        } else if (nsd == 3) {
            masterContribs.at(2) = 0.0;           // d_dev_11
            masterContribs.at(3) = dx.at(2);      // d_dev_22
            masterContribs.at(4) = 0.0;           // d_dev_33
            masterContribs.at(5) = dx.at(3)/2.0;  // gamma_23
            masterContribs.at(6) = 0.0;           // gamma_13
            masterContribs.at(7) = dx.at(1)/2.0;  // gamma_12
        }
    } else if ( id == D_w || id == V_w ) { // 3D only:
        masterContribs.at(1) = dx.at(3)/3.0;  // d_vol
        masterContribs.at(2) = 0.0;           // d_dev_11
        masterContribs.at(3) = 0.0;           // d_dev_22
        masterContribs.at(3) = dx.at(3);      // d_dev_33
        masterContribs.at(4) = dx.at(2)/2.0;  // gamma_23
        masterContribs.at(4) = dx.at(1)/2.0;  // gamma_13
        masterContribs.at(4) = 0.0;           // gamma_12
    } else {
        OOFEM_ERROR("MixedGradientPressureBC :: computeDofTransformation - Incompatible id on subjected dof\n");
    }
}


double MixedGradientPressureBC :: giveUnknown(double vol, const FloatArray &dev, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();

    if ( coords == NULL || coords->giveSize() != this->centerCoord.giveSize() ) {
        OOFEM_ERROR("MixedGradientPressureBC :: give - Size of coordinate system different from center coordinate in b.c.");
    }

    FloatArray dx;
    dx.beDifferenceOf(*coords, this->centerCoord);

    int nsd = dx.giveSize(); // Number of spatial dimensions

    double dev11, dev22, dev33, gam23, gam13, gam12;
    if (nsd == 2) {
        dev11 = dev.at(1);
        dev22 = dev.at(2);
        gam12 = dev.at(3);
    } else if (nsd == 3) {
        dev11 = dev.at(1);
        dev22 = dev.at(2);
        dev33 = dev.at(3);
        gam23 = dev.at(4);
        gam13 = dev.at(5);
        gam12 = dev.at(6);
    }

    double val;
    if ( id == D_u || id == V_u ) {
        val = dx.at(1)/3.0*vol;
        if (nsd == 2) val += dx.at(1)*dev11 + dx.at(2)/2.0*gam12;
        if (nsd == 3) val += dx.at(3)/2.0*gam13;
    } else if ( id == D_v || id == V_v ) {
        val = dx.at(2)/3.0*vol + dx.at(1)/2.0*gam12 + dx.at(2)*dev22;
        if (nsd == 3) val += dx.at(3)/2.0*gam23;
    } else if ( id == D_w || id == V_w ) { // 3D only:
        val = dx.at(3)/3.0*vol + dx.at(1)/2.0*gam13 + dx.at(2)/2.0*gam23 + dx.at(3)*dev33;
    }
    return val;
}


double MixedGradientPressureBC :: domainSize()
{
    ///@todo This code is very general, and necessary for all multiscale simulations with pores, so it should be moved into Domain eventually
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    for (int i = 1; i <= this->domain->giveNumberOfElements(); ++i) {
        //BoundaryElement *e = dynamic_cast< BoundaryElement* >(d->giveElement(i));
        ///@todo Support more than 2D
        Line2BoundaryElement *e = dynamic_cast< Line2BoundaryElement* >(this->domain->giveElement(i));
        if (e) {
            domain_size += e->computeNXIntegral();
        }
    }
    if  (domain_size == 0.0) { // No boundary elements? Assume full density;
        return this->domain->giveArea(); ///@todo Support more than 2D
    } else {
        return domain_size/nsd;
    }
}


void MixedGradientPressureBC :: computeFields(FloatArray &sigmaDev, double &vol, EquationID eid, TimeStep *tStep)
{
    vol = this->giveVolDof()->giveUnknown(eid, VM_Total, tStep);
    int npeq = this->giveDomain()->giveEngngModel()->giveNumberOfPrescribedDomainEquations(this->giveDomain()->giveNumber(), eid);
    sigmaDev.resize(npeq);
    sigmaDev.zero();
    this->domain->giveEngngModel()->assembleVector(sigmaDev, tStep, eid,
        InternalForcesVector, VM_Total, EModelDefaultPrescribedEquationNumbering(), this->domain);
    // Divide by the RVE-volume
    sigmaDev.times(1.0/this->domainSize());
    // Above, full sigma = sigmaDev - p*I is assembled, so to obtain the deviatoric part we add p to the diagonal part.
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    for (int i = 1; i <= nsd; ++i ) {
        sigmaDev.at(i) += this->pressure;
    }
}


void MixedGradientPressureBC :: computeTangents(
    FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, EquationID eid, TimeStep *tStep)
{
    double size = this->domainSize();
    // Fetch some information from the engineering model
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    SparseLinearSystemNM *solver = CreateUsrDefSparseLinSolver(ST_Petsc, 1, this->domain, this->domain->giveEngngModel());// = rve->giveLinearSolver();
    SparseMtrx *Kff, *Kfp, *Kpf, *Kpp;
    SparseMtrxType stype = SMT_PetscMtrx;// = rve->giveSparseMatrixType();
    EModelDefaultEquationNumbering fnum;
    EModelDefaultPrescribedEquationNumbering pnum;

    // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
    Kff = CreateUsrDefSparseMtrx(stype);
    Kfp = CreateUsrDefSparseMtrx(stype);
    Kpf = CreateUsrDefSparseMtrx(stype);
    Kpp = CreateUsrDefSparseMtrx(stype);
    if ( !Kff ) {
        OOFEM_ERROR2("MixedGradientPressureBC :: computeTangents - Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure( rve, 1, eid, fnum, fnum );
    Kfp->buildInternalStructure( rve, 1, eid, fnum, pnum );
    Kpf->buildInternalStructure( rve, 1, eid, pnum, fnum );
    Kpp->buildInternalStructure( rve, 1, eid, pnum, pnum );
    rve->assemble(Kff, tStep, eid, StiffnessMatrix, fnum, fnum, this->domain );
    rve->assemble(Kfp, tStep, eid, StiffnessMatrix, fnum, pnum, this->domain );
    rve->assemble(Kpf, tStep, eid, StiffnessMatrix, pnum, fnum, this->domain );
    rve->assemble(Kpp, tStep, eid, StiffnessMatrix, pnum, pnum, this->domain );

    // Setup up indices and locations
    int neq = Kff->giveNumberOfRows();
    int npeq = Kpf->giveNumberOfRows();

    // Indices and such of internal dofs
    int dvol_eq = this->giveVolDof()->giveEqn();
    int ndev = this->devGradient.giveSize();
    // Matrices and arrays for sensitivities
    FloatMatrix ddev_pert(ndev,npeq); // In fact, npeq should most likely equal ndev
    FloatMatrix rhs_d(neq,npeq); // RHS for d_dev [d_dev11, d_dev22, d_dev12] in 2D
    FloatMatrix s_d(neq,npeq); // Sensitivity fields for d_dev
    FloatArray rhs_p(neq); // RHS for pressure
    FloatArray s_p(neq); // Sensitivity fields for p

    // Unit pertubations for d_dev
    ddev_pert.zero();
    for (int i = 1; i <= ndev; ++i) {
        int eqn = this->devdman->giveDof(i)->__givePrescribedEquationNumber();
        ddev_pert.at(eqn, i) = -1.0; // Minus sign for moving it to the RHS
    }
    Kfp->times(ddev_pert, rhs_d);

    // Sensitivity analysis for p (already in rhs, just set value directly)
    rhs_p.zero();
    rhs_p.at(dvol_eq) = -1.0; // dp = -1.0 (unit size)

    // Solve all sensitivities
    solver->solve(Kff,&rhs_p,&s_p);
    solver->solve(Kff,rhs_d,s_d);

    // Sensitivities for d_vol is solved for directly;
    Cp = s_p.at(dvol_eq);
    Cd.resize(ndev);
    for (int i = 1; i <= ndev; ++i) {
        Cd.at(i) = s_d.at(dvol_eq, i); // Copy over relevant row from solution
    }

    // Sensitivities for d_dev is obtained as reactions forces;
    Kpf->times(s_p, Ep);

    FloatMatrix tmpMat;
    Kpp->times(ddev_pert, tmpMat);
    Kpf->times(s_d, Ed);

    Ed.subtract(tmpMat);

    Ed.times(1.0/size);

    // Not sure if i actually need to do this part, but the obtained tangents are to dsigma/dp, not dsigma_dev/dp, so they need to be corrected;
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    for (int i = 1; i <= nsd; ++i) {
        Ep.at(i) += 1.0;
        Cd.at(i) += 1.0;
    }
    // Makes Ed 4th order deviatoric, in Voigt form (!)
    /*for (int i = 1; i <= Ed.giveNumberOfColumns(); ++i) {
        // Take the mean of the "diagonal" part of each column
        double mean = 0.0;
        for (int j = 1; j <= nsd; ++j) {
            mean += Ed.at(i,j);
        }
        mean /= (double)nsd;
        // And subtract it from that column
        for (int j = 1; j <= nsd; ++j) {
            Ed.at(i,j) -= mean;
        }
    }*/

    delete Kff;
    delete Kfp;
    delete Kpf;
    delete Kpp;
}


double MixedGradientPressureBC :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if (this->isDevDof(dof)) {
        return this->devGradient(dof->giveDofID() - X_1);
    }
    return this->giveUnknown(this->giveVolDof()->giveUnknown(field, mode, tStep), this->devGradient, mode, tStep, dof);
}


double MixedGradientPressureBC :: giveUnknown(EquationID eid, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if (this->isDevDof(dof)) {
        return this->devGradient(dof->giveDofID() - X_1);
    }
    return this->giveUnknown(this->giveVolDof()->giveUnknown(eid, mode, tStep), this->devGradient, mode, tStep, dof);
}


void MixedGradientPressureBC :: setPrescribedDeviatoricGradientFromVoigt(const FloatArray &t)
{
    devGradient = t;
}


double MixedGradientPressureBC :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                    CharType type, ValueModeType mode, const UnknownNumberingScheme &s, Domain *domain)
{
    if (type != ExternalForcesVector)
        return 0.0;

    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance)
        return 0.0;

    int vol_loc = this->giveVolDof()->giveEquationNumber(s);

    double rve_size = this->domainSize();
    answer.at(vol_loc) -= rve_size*pressure; // Note the negative sign (pressure as opposed to mean stress)

    return 0.0;
}


bool MixedGradientPressureBC :: isPrimaryDof(ActiveDof *dof)
{
    return this->isDevDof(dof);
}


double MixedGradientPressureBC :: giveBcValue(ActiveDof *dof, ValueModeType mode, TimeStep *tStep)
{
    if (this->isDevDof(dof)) {
        return this->devGradient(dof->giveDofID() - X_1);
    }
    OOFEM_ERROR("MixedGradientPressureBC :: giveBcValue - Has no prescribed value from bc.");
    return 0.0;
}


bool MixedGradientPressureBC :: hasBc(ActiveDof *dof, TimeStep *tStep)
{
    return this->isDevDof(dof);
}


bool MixedGradientPressureBC :: isDevDof(ActiveDof *dof)
{
    for (int i = 1; i <= this->devdman->giveNumberOfDofs(); ++i) {
        if (devdman->giveDof(i) == dof) {
            return true;
        }
    }
    return false;
}


IRResultType MixedGradientPressureBC :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    GeneralBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->devGradient, IFT_MixedGradientPressureBC_devGradient, "devgradient");
    IR_GIVE_FIELD(ir, this->pressure, IFT_MixedGradientPressureBC_pressure, "pressure");

    IRResultType rt = IR_GIVE_OPTIONAL_FIELD(ir, this->centerCoord, IFT_MixedGradientPressureBC_centerCoords, "ccoord")
    if ( rt != IRRT_OK ) {
        this->centerCoord.resize( domain->giveNumberOfSpatialDimensions() );
        this->centerCoord.zero();
    }

    return IRRT_OK;
}


int MixedGradientPressureBC :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    GeneralBoundaryCondition :: giveInputRecordString(str, keyword);

    sprintf( buff, " devgradient %d ", this->devGradient.giveSize() );
    for ( int i = 1; i <= this->devGradient.giveSize(); i++ ) {
        sprintf( buff, " %e", this->devGradient.at(i) );
        str += buff;
    }

    sprintf( buff, " ccoord %d", this->centerCoord.giveSize() );
    for ( int i = 1; i <= this->centerCoord.giveSize(); i++ ) {
        sprintf( buff, " %e", this->centerCoord.at(i) );
        str += buff;
    }

    return 1;
}
} // end namespace oofem

