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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
#include "classtype.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "loadtime.h"
#include "engngm.h"
#include "node.h"
#include "activedof.h"
#include "masterdof.h"

namespace oofem {

MixedGradientPressureBC :: MixedGradientPressureBC(int n, Domain *d) : ActiveBoundaryCondition(n,d)
{
    // The unknown volumetric strain
    voldman = new Node(100, d);
    voldman->appendDof(new MasterDof(1, voldman, X_1));

    int ndim = d->giveNumberOfSpatialDimensions();
    int components = (ndim+1)*ndim/2;
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

    int ndim = dx.giveSize();

    masterContribs.resize(devdman->giveNumberOfDofs()+1);

    if ( id == D_u || id == V_u ) {
        masterContribs.at(1) = dx.at(1)/ndim; // d_vol
        if (ndim == 2) {
            masterContribs.at(2) = dx.at(1);      // d_dev_11
            masterContribs.at(3) = 0.0;           // d_dev_22
            masterContribs.at(4) = dx.at(2)/2.0;  // gamma_12
        } else if (ndim == 3) {
            masterContribs.at(2) = dx.at(1);      // d_dev_11
            masterContribs.at(3) = 0.0;           // d_dev_22
            masterContribs.at(3) = 0.0;           // d_dev_33
            masterContribs.at(4) = 0.0;           // gamma_23
            masterContribs.at(5) = dx.at(3)/2.0;  // gamma_13
            masterContribs.at(6) = dx.at(2)/2.0;  // gamma_12
        }
    } else if ( id == D_v || id == V_v ) {
        masterContribs.at(1) = dx.at(2)/ndim; // d_vol
        if (ndim == 2) {
            masterContribs.at(2) = 0.0;           // d_dev_11
            masterContribs.at(3) = dx.at(2);      // d_dev_22
            masterContribs.at(4) = dx.at(1)/2.0;  // gamma_12
        } else if (ndim == 3) {
            masterContribs.at(2) = 0.0;           // d_dev_11
            masterContribs.at(3) = dx.at(2);      // d_dev_22
            masterContribs.at(4) = 0.0;           // d_dev_33
            masterContribs.at(5) = dx.at(3)/2.0;  // gamma_23
            masterContribs.at(6) = 0.0;           // gamma_13
            masterContribs.at(7) = dx.at(1)/2.0;  // gamma_12
        }
    } else if ( id == D_w || id == V_w ) { // 3D only:
        masterContribs.at(1) = dx.at(3)/ndim; // d_vol
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

    int ndim = dx.giveSize();

    double dev11, dev22, dev33, gam23, gam13, gam12;
    if (ndim == 2) {
        dev11 = dev.at(1);
        dev22 = dev.at(2);
        gam12 = dev.at(3);
    } else if (ndim == 3) {
        dev11 = dev.at(1);
        dev22 = dev.at(2);
        dev33 = dev.at(3);
        gam23 = dev.at(4);
        gam13 = dev.at(5);
        gam12 = dev.at(6);
    }

    double val;
    switch ( id ) {
        case D_u:
        case V_u:
            val = dx.at(1)/ndim*vol;
            if (ndim == 2) val += dx.at(1)*dev11 + dx.at(2)/2.0*gam12;
            if (ndim == 3) val += dx.at(3)/2.0*gam13;
        case D_v:
        case V_v:
            val = dx.at(2)/ndim*vol + dx.at(1)/2.0*gam12 + dx.at(2)*dev22;
            if (ndim == 3) val += dx.at(3)/2.0*gam23;
        case D_w:
        case V_w:
            val = dx.at(3)/ndim*vol + dx.at(1)/2.0*gam13 + dx.at(2)/2.0*gam23 + dx.at(3)*dev33;
        default:
            return 0.0;
    }
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


void MixedGradientPressureBC :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                    CharType type, ValueModeType mode, const UnknownNumberingScheme &s, Domain *domain)
{
    if (type != LoadVector)
        return;

    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance)
        return;

    int vol_loc = this->giveVolDof()->giveEquationNumber(s);

    double rve_size = this->domain->giveArea();
    answer.at(vol_loc) += rve_size*pressure;
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

