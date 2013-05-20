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

#ifndef trplanrot_h
#define trplanrot_h

#include "trplanstrss.h"

///@name Input fields for TrPlaneStrRot
//@{
#define _IFT_TrPlaneStrRot_Name "trplanestrrot"
#define _IFT_TrPlaneStrRot_niprot "niprot"
//@}

namespace oofem {
/**
 * Class implements an triangular three-node  plane-
 * stress elasticity finite element with independent rotation field.
 * Each node has 3 degrees of freedom.
 */
class TrPlaneStrRot : public TrPlaneStress2d
{
protected:
    int numberOfRotGaussPoints;

public:
    TrPlaneStrRot(int, Domain *);
    virtual ~TrPlaneStrRot() { }

protected:
    virtual void computeGaussPoints();
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual double giveArea();
    virtual void giveNodeCoordinates(FloatArray &x, FloatArray &y);

    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode);

public:
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TrPlaneStrRot_Name; }
    virtual const char *giveClassName() const { return "TrPlaneStrRot"; }
    virtual classType giveClassID() const { return TrPlaneStrRotClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _PlaneStressRot; }
    virtual integrationDomain giveIntegrationDomain() { return _Triangle; }

    virtual int computeNumberOfDofs(EquationID ut) { return 9; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane) { return 0.; }

    FloatArray *GivePitch();
    FloatArray *GiveDerivativeUX(GaussPoint *gp);
    FloatArray *GiveDerivativeVX(GaussPoint *gp);
    FloatArray *GiveDerivativeUY(GaussPoint *gp);
    FloatArray *GiveDerivativeVY(GaussPoint *gp);
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);

};
} // end namespace oofem
#endif //  trplanrot_h
