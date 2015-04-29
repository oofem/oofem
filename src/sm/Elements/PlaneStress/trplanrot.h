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

#ifndef trplanrot_h
#define trplanrot_h

#include "../sm/Elements/PlaneStress/trplanstrss.h"

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
 *
 * This element is based on the following paper:
 *   Ibrahimbegovic, A., Taylor, R.L., Wilson, E. L.: A robust quadrilateral membrane finite element with drilling degrees of freedom
 *   Int. J. Num. Meth. Engng., 30, 445-457, 1990.
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
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);

    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual double giveArea();
    virtual void giveNodeCoordinates(FloatArray &x, FloatArray &y);

    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode);

public:
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TrPlaneStrRot_Name; }
    virtual const char *giveClassName() const { return "TrPlaneStrRot"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _PlaneStressRot; }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int computeNumberOfDofs() { return 9; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane);

    FloatArray GivePitch();
    FloatArray GiveDerivativeUX(const FloatArray &lCoords);
    FloatArray GiveDerivativeVX(const FloatArray &lCoords);
    FloatArray GiveDerivativeUY(const FloatArray &lCoords);
    FloatArray GiveDerivativeVY(const FloatArray &lCoords);
    //virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual int testElementExtension(ElementExtension ext) { return 0; }
};
} // end namespace oofem
#endif //  trplanrot_h
