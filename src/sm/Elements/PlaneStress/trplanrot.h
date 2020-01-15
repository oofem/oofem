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

#include "sm/Elements/PlaneStress/trplanstrss.h"

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
    void computeGaussPoints() override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;

    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    double giveArea() override;
    virtual void giveNodeCoordinates(FloatArray &x, FloatArray &y);

    void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode) override;

public:
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TrPlaneStrRot_Name; }
    const char *giveClassName() const override { return "TrPlaneStrRot"; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _PlaneStressRot; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    int computeNumberOfDofs() override { return 9; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override;

    FloatArray GivePitch();
    FloatArray GiveDerivativeUX(const FloatArray &lCoords);
    FloatArray GiveDerivativeVX(const FloatArray &lCoords);
    FloatArray GiveDerivativeUY(const FloatArray &lCoords);
    FloatArray GiveDerivativeVY(const FloatArray &lCoords);
    //void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

    int testElementExtension(ElementExtension ext) override { return 0; }
};
} // end namespace oofem
#endif //  trplanrot_h
