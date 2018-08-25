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

#ifndef trplanestressrotallman3d_h
#define trplanestressrotallman3d_h

#include "sm/Elements/PlaneStress/trplanestressrotallman.h"

#define _IFT_TrPlanestressRotAllman3d_Name "trplanestressrotallman3d"

namespace oofem {
#ifndef __CHARTENSOR
 #define __CHARTENSOR
enum CharTensor {
    LocalStrainTensor,
    GlobalStrainTensor,
    LocalCurvatureTensor,
    GlobalCurvatureTensor,

    LocalForceTensor,
    GlobalForceTensor,
    LocalMomentTensor,
    GlobalMomentTensor
};
#endif

/**
 * This class represent triangular plane stress element with rotational degree of freedom around normal
 * that can be arbitrary oriented in space, in contrast to base TrPlanestressRotAllman element that is
 * defined in xy-plane.
 * 
 *
 */
class TrPlanestressRotAllman3d : public TrPlanestressRotAllman
{
protected:
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix *GtoLRotationMatrix;

public:
    TrPlanestressRotAllman3d(int n, Domain * d);
    virtual ~TrPlanestressRotAllman3d() {
        delete GtoLRotationMatrix;
    }

protected:
    void computeLocalNodalCoordinates(std::vector< FloatArray > &lxy) override;

    double computeVolumeAround(GaussPoint *gp) override;
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);

    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode) override;

public:
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TrPlanestressRotAllman3d_Name; }
    const char *giveClassName() const override { return "TrPlaneStressRotAllman3d"; }

    int computeNumberOfDofs() override { return 9; }
    int computeNumberOfGlobalDofs() override { return 18; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    int testElementExtension(ElementExtension ext) override
    { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    void printOutputAt(FILE *file, TimeStep *tStep) override;
};
} // end namespace oofem
#endif //  trplanestressrotallman3d_h
