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

#ifndef trplanrot3d_h
#define trplanrot3d_h

#include "../sm/Elements/PlaneStress/trplanrot.h"

#define _IFT_TrPlaneStrRot3d_Name "trplanestrrot3d"

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
    LocalMomentumTensor,
    GlobalMomentumTensor
};
#endif

/**
 * This class represent triangular plane stress element with rotational degree of freedom around normal
 * that can be arbitrary oriented in space, in contrast to base TrPlaneStrRot element that is
 * defined in xy-plane.
 *
 * Each node has 3 degrees of freedom.
 *
 * @author Ladislav Svoboda
 * @date 2010-5-25
 */
class TrPlaneStrRot3d : public TrPlaneStrRot
{
protected:
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix GtoLRotationMatrix;

public:
    TrPlaneStrRot3d(int n, Domain * d);
    virtual ~TrPlaneStrRot3d() { }

protected:
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    virtual void giveNodeCoordinates(FloatArray &x, FloatArray &y);

    virtual double computeVolumeAround(GaussPoint *gp);
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode);

    friend class TR_SHELL01;
    /**
     * @name Surface load support
     */
    //@{
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;
    virtual IntegrationRule *GetSurfaceIntegrationRule(int iSurf);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf);
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    //@}

public:
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TrPlaneStrRot3d_Name; }
    virtual const char *giveClassName() const { return "TrPlaneStrRot3d"; }

    virtual int computeNumberOfDofs() { return 9; }
    virtual int computeNumberOfGlobalDofs() { return 18; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    const FloatMatrix *computeGtoLRotationMatrix();
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    virtual int testElementExtension(ElementExtension ext)
    { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
};
} // end namespace oofem
#endif //  trplanrot3d_h
