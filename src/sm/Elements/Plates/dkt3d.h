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

#ifndef dkt3d_h
#define dkt3d_h

#include "dkt.h"
#include "floatmatrix.h"

#define _IFT_DKTPlate3d_Name "dktplate3d"

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
 * This class represent DKT plate element that can be arbitrary
 * oriented in space, in contrast to base DKT element that is
 * defined in xy plane.
 * @see DKTPlate
 *
 * This class implements a triangular three-node plate DKT finite element.
 * Each node has 3 degrees of freedom.
 *
 */
class DKTPlate3d : public DKTPlate
{
protected:
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency.
     */
    FloatMatrix GtoLRotationMatrix;

public:
    DKTPlate3d(int n, Domain * d);
    virtual ~DKTPlate3d() {}

protected:
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    virtual void giveNodeCoordinates(double &x1, double &x2, double &x3,
                                     double &y1, double &y2, double &y3,
                                     double &z1, double &z2, double &z3);

    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode);

    /**
     * @name Edge load support
     */
    //@{
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    //@}
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
    virtual const char *giveClassName() const { return "DktPlate3d"; }
    virtual const char *giveInputRecordName() const { return _IFT_DKTPlate3d_Name; }

    virtual int computeNumberOfDofs() { return 9; }
    virtual int computeNumberOfGlobalDofs() { return 18; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual const FloatMatrix *computeGtoLRotationMatrix();
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);

    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual int giveLocalCoordinateSystem(FloatMatrix &answer)
    {
        OOFEM_ERROR("calling of this function id not allowed");
        return 0;
    }
    virtual int testElementExtension(ElementExtension ext)
    { return ( ( ext == Element_SurfaceLoadSupport )  ? 1 : 0 ); }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
};
} // end namespace oofem
#endif // dkt3d_h
