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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#ifndef libeam3dboundary_h
#define libeam3dboundary_h

#include "sm/Elements/Beams/libeam3d.h"
#include "sm/CrossSections/fiberedcs.h"
#include "nodalaveragingrecoverymodel.h"

///@name Input fields for LIBeam3dBoundary
//@{
#define _IFT_LIBeam3dBoundary_Name "libeam3dboundary"
#define _IFT_LIBeam3dBoundary_refnode "refnode"
#define _IFT_LIBeam3dBoundary_location "location"
//@}

namespace oofem {
/**
 * This class implements a boundary version of the 3-dimensional mindlin theory Linear Isoparametric
 * beam element, with reduced integration. Useful for prescribing periodicity in multiscale analyses.
 * MACROSCOPIC INPUT: DEFORMATION GRADIENT TENSOR (3D, 9 COMPONENTS: Exx Exy Exz Eyx Eyy Eyz Ezx Ezy Ezz)
 *
 * @author: Adam Sciegaj
 */
class LIBeam3dBoundary : public LIBeam3d
{
protected:
    double length;
    int referenceNode;
    IntArray location;

public:
    LIBeam3dBoundary(int n, Domain *d);
    virtual ~LIBeam3dBoundary() { }

    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return 21; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void recalculateCoordinates(int nodeNumber, FloatArray &coords) override;
    const IntArray giveLocation() override { return location; };

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LIBeam3dBoundary_Name; }
    const char *giveClassName() const override { return "LIBeam3dBoundary"; }

protected:
    virtual void computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep);
    void giveSwitches(IntArray &answer, int location);
    double computeLength() override;

    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // libeam3dboundary_h
