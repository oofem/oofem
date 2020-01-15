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

#ifndef latticelink3dboundary_h
#define latticelink3dboundary_h

#include "latticelink3d.h"

///@name Input fields for Latticelink3dboundary
//@{
#define _IFT_LatticeLink3dBoundary_Name "latticelink3dboundary"
#define _IFT_LatticeLink3dBoundary_location "location"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional lattice link element for the boundaries of a periodic cell.
 * The first two nodes have each 6 degrees of freedom (3 translation and 3 rotations).
 * One of the nodes  lies either on the boundary or outside the specimen.
 * This node is replaced with the corresponding image node. The control node (node 3) is used to impose the macroscopic (average) strain.
 * MACROSCOPIC INPUT for the control node: STRAIN TENSOR (3D, 6 COMPONENTS, VOIGT NOTATION: Exx Eyy Ezz Gyz Gxz Gxy)
 * @author: Peter Grassl
 */

class LatticeLink3dBoundary : public LatticeLink3d
{
protected:
    IntArray location;

public:
    LatticeLink3dBoundary(int n, Domain *);
    virtual ~LatticeLink3dBoundary();

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 18; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    void giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0) override;
    void computeGeometryProperties() override;

    void giveGPCoordinates(FloatArray &coords) override { coords = this->globalCentroid; }
    const char *giveInputRecordName() const override { return _IFT_LatticeLink3dBoundary_Name; }
    const char *giveClassName() const override { return "LatticeLink3dBoundary"; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const IntArray giveLocation() override { return location; };
    void recalculateCoordinates(int nodeNumber, FloatArray &coords) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif

protected:
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    void giveSwitches(IntArray &answer, int location);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN) override;
};
} // end namespace oofem
#endif
