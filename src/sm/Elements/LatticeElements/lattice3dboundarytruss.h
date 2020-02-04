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

#ifndef lattice3dboundarytruss_h
#define lattice3dboundarytruss_h

#include "lattice3dboundary.h"

///@name Input fields for Lattice3dBoundaryTruss
//@{
#define _IFT_Lattice3dBoundaryTruss_Name "lattice3dboundarytruss"
#define _IFT_Lattice3dBoundaryTruss_location "location"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional lattice element for the boundaries of a periodic cell representing a truss.
 * The first two nodes have each 6 degrees of freedom (3 translation and 3 rotations).
 * At least one nodes is located at the image boundary or outside the specimen.
 * These nodes are replaced with the corresponding image nodes and a control node is used to impose the macroscopic (average) strain.
 * The third node has 1 degree of fredom.
 * MACROSCOPIC INPUT: Axial strain (Exx)
 *
 * @author: Peter Grassl
 */

class Lattice3dBoundaryTruss : public Lattice3dBoundary
{
protected:
    IntArray location;

public:
    Lattice3dBoundaryTruss(int n, Domain *);
    virtual ~Lattice3dBoundaryTruss();

    void recalculateCoordinates(int nodeNumber, FloatArray &coords) override;
    const IntArray giveLocation() override { return this->location; }

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 13; }

    double computeVolumeAround(GaussPoint *) override;

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    void  giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0) override;
    void computeGeometryProperties() override;

    void giveGPCoordinates(FloatArray &coords) override { coords = this->globalCentroid; }
    const char *giveInputRecordName() const override { return _IFT_Lattice3dBoundaryTruss_Name; }
    const char *giveClassName() const override { return "Lattice3dBoundaryTruss"; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }
    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    void drawRawCrossSections(oofegGraphicContext &, TimeStep *tStep);
    void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif

protected:
    void giveSwitches(IntArray &answer, int location);
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN) override;
};
} // end namespace oofem
#endif
