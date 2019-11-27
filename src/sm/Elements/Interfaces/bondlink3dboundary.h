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

#ifndef bondlink3dboundary_h
#define bondlink3dboundary_h

#include "bondlink3d.h"

///@name Input fields for BondLink3d
//@{
#define _IFT_BondLink3dBoundary_Name "bondlink3dboundary"
#define _IFT_LatticeLink3dBoundary_location "location"
//@}

namespace oofem {
/**
 * This class implements a bond link for connecting beam (frame) and continuum elements in unstructured meshes.
 * The main idea is to use the rotation of the beam element and the rigid arm from the beam node to the continuum element node
 * to compute the displacement jump along the rebar element (and two components, which are perpendicular to each other and lie
 * in a plane for which the direction along the rebar is normal to.
 * At least one node is located at the image boundary.
 * These nodes are replaced with a periodic mirror nodes and a control node is used to impose the macroscopic (average) strain.
 * MACROSCOPIC INPUT: DEFORMATION GRADIENT TENSOR (3D, 9 COMPONENTS: Exx Exy Exz Eyx Eyy Eyz Ezx Ezy Ezz)
 *
 */

class BondLink3dBoundary : public BondLink3d
{
protected:
    IntArray location;

public:
    BondLink3dBoundary(int n, Domain *);
    virtual ~BondLink3dBoundary();

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 18; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    void computeGeometryProperties() override;

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    const char *giveInputRecordName() const override { return _IFT_BondLink3dBoundary_Name; }
    const char *giveClassName()  const override { return "BondLink3dBoundary"; }
    void initializeFrom(InputRecord &ir) override;

    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif


protected:
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    virtual void computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep);
    void giveSwitches(IntArray &answer, int location);

    /**
     * This computes the geometrical properties of the element. It is called only once.
     */
    void computePropertiesOfCrossSection();

    integrationDomain giveIntegrationDomain() const override { return _Line; }
};
} // end namespace oofem
#endif //bondlink3dboundary_h
