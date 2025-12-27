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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef latticelink3d_h
#define latticelink3d_h

#include "latticestructuralelement.h"

///@name Input fields for LatticeLink3d
//@{
#define _IFT_LatticeLink3d_Name "latticelink3d"
#define _IFT_LatticeLink3d_length "length"
#define _IFT_LatticeLink3d_diameter "diameter"
#define _IFT_LatticeLink3d_dirvector "dirvector"
#define _IFT_LatticeLink3d_l_end "l_end"
//@}

namespace oofem {

/**
 * This class implements a 3-dimensional lattice element
 */

class LatticeLink3d : public LatticeStructuralElement
{
protected:
    double bondLength;

    FloatMatrix localCoordinateSystem;
    double bondDiameter;
    FloatArray directionVector;
    int geometryFlag;
    double bondEndLength;
    FloatArray rigid;
    FloatArray globalCentroid;

    static ParamKey IPK_LatticeLink3d_length;
    static ParamKey IPK_LatticeLink3d_diameter;
    static ParamKey IPK_LatticeLink3d_dirvector;
    static ParamKey IPK_LatticeLink3d_l_end;

public:
    LatticeLink3d(int n, Domain *);
    virtual ~LatticeLink3d();

    double computeVolumeAround(GaussPoint *aGaussPoint) override;

    double giveLength() override;

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    double giveBondLength() override;

    double giveBondDiameter() override;

    double giveBondEndLength() override;

    int computeNumberOfDofs() override { return 12; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    virtual void giveGPCoordinates(FloatArray &coords);

    virtual void computeGeometryProperties();

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeLink3d_Name; }
    const char *giveClassName()  const override { return "LatticeLink3d"; }
    void initializeFrom(InputRecord &ir, int priority) override;
    void postInitialize() override;

    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif


protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;


    /**
     * This computes the geometrical properties of the element. It is called only once.
     */
    void computePropertiesOfCrossSection();

    void computeGaussPoints() override;
    integrationDomain giveIntegrationDomain() const override { return _Line; }
};
} // end namespace oofem
#endif
