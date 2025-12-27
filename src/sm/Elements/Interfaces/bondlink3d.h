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

#ifndef bondlink3d_h
#define bondlink3d_h

#include "../structuralelement.h"

///@name Input fields for BondLink3d
//@{
#define _IFT_BondLink3d_Name "bondlink3d"
#define _IFT_BondLink3d_length "length"
#define _IFT_BondLink3d_diameter "diameter"
#define _IFT_BondLink3d_dirvector "dirvector"
#define _IFT_BondLink3d_length_end "length_end"
//@}

namespace oofem {
class ParamKey;
/**
 * This class implements a bond link for connecting beam (frame) and continuum elements in unstructured meshes.
 * The main idea is to use the rotation of the beam element and the rigid arm from the beam node to the continuum element node
 * to compute the displacement jump along the rebar element (and two components, which are perpendicular to each other and lie
 * in a plane for which the direction along the rebar is normal to.
 * This element differs from the lattice link element, for which both the beam and the lattice nodes have rotational DOFs and
 * therefore the lattice node's rotations can be used to compute the displacement jump at the beam element location.
 *
 * @author: Peter Grassl
 */

class BondLink3d : public StructuralElement
{
protected:
    double bondLength = 0.;

    FloatMatrix localCoordinateSystem;
    double bondDiameter = 0.;
    FloatArray directionVector;
    int geometryFlag = 0;
    double bondEndLength = 0.;
    FloatArray rigid;
    FloatArray globalCentroid;

    static ParamKey IPK_BondLink3d_length;
    static ParamKey IPK_BondLink3d_diameter;
    static ParamKey IPK_BondLink3d_dirvector;
    static ParamKey IPK_BondLink3d_length_end;

public:
    BondLink3d(int n, Domain *);

    double computeVolumeAround(GaussPoint *aGaussPoint) override;

    double giveLength();

    MaterialMode giveMaterialMode() override { return _3dInterface; }

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    /**
     * This function is different from the standard computeGlobalCordinates
     * function as it returns the global coordinates of the gausspoint
     * independent to the value of the lcoords.
     */
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    virtual double giveBondLength();

    virtual double giveBondDiameter();

    virtual double giveBondEndLength();

    int computeNumberOfDofs() override { return 9; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    virtual void giveGPCoordinates(FloatArray &coords);

    virtual void computeGeometryProperties();

    void giveInternalForcesVector(FloatArray &answer,
                                  TimeStep *tStep, int useUpdatedGpRecord) override;

    const char *giveInputRecordName() const override { return _IFT_BondLink3d_Name; }
    const char *giveClassName()  const override { return "BondLink3d"; }
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
    integrationDomain  giveIntegrationDomain() const override { return _Line; }
};
} // end namespace oofem
#endif
