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

#ifndef bondlink2dtruss_h
#define bondlink2dtruss_h

#include "../structuralelement.h"

///@name Input fields for BondLink2dTruss
//@{
#define _IFT_BondLink2dTruss_Name "bondlink2dtruss"
#define _IFT_BondLink2dTruss_length "length"
#define _IFT_BondLink2dTruss_diameter "diameter"
#define _IFT_BondLink2dTruss_dirvector "dirvector"
//@}

namespace oofem {
/**
 * This class implements a 2D bond link for connecting truss and continuum elements with coinciding nodes.
 * It is the 2D counterpart of BondLink3dTruss: node 1 is the truss node, node 2 the (hanging) continuum node.
 * The displacement jump u1 - u2 is expressed in a local system with the first axis along dirvector (slip)
 * and the second axis perpendicular to it in the x-y plane (lateral). The material is evaluated through the
 * ntt interface with the jump padded by a zero third component, so the same interface materials
 * (e.g. linkslip) can be used in 2D and 3D.
 *
 * @author: Peter Grassl
 */

class BondLink2dTruss : public StructuralElement
{
protected:
    double bondLength = 0.;

    FloatMatrix localCoordinateSystem;
    double bondDiameter = 0.;
    FloatArray directionVector;
    int geometryFlag = 0;
    FloatArray rigid;
    FloatArray globalCentroid;

public:
    BondLink2dTruss(int n, Domain *);

    double computeVolumeAround(GaussPoint *aGaussPoint) override;

    double giveLength();

    MaterialMode giveMaterialMode() override { return _2dInterface; }

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    /// Returns the global coordinates of the integration point, independent of lcoords.
    int computeGlobalCoordinates(Coordinates &answer, const FloatArray &lcoords) override;

    virtual double giveBondLength();

    virtual double giveBondDiameter();

    int computeNumberOfDofs() override { return 4; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    virtual void computeGeometryProperties();

    void giveInternalForcesVector(FloatArray &answer,
                                  TimeStep *tStep, int useUpdatedGpRecord) override;

    const char *giveInputRecordName() const override { return _IFT_BondLink2dTruss_Name; }
    const char *giveClassName()  const override { return "BondLink2dTruss"; }
    void initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority) override;

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
