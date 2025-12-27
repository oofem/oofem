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

#ifndef latticebeam3d_h
#define latticebeam3d_h

#include "latticestructuralelement.h"

///@name Input fields for Lattice3d
//@{
#define _IFT_LatticeBeam3d_Name "latticebeam3d"
#define _IFT_LatticeBeam3d_diameter "diameter"
//@}

namespace oofem {
class ParamKey;
/**
 * This class implements a 3-dimensional elastic bernoulli beam element in the lattice framework
 */

class LatticeBeam3d : public LatticeStructuralElement
{
protected:
    double kappa, length, diameter;
    double I1, I2, Ip;
    FloatMatrix localCoordinateSystem;
    double area;
    FloatArray midPoint, globalCentroid, normal;
    int geometryFlag;
    double myPi;

    static ParamKey IPK_LatticeBeam3d_diameter;

public:
    LatticeBeam3d(int n, Domain *);
    virtual ~LatticeBeam3d();


    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override {};

    virtual int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    virtual double giveLength() override;

    virtual double giveArea() override;

    virtual int computeNumberOfDofs() override { return 12; }

    virtual void giveDofManDofIDMask(int inode, IntArray &) const override;

    virtual double computeVolumeAround(GaussPoint *) override;

    virtual void giveGPCoordinates(FloatArray &coords);

    virtual void computeGeometryProperties();

    virtual void computeCrossSectionProperties();

    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeBeam3d_Name; }
    const char *giveClassName() const override { return "LatticeBeam3d"; }
    void initializeFrom(InputRecord &ir, int priority) override;
    void postInitialize() override;

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    virtual Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    virtual void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    virtual void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif


protected:
    virtual bool computeGtoLRotationMatrix(FloatMatrix &) override;
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    /**
     * This computes the geometrical properties of the element. It is called only once.
     */
    void computePropertiesOfCrossSection();

    virtual void computeGaussPoints() override;
    virtual integrationDomain  giveIntegrationDomain() const override{ return _Line; }
};
} // end namespace oofem
#endif
