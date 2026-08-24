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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#ifndef latticestructuralelement_h
#define latticestructuralelement_h

#include "sm/Elements/structuralelement.h"
#include "sm/Materials/LatticeMaterials/latticestructuralmaterial.h"

namespace oofem {
/**
 * This class implements the base of a special lattice element following
 * the concepts orginally developed by John Bolander. In this lattice
 * framework, elements consists of rigid bodies connected by a set of axial
 * tranversal and rotational springs at the position of the GP.
 * In the case of an irregular arrangement of lattice elements,
 * the position of the GP is not
 * placed at the midpoint of the element, but the midpoint of the element's midcross-section.
 * There is no way to relate the element geometry to the position of the GP. Instead, this information
 * is part of the input.
 * In this base class common interfaces of derived elements are defined.
 */
class LatticeStructuralElement : public StructuralElement
{
public:
    LatticeStructuralElement(int n, Domain *d);

    void initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority) override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    /**
     * Returns the cross-sectional area of the lattice element.
     * @return Cross-section area.
     */
    virtual double giveArea(GaussPoint *gp) { return 0; }


  virtual void giveSectionScaleFactors3d(FloatArray &q, GaussPoint *gp);

  virtual void giveSectionScaleFactors2d(FloatArray &q, GaussPoint *gp);

  virtual void convertStressToResultants3d(FloatArray &S, const FloatArray &sigma, GaussPoint *gp);

  virtual void convertStressToResultants2d(FloatArray &S, const FloatArray &sigma, GaussPoint *gp);

  virtual void convertTangentToResultantTangent3d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp);

  virtual void convertTangentToResultantTangent2d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp);


  virtual bool giveRectangularSectionDimensions(double &by, double &bz, GaussPoint *gp) const
  {
    return false;
  }

    /**
     * Returns the element length
     * @return Element length.
     */
    virtual double giveLength() { return 0; }

    /**
     * Returns the length associated with bond
     * @return bond length.
     */
    virtual double giveBondLength() { return 0; }

    /**
     * Returns the length associated with bond
     * @return bond length.
     */
    virtual double giveBondEndLength() { return 0; }

    /**
     * Returns the length associated with bond
     * @return bond length.
     */
    virtual double giveBondDiameter() { return 0; }

    /**
     * Returns the crack flag
     * @return Crack flag.
     */
    virtual int giveCrackFlag() { return 0; }

    /**
     * Returns the number of crossSection nodes
     * @return Number of crosssection nodes.
     */
    virtual int giveNumberOfCrossSectionNodes() { return 0; }

    /**
     * @return Crack width
     */
    virtual double giveCrackWidth() { return 0; }

    /**
     * @return oldCrackWidth
     */
    virtual double giveOldCrackWidth() { return 0; }

    /**
     * Returns the energy dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return dissipation
     */
    virtual double giveDissipation() { return 0; }

    /**
     * Returns pressures
     * @return pressures.
     */
    virtual void givePressures(FloatArray &pressures) { return; }

    /**
     * Returns plastic strains
     * @return plasticStrains
     */
    virtual void givePlasticStrain(FloatArray &plas) { return; }

    /**
     * Returns plastic strains
     * @return plasticStrains
     */
    virtual void giveOldPlasticStrain(FloatArray &plas) { return; }

    /**
     * Usually computes interpolation function, which is not needed for the lattice elements.
     * However, structural element requires implementation.
     */
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { return; }


    /**
     * Returns the increment of dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return increment of dissipation
     */
    virtual double giveDeltaDissipation() { return 0; }

    /**
     * Returns the coupling flag.
     * @return couplingFlag
     */

    virtual int giveCouplingFlag() { return 0; }

    /**
     * This function gives the cross-section coordinates.
     */
    virtual void giveCrossSectionCoordinates(FloatArray &coords) {; }

    /**
     * Returns the coupling numbers
     * @return couplingNumbers.
     */
    virtual void giveCouplingNumbers(IntArray &numbers) { }

    /**
     * Returns the normal stress.
     * @return normalStress
     */
    virtual double giveNormalStress() { return 0; }

    /**
     * Returns the temp (current, not-yet-committed) normal stress.
     * @return tempNormalStress
     */
    virtual double giveTempNormalStress() { return 0; }

    /**
     * Returns the old normal stress.
     * @return oldNormalStress
     */
    virtual double giveOldNormalStress() { return 0; }

    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

    /**
     * Returns flag indicating if status has been updated
     */
    virtual int hasBeenUpdated() { return 0; }


    /**
     * Gives the GP coordinates
     */
  virtual void  giveGpCoordinates(FloatArray &coords, GaussPoint *gp = nullptr) { }

  virtual double giveTributaryWidth(GaussPoint *gp)
{
    return 1.; // current element-constant value
}

    /**
     * Gives the y second moment of area
     */
    virtual double giveI1(GaussPoint *gp) { return 0.; }

    /**
     * Gives the I2 second moment of area
     */
    virtual double giveI2(GaussPoint *gp) { return 0.; }

    /**
     * Gives the torsional Constant
     */
    virtual double giveJ(GaussPoint *gp) { return 0.; }

    /**
     * Gives the shear area 1
     */
    virtual double giveShearArea1(GaussPoint *gp) { return 0.; }

    /**
     * Gives the shear area 2
     */
    virtual double giveShearArea2(GaussPoint *gp) { return 0.; }

protected:
    /// Rodrigues rotation matrix from a rotation (spin) vector. Shared by the
    /// geometrically nonlinear lattice elements (Lattice3dNL, LatticeLink3dNL).
    static void computeGlobalRotationMatrix(FloatMatrix &answer, const FloatArray &psi);
    /// Skew-symmetric matrix of a 3-vector (used by computeGlobalRotationMatrix).
    static void computeSMtrx(FloatMatrix &answer, const FloatArray &vec);
    /// Axial vector of log(R): the rotation vector of a rotation matrix (inverse of
    /// computeGlobalRotationMatrix). Robust at angle -> 0 and angle -> pi.
    static void logRotationVector(const FloatMatrix &R, FloatArray &theta);

};
} // end namespace oofem
#endif
