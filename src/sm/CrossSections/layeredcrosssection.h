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

#ifndef layeredcrosssection_h
#define layeredcrosssection_h

#include "sm/CrossSections/structuralcrosssection.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "interface.h"
#include "gaussintegrationrule.h"
#include "domain.h"

#include <vector>
#include <memory>

///@name Input fields for LayeredCrossSection
//@{
#define _IFT_LayeredCrossSection_Name "layeredcs"
#define _IFT_LayeredCrossSection_nlayers "nlayers"
#define _IFT_LayeredCrossSection_layermaterials "layermaterials"
#define _IFT_LayeredCrossSection_interfacematerials "interfacematerials"
#define _IFT_LayeredCrossSection_layerRotations "rotations"
#define _IFT_LayeredCrossSection_thicks "thicks"
#define _IFT_LayeredCrossSection_widths "widths"
#define _IFT_LayeredCrossSection_midsurf "midsurf"
#define _IFT_LayeredCrossSection_nintegrationpoints "nintegrationpoints"
#define _IFT_LayeredCrossSection_initiationlimits "initiationlimits"
#define _IFT_LayeredCrossSection_shearcoeff_xz "beamshearcoeffxz"
//@}

namespace oofem {
class StructuralMaterial;

/**
 * This class implements a layered cross section in a finite element problem. A cross
 * section  is an attribute of a domain. It is usually also attribute of many
 * elements.
 *
 * The attribute 'propertyDictionary' contains all the properties of a
 * layered cross section, like thickness and width of each layer.
 * The attribute 'layerMaterials' contains an array of Materials corresponding
 * to each layer.
 *
 * It uses master - slave GaussPoint approach, where master gp has more slaves gp.
 * slave gp represent for each layer material point. It's coordinate sections
 * contains z-coordinate (-1,1) from mid-section. The slaves are manage completely
 * ( created, saved their context.,,,) from this class. Master gp only deletes
 * slaves in destructor.
 *
 * Tasks:
 * - Returning standard material stiffness matrices (like 3d stress-strain, 2d plane ,
 *   plate, 3d beam, 2d beam ..) according to current state determined by parameter
 *   StressMode by calling material and by
 *   possible modifying returned matrix (for example in layered mode approach
 *   each layer is asked for 3dMaterialStiffness and this is integrated for example
 *   over thickness for plate bending problems).
 * - Returning RealStress state in Gauss point and for given Stress mode.
 * - Returning a properties of cross section like thickness or area.
 */
class LayeredCrossSection : public StructuralCrossSection
{
protected:
    IntArray layerMaterials; ///< Material of each layer.
    IntArray interfacerMaterials; ///< Interface (cohesive zone) material for each interface.
    FloatArray layerThicks; ///< Thickness for each layer.
    FloatArray layerWidths; ///< Width for each layer.
    FloatArray layerMidZ;   ///< z-coord of the mid plane for each layer
    FloatArray layerRots;   ///< Rotation of the material in each layer.
    int numberOfLayers = 0;
    int numberOfIntegrationPoints = 1; ///< num integration points per layer
    double midSurfaceZcoordFromBottom = 0.;
    double midSurfaceXiCoordFromBottom = 0.;
    double totalThick = 0.;
    double area = 0.;
  double beamShearCoeffxz = 1.0;
public:
    LayeredCrossSection(int n, Domain * d) : 
        StructuralCrossSection(n, d)
    { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void createMaterialStatus(GaussPoint &iGP) override;

    //Create slave integration points
    int setupIntegrationPoints(IntegrationRule &irule, int npoints, Element *element) override;
    /**
     * Sets up integration rule for the given element.
     * Default behavior is just to call the Gauss integration rule, but for example the layered and fibered crosssections need to do their own thing.
     * @param irule Integration rule to set up.
     * @param npointsXY Number of integration points in xi-eta.
     * @param npointsZ Number of integration points in zeta.
     * @param element Element which the integration rule belongs to.
     * @return Number of integration points.
     */
    int setupIntegrationPoints(IntegrationRule &irule, int npointsXY, int npointsZ, Element *element) override;

    FloatArrayF<6> giveRealStress_3d(const FloatArrayF<6> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<6> giveRealStress_3dDegeneratedShell(const FloatArrayF<6> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<4> giveRealStress_PlaneStrain(const FloatArrayF<4> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveRealStress_PlaneStress(const FloatArrayF<3> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<1> giveRealStress_1d(const FloatArrayF<1> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<2> giveRealStress_Warping(const FloatArrayF<2> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<6,6> giveStiffnessMatrix_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> giveStiffnessMatrix_PlaneStress(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> giveStiffnessMatrix_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<1,1> giveStiffnessMatrix_1d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<3> giveGeneralizedStress_Beam2d(const FloatArrayF<3> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<6> giveGeneralizedStress_Beam3d(const FloatArrayF<6> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<5> giveGeneralizedStress_Plate(const FloatArrayF<5> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<8> giveGeneralizedStress_Shell(const FloatArrayF<8> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<4> giveGeneralizedStress_MembraneRot(const FloatArrayF<4> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveGeneralizedStress_PlateSubSoil(const FloatArrayF<3> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;

    void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode mode) const override;

    FloatMatrixF<3,3> give2dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> give3dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<5,5> give2dPlateStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<8,8> give3dShellStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> give3dDegeneratedShellStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> giveMembraneRotStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give2dPlateSubSoilStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *) override;
    FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *) override;

    double give(CrossSectionProperty a, GaussPoint *gp) const override;
    double give(CrossSectionProperty a, const FloatArray &coords, Element *elem, bool local) const override;
    int giveNumberOfLayers() const;
	int giveLayer(GaussPoint *gp) const;

    /// Returns the total thickness of all layers.
    double computeIntegralThick() const;
    void setupLayerMidPlanes();

    int giveLayerMaterial(int layer) const {
        return this->layerMaterials.at(layer);
    }

    Material *giveMaterial(IntegrationPoint *ip) const override;

    int giveInterfaceMaterialNum(int interface) {
        return this->interfacerMaterials.at(interface);
    }

    Material *giveInterfaceMaterial(int interface) {
        int matNum = this->giveInterfaceMaterialNum(interface);
        if ( matNum ) {
            return this->giveDomain()->giveMaterial( this->interfacerMaterials.at(interface) );
        } else {
            return nullptr;
        }
    }

    int checkConsistency() override;

    double giveLayerMidZ(int layer) const {
        // Gives the z-coord measured from the geometric midplane of the (total) cross section.
        return this->layerMidZ.at(layer);
    }
    double giveLayerThickness(int layer) const {
        return this->layerThicks.at(layer);
    }
    int giveNumIntegrationPointsInLayer() const {
        return this->numberOfIntegrationPoints;
    }
    double giveMidSurfaceZcoordFromBottom() const {
        return this->midSurfaceZcoordFromBottom;
    }
    double giveMidSurfaceXiCoordFromBottom() const {
        return this->midSurfaceXiCoordFromBottom;
    }
    void giveInterfaceXiCoords(FloatArray &answer) const;

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_LayeredCrossSection_Name; }
    const char *giveClassName() const override { return "LayeredCrossSection"; }
    void printYourself() override;

    static MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode mode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *gp, int layer, int igp) const;

    void saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;
    void restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;


    void mapLayerGpCoordsToShellCoords(std :: vector< std :: unique_ptr< IntegrationRule > > &layerIntegrationRulesArray);

    void setupLayeredIntegrationRule(std :: vector< std :: unique_ptr< IntegrationRule > > &layerIntegrationRulesArray, Element *el, int numInPlanePoints);

    int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep) override;
    double give(int aProperty, GaussPoint *gp) const override;

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    int estimatePackSize(DataStream &buff, GaussPoint *ip) override
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }


    void giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    void giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    void giveStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    void giveStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }

protected:
    double giveArea() const;
};

/**
 * The element interface required by LayeredCrossSection.
 */
class LayeredCrossSectionInterface : public Interface
{
public:
    LayeredCrossSectionInterface() { }

    /**
     * Computes full 3D strain vector in element layer.
     * This function is necessary if layered cross section is specified..
     * @param answer Full layer strain vector.
     * @param masterGpStrain Generalized strain at master gauss point.
     * @param masterGp Element integration point.
     * @param slaveGp Slave integration point representing particular layer.
     * @param tStep Time step.
     */
    virtual void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain, GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep) = 0;
};

class LayeredIntegrationRule : public IntegrationRule
{
public:
    LayeredIntegrationRule(int n, Element * e, int startIndx, int endIndx, bool dynamic = false);
    LayeredIntegrationRule(int n, Element * e);

    const char *giveClassName() const override { return "LayeredIntegrationRule"; }

    //virtual int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder);

    // Stores the ip numbers of the points lying on the lower and upper surface of the element.
    // Thus they will correspond to points lying on the interface between layers.
    IntArray lowerInterfacePoints, upperInterfacePoints;

    //int SetUpPointsOnLine(int, MaterialMode) override;      // could be used for beams
    //int SetUpPointsOnCube(int, MaterialMode mode) override; // could be used for plates/shells/solids
    int SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode) override;
};
} // end namespace oofem
#endif // layeredcrosssection_h
