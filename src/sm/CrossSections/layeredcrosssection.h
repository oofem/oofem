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

#include "../sm/CrossSections/structuralcrosssection.h"
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
 * contains z-coordinate (-1,1) from mid-section. the slaves are manage completely
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
    int numberOfLayers;
    int numberOfIntegrationPoints; ///< num integration points per layer
    double midSurfaceZcoordFromBottom;
    double midSurfaceXiCoordFromBottom;
    double totalThick;
    double area;

public:
    LayeredCrossSection(int n, Domain * d) : StructuralCrossSection(n, d), layerMaterials(), layerThicks(), layerWidths()
    {
        numberOfLayers = 0;
        totalThick = 0.;
        area = -1.0;
    }

    virtual ~LayeredCrossSection() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void createMaterialStatus(GaussPoint &iGP); // ES

    virtual int setupIntegrationPoints(IntegrationRule &irule, int npoints, Element *element);

    virtual void giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveStiffnessMatrix_PlaneStress(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveStiffnessMatrix_PlaneStrain(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveStiffnessMatrix_1d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    virtual void giveGeneralizedStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_PlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);

    virtual void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode mode);

    virtual void give2dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give2dPlateStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dShellStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *);
    virtual FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *);

    virtual double give(CrossSectionProperty a, GaussPoint *gp);
    virtual double give(CrossSectionProperty a, const FloatArray &coords, Element *elem, bool local);
    int giveNumberOfLayers();

    /// Returns the total thickness of all layers.
    double computeIntegralThick();
    void setupLayerMidPlanes();

    int giveLayerMaterial(int layer) {
        return this->layerMaterials.at(layer);
    }

    virtual Material *giveMaterial(IntegrationPoint *ip);

    int giveInterfaceMaterialNum(int interface) {
        return this->interfacerMaterials.at(interface);
    }

    Material *giveInterfaceMaterial(int interface) {
        int matNum = this->giveInterfaceMaterialNum(interface);
        if ( matNum ) {
            return this->giveDomain()->giveMaterial( this->interfacerMaterials.at(interface) );
        } else {
            return NULL;
        }
    }

    virtual int checkConsistency();

    double giveLayerMidZ(int layer) {
        // Gives the z-coord measured from the geometric midplane of the (total) cross section.
        return this->layerMidZ.at(layer);
    }
    double giveLayerThickness(int layer) {
        return this->layerThicks.at(layer);
    }
    int giveNumIntegrationPointsInLayer() {
        return this->numberOfIntegrationPoints;
    }
    double giveMidSurfaceZcoordFromBottom() {
        return this->midSurfaceZcoordFromBottom;
    }
    double giveMidSurfaceXiCoordFromBottom() {
        return this->midSurfaceXiCoordFromBottom;
    }
    void giveInterfaceXiCoords(FloatArray &answer);

    // identification and auxiliary functions
    virtual const char *giveInputRecordName() const { return _IFT_LayeredCrossSection_Name; }
    virtual const char *giveClassName() const { return "LayeredCrossSection"; }
    virtual void printYourself();

    MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode mode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *gp, int slaveIndex);

    virtual contextIOResultType saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);
    virtual contextIOResultType restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);


    void mapLayerGpCoordsToShellCoords(std :: vector< std :: unique_ptr< IntegrationRule > > &layerIntegrationRulesArray);

    void setupLayeredIntegrationRule(std :: vector< std :: unique_ptr< IntegrationRule > > &layerIntegrationRulesArray, Element *el, int numInPlanePoints);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep);
    virtual double give(int aProperty, GaussPoint *gp);

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    virtual int estimatePackSize(DataStream &buff, GaussPoint *ip)
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }


    virtual void giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep)
    { OOFEM_ERROR("not implemented"); }
    virtual void giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep)
    { OOFEM_ERROR("not implemented"); }
    virtual void giveStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
    { OOFEM_ERROR("not implemented"); }
    virtual void giveStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
    { OOFEM_ERROR("not implemented"); }


protected:
    double giveArea();
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
     * @param masterGpStrain Generalized strain at maters gauss point.
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
    virtual ~LayeredIntegrationRule();

    virtual const char *giveClassName() const { return "LayeredIntegrationRule"; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    //virtual int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder);

    // Stores the ip numbers of the points lying on the lower and upper surface of the element.
    // Thus they will correspond to points lying on the interface between layers.
    IntArray lowerInterfacePoints, upperInterfacePoints;

    //virtual int SetUpPointsOnLine(int, MaterialMode);      // could be used for beams
    //virtual int SetUpPointsOnCube(int, MaterialMode mode); // could be used for plates/shells/solids
    virtual int SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode);
};
} // end namespace oofem
#endif // layeredcrosssection_h
