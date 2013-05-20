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

#ifndef layeredcrosssection_h
#define layeredcrosssection_h

#include "structuralcrosssection.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "interface.h"

///@name Input fields for LayeredCrossSection
//@{
#define _IFT_LayeredCrossSection_Name "layeredcs"
#define _IFT_LayeredCrossSection_nlayers "nlayers"
#define _IFT_LayeredCrossSection_layermaterials "layermaterials"
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
    FloatArray layerThicks; ///< Thickness for each layer.
    FloatArray layerWidths; ///< Width for each layer.
    FloatArray layerMidZ;   ///< z-coord of the mid plane for each layer
    int numberOfLayers;
    int numberOfIntegrationPoints; ///< num integration points per layer
    double midSurfaceZcoordFromBottom;
    double midSurfaceXiCoordFromBottom;
    double totalThick;
    double area;

public:
    LayeredCrossSection(int n, Domain *d) : StructuralCrossSection(n, d), layerMaterials(), layerThicks(), layerWidths()
    {
        numberOfLayers = 0;
        totalThick = 0.;
        area = -1.0;
    }

    virtual ~LayeredCrossSection() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveRealStresses(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                  const FloatArray &reducedStrainIncrement, TimeStep *tStep);

    virtual void giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseMode mode,
                                                 GaussPoint *gp,
                                                 TimeStep *tStep);

    // next function is intended to be used if we would like to obtain
    // char matrix form different material which is not associated with gp and its element.
    // (mainly for obtaining linear elastic matrix)
    // stress-strain mode is taken from gp.
    // NORMALLY - PLEASE USE GiveCharMaterialStiffnessMatrix function
    virtual void giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                   MatResponseForm form, MatResponseMode rMode,
                                                   GaussPoint *, StructuralMaterial *,
                                                   TimeStep *tStep);



    virtual void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &charVector3d);
    virtual void giveFullCharacteristicVector(FloatArray &answer,
                                              GaussPoint *gp, const FloatArray &strainVector);

    virtual FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *);
    virtual FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *);

    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                      MaterialMode mmode, StructuralMaterial *mat) const;
    virtual void giveLayerMaterialStiffnessMatrix(FloatMatrix &layerMatrix, MatResponseForm form,
                                                  MatResponseMode rMode, GaussPoint *layerGp,
                                                  TimeStep *tStep);

    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual double give(CrossSectionProperty a);
    int giveNumberOfLayers();

    /// Returns the total thickness of all layers.
    double computeIntegralThick();
    void setupLayerMidPlanes();

    int giveLayerMaterial(int layer) { 
        return this->layerMaterials.at(layer); 
    }
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

    // identification and auxiliary functions
    virtual const char *giveInputRecordName() const { return _IFT_LayeredCrossSection_Name; }
    virtual const char *giveClassName() const { return "LayeredCrossSection"; }
    virtual classType giveClassID() const { return LayeredCrossSectionClass; }
    virtual void printYourself();

    MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *gp, int slaveIndex);

    virtual contextIOResultType saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);
    virtual contextIOResultType restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);


    void mapLayerGpCoordsToShellCoords(LayeredCrossSection *layeredCS, IntegrationRule **layerIntegrationRulesArray);

#ifdef __PARALLEL_MODE
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    {
        _error("packUnknowns: not implemented");
        return 0;
    }

    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    {
        _error("unpackAndUpdateUnknowns: not implemented");
        return 0;
    }

    int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
    {
        _error("estimatePackSize: not implemented");
        return 0;
    }
#endif

protected:
    virtual void giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                               MatResponseForm form,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               StructuralMaterial *mat,
                                               TimeStep *tStep);
    void giveDerivedMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form,
                                            MatResponseMode mode,
                                            GaussPoint *, StructuralMaterial *mat,
                                            TimeStep *tStep);

    void give2dPlateMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            StructuralMaterial *mat,
                                            TimeStep *tStep);
    void give3dShellMaterialStiffness(FloatMatrix &answer,
                                      MatResponseForm form,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      StructuralMaterial *mat,
                                      TimeStep *tStep);
    void give2dBeamMaterialStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           StructuralMaterial *mat,
                                           TimeStep *tStep);

    void giveIntegrated3dShellStress(FloatArray &answer, GaussPoint *gp);

    double giveArea();

    friend class Material;
};

/**
 * The element interface required by LayeredCrossSection.
 */
class LayeredCrossSectionInterface : public Interface
{
public:
    LayeredCrossSectionInterface() { }

    /**
     * Computes full 3D strain vector in element layer. This function is necessary
     * if layered cross section is specified. If it is implemented, the testElementExtension
     * service should return nonzero for Element_LayeredSupport parameter. This service is used by
     * layered cross section models.
     * @param answer Full layer strain vector.
     * @param masterGp Element integration point.
     * @param slaveGp Slave integration point representing particular layer.
     * @param tStep Time step.
     */
    virtual void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
        GaussPoint *slaveGp, TimeStep *tStep) {} ;
};
} // end namespace oofem
#endif // layeredcrosssection_h
