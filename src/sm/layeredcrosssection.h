/* $Header: /home/cvs/bp/oofem/sm/src/layeredcrosssection.h,v 1.5 2003/04/14 16:01:01 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ***********************************
//   *** CLASS LAYERED CROSSSECTOION ***
//   ***********************************

#ifndef layeredcrosssection_h
#define layeredcrosssection_h

#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "element.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "interface.h"


class GaussPoint;
class LayeredCrossSectionModelInterface;

class LayeredCrossSection : public StructuralCrossSection
{
    /*
     * This class implements a layered cross section in a finite element problem. A cross
     * section  is an attribute of a domain. It is usually also attribute of many
     * elements.
     *
     * DESCRIPTION
     * The attribute 'propertyDictionary' contains all the properties of a
     * layered cross section, like thickness and width of each layer.
     * The atribute 'layerMaterials' contains an array of Materials corresponding
     * to each layer.
     *
     * it uses master - slave GaussPoint approach, where master gp has more slaves gp.
     * slave gp represent for each layer material point. It's coordinate sections
     * conteins z-coordinate (-1,1) from mid-section. the slaves are manageg completely
     * ( created, saved their context.,,,) from this class. Master gp only deletes
     * slaves in destructor.
     *
     * TASK
     * - Returning standard material stiffness marices (like 3dstress-strain, 2d plane ,
     * plate, 3dbeam, 2d beam ..) according to current state determined by parametr
     * StressMode by calling gp->material->GiveMaterialStiffnessMatrix (....) and by
     * possible modifiing returned matrix. (for example in layerde mode aproach
     * each layer  is asked for 3dMatrialStiffnes and this is integrated for example
     * over thickness for plate bending broblems)
     * - Returning RealStress state in gauss point and for given Stress mode.
     * - Returning a properties of cross section like thickness or area.
     */

protected:
    IntArray layerMaterials; // material of each layer
    FloatArray layerThicks; // thikness for each layer
    FloatArray layerWidths; // width for each layer
    int numberOfLayers;
    double midSurfaceZcoordFromBottom, totalThick;
    double area;

public:

    LayeredCrossSection(int n, Domain *d) : StructuralCrossSection(n, d), layerMaterials(), layerThicks(), layerWidths()
    { numberOfLayers = 0;
      totalThick = 0.;
      area = -1.0; }

    ~LayeredCrossSection()  { }

    void giveRealStresses(FloatArray & answer, MatResponseForm, GaussPoint *,
                          const FloatArray &, TimeStep * tStep);

    // updates gp - record
    // stressMode is stored in gp

    void giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseMode rMode,
                                         GaussPoint *,
                                         TimeStep *tStep);
    // next function is intendet to be used if we would like to obtain
    // char matrix form different material which is not associated with gp and its element.
    // (mainly for obtaining linear elastic matrix)
    // stress-strain mode is taken from gp.
    // NORMALLY - PLEASE USE GiveCharMaterialStiffnessMatrix function
    //
    virtual void giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                   MatResponseForm form, MatResponseMode rMode,
                                                   GaussPoint *, StructuralMaterial *,
                                                   TimeStep *tStep);



    void    giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *,
                                            const FloatArray &charVector3d);
    void    giveFullCharacteristicVector(FloatArray &answer,
                                         GaussPoint *, const FloatArray &);
    FloatArray *imposeStressConstrainsOnGradient(GaussPoint *, FloatArray *);
    FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *, FloatArray *);
    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                      MaterialMode mmode, StructuralMaterial *mat) const;
    virtual void giveLayerMaterialStiffnessMatrix(FloatMatrix &layerMatrix, MatResponseForm FullForm,
                                                  MatResponseMode rMode, GaussPoint *layerGp,
                                                  TimeStep *tStep);

    /**
     * Computes strain vector not dependent on sresses in given integration point. Returned vector is
     * generated by temperature or shrinkage effects, for example.
     * The load mode (Incremental or Total Load form) passed as parameter is taken into account.
     * Depend on load form, tre resulting strain is total strain or its increment from previous
     * step. Overloaded to support beams, paltes and shells.
     * @param answer stress independent strain vector
     * @param gp integration point
     * @param mode determines load mode
     * @param stepN time step (most models are able to respond only when atTime is current time step)
     * @param mode determines the response mode.
     */
    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *stepN, ValueModeType mode);


    double   give(int);

    // identification and auxiliary functions
    const char *giveClassName() const { return "LayeredCrossSection"; }
    classType giveClassID()         const { return LayeredCrossSectionClass; }
    IRResultType initializeFrom(InputRecord *ir);
    void     printYourself();
    double   computeIntegralThick();
    MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *, int);


    // store & restore context functions
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

#ifdef __PARALLEL_MODE
    /**
     * Pack all necessary data of integration point (according to element parallel_mode)
     * into given communication buffer. The corresponding material model service for particular integration point
     * is invoked. The nature of packed data is material model dependent.
     * Typically, for material of "local" response (response depeneds only on integration point local state)
     * no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
     * undergo averaging is performed between local and corressponding remote elements.
     * @param buff communication buffer
     * @param stepN solution step
     * @param ip integration point
     */
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    { _error("packUnknowns: not implemented");
      return 0; }
    /**
     * Unpack and updates all necessary data of given integration point (according to element parallel_mode)
     * into given communication buffer.
     * @see packUnknowns service.
     * @param buff communication buffer
     * @param stepN solution step.
     * @param ip integration point
     */
    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    { _error("unpackAndUpdateUnknowns: not implemented");
      return 0; }
    /**
     * Estimates the necessary pack size to hold all packed data of receiver.
     * The corresponding material model  service is invoked. The
     * nature of packed data is typically material model dependent.
     */
    int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
    { _error("estimatePackSize: not implemented");
      return 0; }
#endif

protected:
    virtual void giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                               MatResponseForm form,
                                               MatResponseMode rMode,
                                               GaussPoint *gp,
                                               StructuralMaterial *mat,
                                               TimeStep *tStep);
    void giveDerivedMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form,
                                            MatResponseMode rMode,
                                            GaussPoint *, StructuralMaterial *mat,
                                            TimeStep *tStep);

    void give2dPlateMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form,
                                            MatResponseMode rMode,
                                            GaussPoint *gp,
                                            StructuralMaterial *mat,
                                            TimeStep *tStep);
    void give3dShellMaterialStiffness(FloatMatrix &answer,
                                      MatResponseForm form,
                                      MatResponseMode rMode,
                                      GaussPoint *gp,
                                      StructuralMaterial *mat,
                                      TimeStep *tStep);
    void give2dBeamMaterialStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode rMode,
                                           GaussPoint *gp,
                                           StructuralMaterial *mat,
                                           TimeStep *tStep);

    FloatArray *GiveIntegrated3dShellStress(GaussPoint *gp);
    double       giveArea();



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
     * Computes full 3d strain vector in element layer. This function is necesary
     * if layered cross section is specified. If it is implemented, the testElementExtension
     * servise should return nonzero for Element_LayeredSupport parameter. This service is used by
     * layered cross section models.
     * @param answer full layer starin vector
     * @param masterGp element integration point
     * @param slaveGp slave integration point representing particular layer
     * @tStep time step
     */
    virtual void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                            GaussPoint *slaveGp, TimeStep *tStep) = 0;
};


#endif // layeredcrosssection_h


