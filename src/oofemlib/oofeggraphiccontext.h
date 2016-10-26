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


#ifndef oofeggraphiccontext_h
#define oofeggraphiccontext_h

#include "oofemcfg.h"
#include "intarray.h"
#include "internalstatemode.h"
#include "internalstatetype.h"

#include <list>

//
// for c++ compiler to be successful on some c files
//

extern "C" {
 #define new __new
 #define class __class
 #define inline __inline
 #define min __min
 #define max __max
 #define sgn __sgn
 #define macbra __macbra
 #define Status __Status
 #define Request __Request
 #define Complex __Complex

 #include <Esimple.h>

 #undef new
 #undef class
 #undef inline
 #undef min
 #undef max
 #undef sgn
 #undef macbra
 #undef Status
 #undef Request
 #undef PI
 #undef Complex
 #undef None
};



// not declared in any *.h Elixir file

extern "C" {
    void EVFastRedraw(EView *v_p);
};

namespace oofem {
// width definition variables
 #define OOFEG_RAW_GEOMETRY_WIDTH 0
 #define OOFEG_DEFORMED_GEOMETRY_WIDTH 0
 #define OOFEG_CRACK_PATTERN_WIDTH 2
 #define OOFEG_ISO_LINE_WIDTH 4
 #define OOFEG_SPARSE_PROFILE_WIDTH 0

// layer definition variables
 #define OOFEG_RAW_GEOMETRY_LAYER         0
 #define OOFEG_DEFORMED_GEOMETRY_LAYER    1
 #define OOFEG_NODE_ANNOTATION_LAYER      2
 #define OOFEG_ELEMENT_ANNOTATION_LAYER   3
 #define OOFEG_VARPLOT_PATTERN_LAYER      4
 #define OOFEG_CRACK_PATTERN_LAYER        5
 #define OOFEG_BCIC_ANNOTATION_LAYER      6
 #define OOFEG_NATURALBC_LAYER            7
 #define OOFEG_SPARSE_PROFILE_LAYER       8
 #define OOFEG_DEBUG_LAYER                9
 #define OOFEG_RAW_CROSSSECTION_LAYER      10
 #define OOFEG_LAST_LAYER                11
class EngngModel;
class Element;
class Range;

enum OGC_PlotModeType {
    OGC_unknown, OGC_rawGeometry, OGC_deformedGeometry, OGC_eigenVectorGeometry,
    OGC_nodeGeometry, OGC_nodeAnnotation, OGC_elementAnnotation, OGC_essentialBC, OGC_naturalBC,
    OGC_nodeScalarPlot, OGC_nodeVectorPlot,
    OGC_scalarPlot, OGC_vectorPlot, OGC_tensorPlot,
    OGC_elemSpecial
};


enum ScalarAlgorithmType { SA_ISO_SURF, SA_ISO_LINE, SA_ZPROFILE, SA_COLORZPROFILE };
enum SmootherType { Smoother_NA, Smoother_ZZ, Smoother_SPR };
enum ScaleMode { SM_Autoscale, SM_UserDefined };

 #define OOFEG_YIELD_STEPS 3

class OOFEM_EXPORT oofegGraphicContext
{
protected:

    /*
     * Common members to all contexts
     */
    static EngngModel *emodel;
    static EFringeTable ft;

    static EPixel meshFillColor;
    static EPixel remoteMeshFillColor;
    static EPixel edgeColor;
    static EPixel remoteEdgeColor;

    static EPixel nodeColor;
    static EPixel bcicColor;
    static EPixel bcForceColor;
    static EPixel deformedElementColor;
    static EPixel crackPatternColor;
    static EPixel activeCrackColor;
    static EPixel yieldPlotColors [ OOFEG_YIELD_STEPS ];
    static EPixel standardSparseProfileColor, extendedSparseProfileColor;
    static EPixel geometryColor;
    static EPixel crossSectionColor;

    static int activeStep, activeStepVersion;
    static double defScale;
    static double zprofilescale; // for landscape plots in 2d
    static int activeEigVal;
    static int activeYieldStep;

    // Material Model (region) filter
    static IntArray matRegFilter;

    // element filter
    static std :: list< Range >element_filter;

    static ScalarAlgorithmType scalarAlgo;
    // smoother type
    static SmootherType smootherType;

    // deformed geometry internal variables plot flag
    // nonzero indicate to use deformed shape
    static int intVarDefGeoFlag;

    // Sparse matrix profile mode (1=marker mode, 0=grid mode)
    static int sparseProfileMode;

    // EFringeTable ft;
    static int activeProblem;
    static int activeDomain;

    // scale mode and color scale values
    static ScaleMode smode;
    static double emin, emax;
    static int scaleInitFlag;

    static bool staticVarsInitFlag;

    // mode of value
    static InternalStateMode varMode;

    /*
     * Attributes
     */
    // current value to display
    InternalStateType varType;
    // component of value
    int component;
    // plot mode (scalar, vector, tensor plot)
    OGC_PlotModeType plotMode;
    // on/off flag
    bool isActiveFlag;

public:
    oofegGraphicContext();
    ~oofegGraphicContext();

    void init(EngngModel *);


    EPixel getElementColor() { return meshFillColor; }
    EPixel getElementEdgeColor() { return edgeColor; }
    EPixel getRemoteElementColor() { return remoteMeshFillColor; }
    EPixel getRemoteElementEdgeColor() { return remoteEdgeColor; }
    EPixel getNodeColor()    { return nodeColor; }
    EPixel getBcIcColor()    { return bcicColor; }
    EPixel getBcForceColor()    { return bcForceColor; }
    EPixel getDeformedElementColor() { return deformedElementColor; }
    EPixel getCrackPatternColor() { return crackPatternColor; }
    EPixel getActiveCrackColor() { return activeCrackColor; }
    EPixel getYieldPlotColor(double ratio)
    { return this->GR_giveColorFromUserColorTable(yieldPlotColors, OOFEG_YIELD_STEPS, ratio); }
    EPixel getCrossSectionColor()    { return crossSectionColor; }
    EPixel getStandardSparseProfileColor() { return standardSparseProfileColor; }
    EPixel getExtendedSparseProfileColor() { return extendedSparseProfileColor; }
    EPixel getGeometryColor() { return geometryColor; }
    int    getSparseProfileMode() { return sparseProfileMode; }


    // DrawMode getDrawMode () {return mode;}
    // EFringeTable getFringeTable () {return ft;}
    int    getActiveStep() { return activeStep; }
    int    getActiveStepVersion() { return activeStepVersion; }
    double getDefScale()   { return defScale; }
    double getLandScale()   { return zprofilescale; }
    int    getActiveEigVal() { return activeEigVal; }
    int    getActiveYieldStep() { return activeYieldStep; }
    int    getInternalVarsDefGeoFlag() { return intVarDefGeoFlag; }
    int    getActiveDomain() { return activeDomain; }
    int    getActiveProblemIndx() { return activeProblem; }
    EngngModel *getActiveProblem();
    EFringeTable getFringeTable() { return ft; }

    void setElementColor(EPixel color) { meshFillColor = color; }
    void setElementEdgeColor(EPixel color) { edgeColor = color; }
    void setNodeColor(EPixel color) { nodeColor    = color; }
    void setDeformedElementColor(EPixel color) { deformedElementColor = color; }
    void setCrackPatternColor(EPixel color) { crackPatternColor = color; }
    void setActiveCrackColor(EPixel color) { activeCrackColor  = color; }
    void setGeometryColor(EPixel color) { geometryColor = color; }
    void setCrossSectionColor(EPixel color) { crossSectionColor = color; }
    void setActiveStep(int n) { activeStep = n; }
    void setActiveStepVersion(int n) { activeStepVersion = n; }
    void setDefScale(double n)   { defScale = n; }
    void setLandScale(double n)   { zprofilescale = n; }
    void setActiveEigVal(int n) { activeEigVal = n; }
    void setActiveYieldStep(int n) { activeYieldStep = n; }
    void setInternalVarsDefGeoFlag(int n) { intVarDefGeoFlag = n; }
    void setSparseProfileMode(int n) { sparseProfileMode = n; }
    int  setActiveDomain(int a) {
        activeDomain = a;
        return activeDomain;
    }
    int  setActiveProblem(int a);

    void setPlotMode(OGC_PlotModeType mode) { plotMode = mode; }
    void setInternalStateMode(InternalStateMode mode) { varMode = mode; }
    void setInternalStateType(InternalStateType type) { varType = type; }
    void setIntVarIndx(int indx) { component = indx; }
    InternalStateType giveIntVarType() { return this->varType; }
    InternalStateMode giveIntVarMode() { return this->varMode; }
    OGC_PlotModeType giveIntVarPlotMode() { return this->plotMode; }
    int               giveIntVarIndx() { return this->component; }

    ScalarAlgorithmType getScalarAlgo() { return scalarAlgo; }
    void setScalarAlgo(ScalarAlgorithmType a) { scalarAlgo = a; }

    SmootherType giveSmootherType()  { return smootherType; }
    void setSmootherType(SmootherType type) { this->smootherType = type; }

    ScaleMode getScaleMode() { return smode; }
    void      setScaleMode(ScaleMode s) { smode = s; }
    double    getScaleMin() { return emin; }
    double    getScaleMax() { return emax; }
    void      setScaleVals(double smin, double smax) {
        emin = smin;
        emax = smax;
    }
    void      resetScaleVals() {
        if ( smode == SM_Autoscale ) {
            scaleInitFlag = 1;
            emin = 1.0;
            emax = -1.0;
        }
    }
    void      updateFringeTableMinMax(double *s, int size);

    // component filters
    // element filter
    /**
     * Test if particular element passed fulfills various filtering criteria for its graphics output.
     * @return Nonzero if output requested, zero otherwise..
     */
    int testElementGraphicActivity(Element *);
    /**
     * Returns the state of material model (region) filter for particular mod.el
     */
    int getMaterialModelFilterState(int i);
    /**
     * Sets the state of material model (region) filter for particular model.
     */
    void setMaterialModelFilterState(int i, int state);
    /**
     * Sets the state of element filter for particular model
     * @param initString String containing valid range string representation, with element_filter keyword.
     */
    void setElementFilterState(char *initString);

    /** tests if context is active */
    bool isActive() { return this->isActiveFlag; }
    /// sets activity flag
    void setActivityFlag(bool flag) { isActiveFlag = flag; }


protected:
    EPixel GR_giveColorFromUserColorTable(EPixel *table, int tableSize, double relVal);
    // void GR_setupUserColors ();
    // void GR_deleteUserColorTables();
};
} // end namespace oofem

/* these reside in global name space */
extern oofem :: oofegGraphicContext gc [ OOFEG_LAST_LAYER ];
extern EView *myview;
extern void deleteLayerGraphics(int iLayer);

#endif // oofeggraphiccontext_h
