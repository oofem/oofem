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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifdef __OOFEG
#include "compiler.h"
#include "oofeggraphiccontext.h"
#include "engngm.h"
#include "element.h"
#include "material.h"
#include "mathfem.h"
// for Range class definition outputmanager.h included
#include "outputmanager.h"
#include "strreader.h"
#include "util.h"

namespace oofem {
EngngModel *oofegGraphicContext :: emodel = NULL;
EFringeTable oofegGraphicContext :: ft;
EPixel oofegGraphicContext :: meshFillColor;
EPixel oofegGraphicContext :: edgeColor;
EPixel oofegGraphicContext :: nodeColor;
EPixel oofegGraphicContext :: bcicColor;
EPixel oofegGraphicContext :: bcForceColor;
EPixel oofegGraphicContext :: deformedElementColor;
EPixel oofegGraphicContext :: crackPatternColor;
EPixel oofegGraphicContext :: activeCrackColor;
EPixel oofegGraphicContext :: yieldPlotColors [ OOFEG_YIELD_STEPS ];
EPixel oofegGraphicContext :: standardSparseProfileColor, oofegGraphicContext :: extendedSparseProfileColor;
EPixel oofegGraphicContext :: geometryColor;
int oofegGraphicContext :: activeStep = -1;
int oofegGraphicContext :: activeStepVersion = 0;
double oofegGraphicContext :: defScale = 1.0;
double oofegGraphicContext :: zprofilescale = 0.0;
int oofegGraphicContext :: activeEigVal = 1;
int oofegGraphicContext :: activeYieldStep;
IntArray oofegGraphicContext :: matRegFilter;
dynaList< Range >oofegGraphicContext :: element_filter;
SmootherType oofegGraphicContext :: smootherType;
ScalarAlgorithmType oofegGraphicContext :: scalarAlgo = SA_ISO_SURF;
int oofegGraphicContext :: intVarDefGeoFlag = 0;
int oofegGraphicContext :: sparseProfileMode;
int oofegGraphicContext :: activeProblem = 0;
int oofegGraphicContext :: activeDomain = 0;
ScaleMode oofegGraphicContext :: smode;
double oofegGraphicContext :: emin, oofegGraphicContext :: emax;
int oofegGraphicContext :: scaleInitFlag;

InternalStateMode oofegGraphicContext :: varMode = ISM_recovered;
bool oofegGraphicContext :: staticVarsInitFlag = 0;


oofegGraphicContext :: oofegGraphicContext()
{
    isActiveFlag = false;
    plotMode = OGC_unknown;
}


void
oofegGraphicContext :: init(EngngModel *d) {
    if ( staticVarsInitFlag == 0 ) {
        BOOLEAN suc;
        int i, nmat;

        emodel = d;
        meshFillColor = ColorGetPixelFromString(oofem_tmpstr("black"), & suc);
        edgeColor = ColorGetPixelFromString(const_cast< char * >("black"), & suc);
        deformedElementColor = ColorGetPixelFromString(const_cast< char * >("BlueViolet"), & suc);
        nodeColor  = ColorGetPixelFromString(const_cast< char * >("black"), & suc);
        bcicColor  = ColorGetPixelFromString(const_cast< char * >("orange"), & suc);
        bcForceColor  = ColorGetPixelFromString(const_cast< char * >("red"), & suc);
        crackPatternColor = ColorGetPixelFromString(const_cast< char * >("gray66"), & suc);
        activeCrackColor  = ColorGetPixelFromString(const_cast< char * >("red2"), & suc);
        standardSparseProfileColor  = ColorGetPixelFromString(const_cast< char * >("blue"), & suc);
        extendedSparseProfileColor  = ColorGetPixelFromString(const_cast< char * >("red"), & suc);
        geometryColor = ColorGetPixelFromString(const_cast< char * >("yellow"), & suc);

        yieldPlotColors [ 0 ] = ColorGetPixelFromString(const_cast< char * >("pink"), & suc);
        yieldPlotColors [ 1 ] = ColorGetPixelFromString(const_cast< char * >("PaleVioletRed"), & suc);
        yieldPlotColors [ 2 ] = ColorGetPixelFromString(const_cast< char * >("maroon"), & suc);

        activeDomain = 1;

        ft = ColorCreateFringeTable();

        nmat = 0;
        for ( int id = 1; id <= d->giveNumberOfDomains(); id++ ) {
            nmat = max( nmat, d->giveDomain(id)->giveNumberOfMaterialModels() );
        }

        // ft = ColorCreateFringeTable();
        matRegFilter.resize(nmat);
        for ( i = 1; i <= nmat; i++ ) {
            matRegFilter.at(i) = 1;
        }

        staticVarsInitFlag = 1;
    }
}

oofegGraphicContext :: ~oofegGraphicContext()
{ }

EPixel
oofegGraphicContext :: GR_giveColorFromUserColorTable(EPixel *table, int tableSize, double relVal)
{
    //
    // returns color from given color table of size tableSize,
    // relVal is relative number (0..1) saying which color to return
    //
    if ( relVal > 1 ) {
        relVal = 1.;
    }

    if ( relVal < 0 ) {
        relVal = 0.;
    }

    int indx = ( int ) nearest(relVal * tableSize);
    return table [ indx - 1 ];
}

int
oofegGraphicContext :: testElementGraphicActivity(Element *e)
{
    int matFilterState = ( this->getMaterialModelFilterState( e->giveMaterial()->giveNumber() ) );
    int elemFiltState = 0;
    if ( element_filter.isEmpty() ) {
        return matFilterState;
    } else {
        dynaList< Range > :: iterator rangeIter;
        for ( rangeIter = this->element_filter.begin(); rangeIter != this->element_filter.end(); ++rangeIter ) {
            if ( ( * rangeIter ).test( e->giveNumber() ) ) {
                elemFiltState = 1;
                break;
            }
        }

        return ( matFilterState && elemFiltState );
    }
}

int
oofegGraphicContext :: getMaterialModelFilterState(int i)
{
    if ( ( i <= 0 ) || ( i > matRegFilter.giveSize() ) ) {
        return 0;
    }

    return matRegFilter.at(i);
}

void
oofegGraphicContext :: setMaterialModelFilterState(int i, int state)
{
    if ( ( i <= 0 ) || ( i > matRegFilter.giveSize() ) ) {
        return;
    }

    matRegFilter.at(i) = state;
}

void
oofegGraphicContext :: setElementFilterState(char *initString)
{
    ///@todo Anyone who uses OOFEG should check to see if StringReader can be removed in favor of the new OOFEMTXTInputRecord parser (bug ticket 24)
#if 1
    StringReader reader;

    element_filter.clear();
    reader.readRangeList(element_filter, initString, "element_filter");
#else    
    OOFEMTXTInputRecord parser(initString);
    element_filter.clear();
    parser.giveField(element_filter, IFT_Unknown, "element_filter");
#endif
}

int
oofegGraphicContext :: setActiveProblem(int a)
{
    EngngModel *slave = emodel->giveSlaveProblem(a);
    int nmat = 0;
    for ( int id = 1; id <= slave->giveNumberOfDomains(); id++ ) {
        nmat = max( nmat, slave->giveDomain(id)->giveNumberOfMaterialModels() );
    }

    // ft = ColorCreateFringeTable();
    matRegFilter.resize(nmat);
    for ( int i = 1; i <= nmat; i++ ) {
        matRegFilter.at(i) = 1;
    }

    activeProblem = a;
    return a;
}

EngngModel *
oofegGraphicContext :: getActiveProblem() {
    if ( activeProblem ) {
        return emodel->giveSlaveProblem(activeProblem);
    } else {
        return emodel;
    }
}


void
oofegGraphicContext :: updateFringeTableMinMax(double *s, int size)
{
    if ( this->getScaleMode() == SM_Autoscale ) {
        double smin = this->getScaleMin();
        double smax = this->getScaleMax();

        for ( int i = 0; i < size; i++ ) {
            if ( scaleInitFlag ) {
                smin = smax = s [ i ];
                scaleInitFlag = 0;
            } else {
                if ( smax < s [ i ] ) {
                    smax = s [ i ];
                }

                if ( smin > s [ i ] ) {
                    smin = s [ i ];
                }
            }
        }

        setScaleVals(smin, smax);
        if ( fabs(smax - smin) < 1.e-12 ) {
            smax += 1.e-12;
        }

        ColorSetupFringeTableByMinMax(this->getFringeTable(), ( FPNum ) smin, ( FPNum ) smax);
    }
}
} // end namespace oofem
#endif

