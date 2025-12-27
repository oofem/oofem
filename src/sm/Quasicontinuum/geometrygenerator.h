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

#ifndef geometrygenerator_h
#define geometrygenerator_h

#include "floatarray.h"
#include "element.h"
#include "node.h"

#define _IFT_GeometryGenerator_numOfParticles "gg_nop"
#define _IFT_GeometryGenerator_numOfIterations "gg_noi"
#define _IFT_GeometryGenerator_numOfItOnePar "gg_noiop"
#define _IFT_GeometryGenerator_particleRadius "gg_rp"
//...

namespace oofem {

/**
 * Generate random geometry of particles and links for CQ simulation.
 * 
 */
class GeometryGenerator 
{
protected:
    // global
    int nop; // number of particles
    int nol; // number of links
    std::vector<FloatArray> Particles;
    std::vector<IntArray> Links;

    // particleGenerator
    double ParticleRadius;  // minimal distance of two particles
    int maxNumOfParticles;  // maximal number of generated particles
    int maxNumOfIterations; // maximal number of iterations during generation
    int maxNumOfItOnePar;   // maximal number of generation of one particle

public:
    GeometryGenerator();
    virtual ~GeometryGenerator();

    void initializeParticleGenerator(InputRecord &ir);
    void generateParticles();
    void loadParticles();

    bool CheckDistances(double R, FloatArray coords, int n);


    void initializeLinkGenerator(InputRecord &ir);
    void generateLinks();
    void loadLinks();

    virtual const char *giveClassName() const { return "QCFullsolveddomain"; }
};
} // end namespace oofem
#endif // geometrygenerator_h
