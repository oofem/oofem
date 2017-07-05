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
 *  License as published by the Free Software Foundation; eitherc
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

#include "../sm/Quasicontinuum/geometrygenerator.h"
#include <stdlib.h>
#include <time.h>


namespace oofem {
GeometryGenerator :: GeometryGenerator()

// Constructor.
{ }

GeometryGenerator :: ~GeometryGenerator()
// Destructor
{ }

IRResultType
GeometryGenerator :: initializeParticleGenerator(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, maxNumOfParticles, _IFT_GeometryGenerator_numOfParticles);
    IR_GIVE_FIELD(ir, maxNumOfIterations, _IFT_GeometryGenerator_numOfIterations);
    IR_GIVE_FIELD(ir, maxNumOfItOnePar, _IFT_GeometryGenerator_numOfItOnePar);
    IR_GIVE_FIELD(ir, ParticleRadius, _IFT_GeometryGenerator_particleRadius);


    //...

    // check input format
#ifdef DEBUG
    /*
     * if ( ) {
     * OOFEM_ERROR("invalid format of ... ??? ");
     * }
     */
#endif

    return IRRT_OK;
}

void
GeometryGenerator :: generateParticles()
{
    /*
     * int NumOfParticles,  NumOfIterations, NumOfItOnePar = 1;
     * std::srand(std::time(0)); // randomize
     * FloatArray coords(3);
     * // minimal and maximal coordinates
     *  double x1= ??? ... x2 y1 y2 z1 z2
     * //nog=1; tnog=1; % pocet generaci je na zacatku
     * while (NumOfParticles<maxNumOfParticles) {
     * coords.at(1) = ((double) std::rand() / (RAND_MAX))*(x2-x1) + x1;
     * coords.at(1) = ((double) std::rand() / (RAND_MAX))*(x2-x1) + x1;
     * coords.at(1) = ((double) std::rand() / (RAND_MAX))*(x2-x1) + x1;
     *
     * if ( CheckDistances(ParticleRadius, coords, n) ) {
     *  // add this particle
     *
     * }
     *
     * if () {
     * // reached maximal number of iteration of one particles
     *
     * }
     * // reached maximal number of all iteration
     * if () {
     *
     * }
     *
     * }
     *
     *
     *    x = x1 + (x2-x1).*rand(1,1);
     *    y = y1 + (y2-y1).*rand(1,1);
     *    z = z1 + (z2-z1).*rand(1,1);
     *
     *    % check zaroven upravuje geometrii do zuzeni
     *    OK = CheckDistances(x, y, z, R, n, particles, quadtree, x1, x2, y1, y2, z1, z2, nod);
     *    if(OK)
     *        particles(n+1,1) = x;
     *        particles(n+1,2) = y;
     *        particles(n+1,3) = z;
     *        quadtree  = addParticleInQuadTree(x, y, n+1, quadtree, x1, x2, y1, y2, nod );
     *        //n=n+1;
     *        //nog=1; % castice se pridala - lokalni citac generaci se nuluje
     *        NumOfParticles++;
     *        nog=1;
     *    end
     *    nog=nog+1; tnog=tnog+1;
     *    if (nog>=maxnog) % prekrocen pocet generaci pro prijeti jedne castice
     *        fprintf('v %d krocich vygenerovano:\n%d z %d castic\n', tnog, n, nop);
     *        fprintf('dosazen maximalni pocet %d generaci pro prijeti jedne castice\n', maxnog);
     *        break
     *
     *    elseif (tnog>=tmaxnog) % prekrocet celkovy pocet generaci
     *        fprintf('v %d krocich vygenerovano:\n%d z %d castic\n', tnog, n, nop);
     *        fprintf('dosazen maximalni celkovy pocet %d generaci\n', tmaxnog);
     *        break
     *    elseif (n>=nop) % dosazen pozadovany pocet castic
     *        fprintf('v %d krocich vygenerovano vsech %d pozadovanych castic\n', tnog, n);
     *    end
     *    % kontrolni vypis
     *    if (floor(tnog/1000)*1000==tnog) % vypisuj kazdou 1000.
     *        fprintf('  aktualni pocet generaci: %d, vygenerovano %d castic\n', tnog, n);
     *    end
     */
}

bool
GeometryGenerator :: CheckDistances(double R, FloatArray coords, int n)
// retunrs true if the distance of particle with coords from first n particles is greater or equal than R
{
    if ( n == 0 ) {
        return true;
    }
    if ( (unsigned int) n > Particles.size() ) {
        n = Particles.size();
    }

    double R2 = R * R;
    for ( int i = 1; i <= n; i++ ) {
        double distance2 = coords.distance_square(Particles [ i - 1 ]);
        if ( distance2 < R2 ) {
            return false;
        }
    }

    return false;
}


void
GeometryGenerator :: loadParticles()
{}

IRResultType
GeometryGenerator :: initializeLinkGenerator(InputRecord *ir)
{
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    return IRRT_OK;

}

void
GeometryGenerator :: generateLinks()
{}


void
GeometryGenerator :: loadLinks()
{}
} // end namespace oofem
