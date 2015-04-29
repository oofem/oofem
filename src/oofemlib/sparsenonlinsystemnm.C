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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "sparsenonlinsystemnm.h"
#include "error.h"
#include "domain.h"
#include "node.h"
#include "unknownnumberingscheme.h"

namespace oofem {

IRResultType
SparseNonLinearSystemNM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    randPertAmplitude = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, randPertAmplitude, _IFT_NonLinearStatic_randPertAmplitude);
    if ( randPertAmplitude < 0. ) {
        OOFEM_WARNING("Random pertubation amplitude can not be negative");
        return IRRT_BAD_FORMAT;
    }
    randSeed = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, randSeed, _IFT_NonLinearStatic_randSeed);
 
    // optional parameters related to perturbations of the initial guess (first iteration)
    igp_PertDmanDofSrcArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, igp_PertDmanDofSrcArray, _IFT_NonLinearStatic_pert);
    igp_PertWeightArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, igp_PertWeightArray, _IFT_NonLinearStatic_pertw);
    if ( igp_PertDmanDofSrcArray.giveSize() ) {
        if ( ( igp_PertDmanDofSrcArray.giveSize() % 2 ) != 0 ) {
            OOFEM_WARNING("Pert map size must be an even number, it contains pairs <node, nodeDof>");
            return IRRT_BAD_FORMAT;
        }
        int nsize = igp_PertDmanDofSrcArray.giveSize() / 2;
        if ( igp_PertWeightArray.giveSize() != nsize ) {
            OOFEM_WARNING("Pert map size and weight array size mismatch");
            return IRRT_BAD_FORMAT;
        }
        pert_init_needed = true;
    } else {
        pert_init_needed = false;
    }
    return IRRT_OK;
}

void 
SparseNonLinearSystemNM :: convertPertMap()
{
    EModelDefaultEquationNumbering dn;
    int count = 0, ndofman = this -> domain -> giveNumberOfDofManagers();
    int size = igp_PertDmanDofSrcArray.giveSize() / 2;
    igp_Map.resize(size);
    igp_Weight.resize(size);

    for ( int j = 1; j <= ndofman; j++ ) {
        int jglobnum = this->domain->giveNode(j)->giveLabel();
        for ( int i = 1; i <= size; i++ ) {
            int inode = igp_PertDmanDofSrcArray.at(2 * i - 1);
            int idof  = igp_PertDmanDofSrcArray.at(2 * i);
            if ( inode == jglobnum ) {
                igp_Map.at(++count) = this -> domain ->giveNode(j)->giveDofWithID(idof)->giveEquationNumber(dn);
                igp_Weight.at(count) = igp_PertWeightArray.at(i);
                continue;
            }
        }
    }
}
    
void 
SparseNonLinearSystemNM :: applyPerturbation(FloatArray* displacement)
{
    int nsize;

    // First type of perturbation - random perturbation with uniform probability density 
    // over the interval (-randPertAmplitude,randPertAmplitude) applied on each unknown.
    if (randPertAmplitude > 0.) {
        nsize = displacement -> giveSize();
        srand(randSeed);
        for ( int i = 1; i <= nsize; i++ ) {
            double pert = randPertAmplitude * (2. * rand() / RAND_MAX - 1.);
            displacement->at(i) += pert;
        }
    }

    // Second type of perturbation - only selected unknowns perturbed by the amount specified by the user.
    if ( pert_init_needed ) {
        convertPertMap();
        pert_init_needed = false;
    }
      
    nsize = igp_Weight.giveSize();
    for ( int i = 1; i <= nsize; i++ ) {
        int iDof = igp_Map.at(i);
        double w = igp_Weight.at(i);
        displacement->at(iDof) += w;
    }
}

} // end namespace oofem
