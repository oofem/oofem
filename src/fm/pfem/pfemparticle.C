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


#include "pfemparticle.h"
#include "timestep.h"
#include "classfactory.h"
#include "fmode.h"
#include "domain.h"
#include "engngm.h"
#include "mathfem.h"


namespace oofem {
REGISTER_DofManager(PFEMParticle);
/**
 * Constructor. Creates a particle with number n, belonging to aDomain.
 */
  PFEMParticle :: PFEMParticle(int n, Domain *aDomain) : Node(n, aDomain), freeFlag(true), alphaShapeFlag(false), activeFlag(true)
{ }

/**
 * Gets from the source line from the data file all the data of the receiver.
 */
IRResultType
PFEMParticle :: initializeFrom(InputRecord *ir)
{
    return Node :: initializeFrom(ir);
}

/**
 * Checks internal data consistency in node.
 */
int
PFEMParticle :: checkConsistency()
{
    int result = 1;

    result = result && Node :: checkConsistency();

    return result;
}

void
PFEMParticle :: updateYourself(TimeStep *tStep)
{
    // TODO the implementation of free particle movement should be to updateYourself(),
    // now it is in pfem::solveYourselfAt()
    Node :: updateYourself(tStep);
}

void
PFEMParticle :: printOutputAt(FILE *stream, TimeStep *stepN)
{
    DofManager :: printOutputAt(stream, stepN);
}

#ifdef __OOFEG
void PFEMParticle :: drawScalar(oofegGraphicContext &gc)
{
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    if ( gc.giveIntVarType() == IST_Pressure ) {
        WCRec p [ 1 ]; /* point */
        p [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveCoordinate(3);

        int dofindx;
        double pressVal;
        if ( ( dofindx = this->findDofWithDofId(P_f) ) ) {
            pressVal = this->giveDof(dofindx)->giveUnknown(VM_Total, tStep);
        }

        EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
        EASValsSetColor( ColorFringeRangeToColor( ColorFringeValueToRange(gc.getFringeTable(), pressVal) ) );
        EASValsSetMType(FILLED_CIRCLE_MARKER);

        EASValsSetMSize(8);
        go = CreateMarker3D(p);
        EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}
#endif
} // end namespace oofem
