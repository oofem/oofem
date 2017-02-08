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

#include "../sm/Elements/Beams/beambaseelement.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "bctracker.h"

#include "bodyload.h"
#include "boundaryload.h"

namespace oofem {


BeamBaseElement :: BeamBaseElement (int n, Domain *aDomain) : StructuralElement(n, aDomain)
{}

BeamBaseElement :: ~BeamBaseElement()
{}

void
BeamBaseElement :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// Why is this function taken separately ?
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further subtract part corresponding to non-nodal loading.
{
    FloatArray helpLoadVector(1);
    answer.clear();

    // loop over body load array first
    int nBodyLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nBodyLoads; i++ ) {
        int id = bodyLoadArray.at(i);
        Load *load = domain->giveLoad(id);
        bcGeomType ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            this->computeBodyLoadVectorAt(helpLoadVector, load, tStep, mode);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            if ( load->giveBCValType() != TemperatureBVT && load->giveBCValType() != EigenstrainBVT ) {
                // temperature and eigenstrain is handled separately at computeLoadVectorAt subroutine
                OOFEM_ERROR("body load %d is of unsupported type (%d)", id, ltype);
            }
        }
    }

    // loop over boundary load array
    int nBoundaryLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nBoundaryLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        Load *load = domain->giveLoad(n);
	BoundaryLoad* bLoad;
	if ((bLoad = dynamic_cast<BoundaryLoad*> (load))) {
	  bcGeomType ltype = load->giveBCGeoType();
	  if ( ltype == EdgeLoadBGT ) {
	    this->computeBoundaryEdgeLoadVector(helpLoadVector, bLoad, id, ExternalForcesVector, mode, tStep, false);
            if ( helpLoadVector.giveSize() ) {
	      answer.add(helpLoadVector);
            }
	  } else if ( ltype == SurfaceLoadBGT ) {
	    this->computeBoundarySurfaceLoadVector(helpLoadVector, bLoad, id, ExternalForcesVector, mode, tStep, false);
            if ( helpLoadVector.giveSize() ) {
	      answer.add(helpLoadVector);
            }
	  } else if ( ltype == PointLoadBGT ) {
            // id not used
	    this->computePointLoadVectorAt(helpLoadVector, load, tStep, mode, false);
            if ( helpLoadVector.giveSize() ) {
	      answer.add(helpLoadVector);
            }
	  } else {
            OOFEM_ERROR("boundary load %d is of unsupported type (%d)", id, ltype);
	  }
	}
    }


    // add exact end forces due to nonnodal loading applied indirectly (via sets)
    BCTracker *bct = this->domain->giveBCTracker();
    BCTracker::entryListType bcList = bct->getElementRecords(this->number);
    FloatArray help;
    
    for (BCTracker::entryListType::iterator it = bcList.begin(); it != bcList.end(); ++it) {
      GeneralBoundaryCondition *bc = this->domain->giveBc((*it).bcNumber);
      BodyLoad *bodyLoad;
      BoundaryLoad *boundaryLoad;
      if (bc->isImposed(tStep)) {
        if ((bodyLoad = dynamic_cast<BodyLoad*>(bc))) { // body load
          this->computeBodyLoadVectorAt(help,bodyLoad, tStep, VM_Total); // this one is local
          answer.add(help);
        } else if ((boundaryLoad = dynamic_cast<BoundaryLoad*>(bc))) {
          // compute Boundary Edge load vector in GLOBAL CS !!!!!!!
          this->computeBoundaryEdgeLoadVector(help, boundaryLoad, (*it).boundaryId,
					      ExternalForcesVector, VM_Total, tStep, false);
          // get it transformed back to local c.s.
          // this->computeGtoLRotationMatrix(t);
          // help.rotatedWith(t, 'n');
          answer.add(help);
        }
      }
    }
}

} // end namespace oofem
