/* $Header: $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2005   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

#include "levelsetpcs.h"
#include "mathfem.h"
#include "timestep.h"
#include "node.h"
#include "datastream.h"
#include "conTable.h"
#include "spatiallocalizer.h"
#include "geotoolbox.h"
#include "error.h"

void
LevelSetPCS::initialize()
{
  if (0) {
    if (previousLevelSetValues.giveSize() != domain->giveNumberOfDofManagers()) 
      OOFEM_ERROR ("LevelSetPCS::initialize size of levelSetValues does not match number of dof managers");
  } else {
  
    if (initialRefMatFlag) {
      int nnodes = domain->giveNumberOfDofManagers();
      previousLevelSetValues.resize(nnodes);
      for (int i=1; i<=nnodes; i++) {
        previousLevelSetValues.at(i) = (-1.0)*initialRefMatVol.pointDistance(domain->giveNode(i)->giveCoordinate(1),
                                                                             domain->giveNode(i)->giveCoordinate(2));
      }
    }

    /*
      int nnodes = domain->giveNumberOfDofManagers();
      previousLevelSetValues.resize(nnodes);
      FloatArray center(3); center.at(1)=1.5; center.at(2)=2.0; center.at(3)=0.0;
      //printf ("\n");
      for (int i=1; i<=nnodes; i++) {
      previousLevelSetValues.at(i) = 1.0 - domain->giveNode(i)->giveCoordinates()->distance(center);
      }
    */
    levelSetValues = previousLevelSetValues;
    previousLevelSetValues.printYourself();
  }
}

IRResultType 
LevelSetPCS::initializeFrom (InputRecord* ir) 
{
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;   
  
  IR_GIVE_OPTIONAL_FIELD (ir, previousLevelSetValues, IFT_LSPCS_levelSetValues, "levelset");
  if (!previousLevelSetValues.giveSize()) {
    FloatArray refmatpoly_x, refmatpoly_y;
    IR_GIVE_OPTIONAL_FIELD (ir, refmatpoly_x, IFT_LSPCS_refmatpoly_x, "refmatpolyx");
    IR_GIVE_OPTIONAL_FIELD (ir, refmatpoly_y, IFT_LSPCS_refmatpoly_y, "refmatpolyy");
    int nvert = refmatpoly_x.giveSize();
    if (nvert) {
      Vertex v; 
      for (int i=1; i<=nvert; i++) {
        // create polygonal representation
        v.setCoords (refmatpoly_x.at(i), refmatpoly_y.at(i));
        initialRefMatVol.addVertex (v);
      }
      // close polygon (add first vertex at the end
      v.setCoords (refmatpoly_x.at(1), refmatpoly_y.at(1));
      initialRefMatVol.addVertex (v);
    }
    initialRefMatFlag = true;
  }
  return IRRT_OK;
}


void
LevelSetPCS::updatePosition (TimeStep* atTime)
{
  int i,j,l,inodes,inode,nsd = 2;
  int ndofman = domain->giveNumberOfDofManagers();
  bool twostage = true;

  double help,dt,volume, gfi_norm;

  FloatArray fs(ndofman), w(ndofman);
  FloatMatrix dN;
  FloatArray fi(4), gfi(nsd), n(nsd), k(4), dfii(4), alpha(4), un;
  LevelSetPCSElementInterface* interface;
  Element* ielem;
  ConnectivityTable* contable = domain->giveConnectivityTable();

  // needed for multistep update
  FloatArray ls_n; 
  int __step = 0, __nstep = 10;

  levelSetValues = previousLevelSetValues;
  dt=atTime->giveTimeIncrement()/ __nstep;

  do {
    ls_n = levelSetValues;

    this->pcs_stage1 (fs, w, atTime, PCS_levelSetUpdate);
    
    // update level set values
    for (inode=1; inode<=ndofman; inode++) {
      if (w.at(inode) > 0.0) {
        // single stage integration
        levelSetValues.at(inode)= ls_n.at(inode)-dt*fs.at(inode)/w.at(inode);
      } else {
        // -------------------------
        // inflow into boundary node
        // -------------------------
        const IntArray* elems; IntArray mask(2);
        double v;
        // get velocity in inode
        mask.at(1) = V_u; mask.at(2) = V_v;
        domain->giveDofManager(inode)->giveUnknownVector (un, mask, EID_MomentumBalance,VM_Total,atTime->givePreviousStep());
        elems = contable->giveDofManConnectivityArray(inode);
        // loop over shared elements
        volume=0.0; help=0.0;
        for (l=1;l<=elems->giveSize(); l++) {
          // get element level set gradient
          ielem = domain->giveElement(elems->at(l));
          inodes = ielem->giveNumberOfNodes();
          interface = (LevelSetPCSElementInterface*) 
            ielem->giveInterface(LevelSetPCSElementInterfaceType);
          
          if (interface) {
            
            interface->LS_PCS_computedN(dN);
            // assemble element vector with  level set values
            for (i=1; i<= inodes; i++) fi.at(i) = ls_n.at(ielem->giveDofManagerNumber(i));
            // compute gradient of level set
            for (j=1; j<=nsd; j++) {
              gfi.at(j) = 0.0;
              for (i=1; i<=inodes; i++) gfi.at(j)+=dN.at(i,j)*fi.at(i);
            }
            volume += (v = interface->LS_PCS_computeVolume());
            gfi_norm = sqrt(dotProduct(gfi,gfi,nsd));
            if (gfi_norm > 1.e-6) 
              help += dotProduct(un,gfi,nsd)*v/sqrt(dotProduct(gfi,gfi,nsd));
          }
        } // end loop over shared nodes
        levelSetValues.at(inode)= ls_n.at(inode) - dt*help/volume;
      }
    } // end loop over nodes
    
    if (twostage) {
      this->pcs_stage1 (fs, w, atTime, PCS_levelSetUpdate);
      
      for (inode=1; inode<=ndofman; inode++) {
        if (w.at(inode) > 0.0) {
          //two stage integration
          // update
          levelSetValues.at(inode) = 0.5*(ls_n.at(inode)+levelSetValues.at(inode))-
            0.5*dt*fs.at(inode)/w.at(inode);
        }
      }
    }
    printf (".");
  } while (++__step <  __nstep);
  printf ("\n");
  // print level set values to stdout (debug only)
  /*
    printf ("Node: Level Set Value\n");
    for (inode=1; inode<=ndofman; inode++) {
    printf ("%5d %le\n",inode, levelSetValues.at(inode));
    }
  */
  // redistance
  this->redistance (atTime);

}


double 
LevelSetPCS::computeCriticalTimeStep (TimeStep* tStep)
{
  return 1.e6;
}


void
LevelSetPCS::giveMaterialMixtureAt (FloatArray& answer, FloatArray& position)
{
  double ls;
  int i;
  FloatArray N(3);
  answer.resize(2);
  
  Element* elem = domain->giveSpatialLocalizer()->giveElementContainingPoint (position);
  LevelSetPCSElementInterface* interface = (LevelSetPCSElementInterface*) elem->giveInterface(LevelSetPCSElementInterfaceType);
  if (interface) {
    if (elem->computeLocalCoordinates(N, position)) {
      int inodes = elem->giveNumberOfNodes();
      for (ls=0.0, i=1; i<=inodes;i++) {
        ls += N.at(i) * levelSetValues.at(elem->giveDofManagerNumber(i));
      } 
      if (ls > 0.0) {
        answer.at(1)=1.0; answer.at(2) = 0.0;
      } else {
        answer.at(1)=0.0; answer.at(2) = 1.0;
      }
    } else {
      OOFEM_ERROR ("LevelSetPCS::giveMaterialMixtureAt: computeLocalCoordinates failed");
    }
  } else {
    answer.at(1) = 1.0;
    answer.at(2) = 0.0;
  }
  
}

void
LevelSetPCS::giveElementMaterialMixture (FloatArray& answer, int ie)
{
  Element *ielem = domain->giveElement(ie);
  int i, inodes = ielem->giveNumberOfNodes();
  LevelSetPCSElementInterface* interface = (LevelSetPCSElementInterface*) ielem->giveInterface(LevelSetPCSElementInterfaceType);
  FloatArray fi(inodes);

  for (i=1; i<= inodes; i++) fi.at(i) = levelSetValues.at(ielem->giveDofManagerNumber(i));
  interface->LS_PCS_computeVOFFractions (answer, fi);
}


void
LevelSetPCS::redistance (TimeStep* atTime)  
{
  int nite=0, inode, nsd = 2, l,i,j,inodes;
  int ndofman = domain->giveNumberOfDofManagers();
  bool twostage = false;
  double dt, c, cm, v, volume, help, gfi_norm;

  FloatArray fs(ndofman), w(ndofman), b;
  FloatMatrix dN;
  FloatArray fi(4), gfi(nsd), n(nsd);
  ConnectivityTable* contable = domain->giveConnectivityTable();
  Element* ielem;
  LevelSetPCSElementInterface* interface;

  dt=atTime->giveTimeIncrement();

  do {

    b = levelSetValues;
    pcs_stage1 (fs, w, atTime, PCS_levelSetRedistance);    
    
    // update level set values
    // single stage integration
    cm = 0.0;
    for (inode=1; inode<=ndofman; inode++) {
      if (w.at(inode) > 0.0) {
        c = dt*fs.at(inode)/w.at(inode);
        levelSetValues.at(inode) = b.at(inode) - c;
        cm = max (cm, fabs(c));
      }
    } 

    if (1) {
      double inodex, inodey, _nodex, _nodey;
      int _node, _count;
      for (inode=1; inode<=ndofman; inode++) {
        if (w.at(inode) == 0.0) {
          // -------------------------
          // inflow into boundary node
          // -------------------------
          
          const IntArray* elems;
          elems = contable->giveDofManConnectivityArray(inode);
          // loop over shared elements
          volume=0.0; help=0.0; _count=0;
          gfi.zero();
          
          for (l=1;l<=elems->giveSize(); l++) {
            // get element level set gradient
            ielem = domain->giveElement(elems->at(l));
            inodes = ielem->giveNumberOfNodes();
            interface = (LevelSetPCSElementInterface*) 
              ielem->giveInterface(LevelSetPCSElementInterfaceType);
            
            if (interface) {
              
              interface->LS_PCS_computedN(dN);
              // assemble element vector with  level set values
              for (i=1; i<= inodes; i++) fi.at(i) = b.at(ielem->giveDofManagerNumber(i));
              // compute average gradient of level set
              v = interface->LS_PCS_computeVolume();
              for (j=1; j<=nsd; j++) {
                for (i=1; i<=inodes; i++) gfi.at(j)+=dN.at(i,j)*fi.at(i);
              }
              gfi_norm = sqrt(dotProduct(gfi,gfi,nsd));
              volume += v;
              inodex=domain->giveNode(inode)->giveCoordinate(1);
              inodey=domain->giveNode(inode)->giveCoordinate(2);
              for (i=1; i<= inodes; i++) {
                _node = ielem->giveDofManagerNumber(i);
                if (_node == inode) continue;
                if (w.at(_node) == 0.0) continue;
                _node = ielem->giveDofManagerNumber(i);
                _nodex=domain->giveNode(_node)->giveCoordinate(1);
                _nodey=domain->giveNode(_node)->giveCoordinate(2);

                help += (levelSetValues.at(_node)+((inodex-_nodex)*gfi.at(1)+(inodey-_nodey)*gfi.at(2))/gfi_norm);
		_count++;
              }
            }
	  } // end loop over shared elements of inflow node
	  levelSetValues.at(inode) = help/_count;
	  cm = max (cm, levelSetValues.at(inode)-b.at(inode));
          
          
        }
      }
    } else {
      // -------------------------
      // inflow into boundary node
      // -------------------------
      const IntArray* elems;
      elems = contable->giveDofManConnectivityArray(inode);
      // loop over shared elements
      volume=0.0; help=0.0;
      gfi.zero();
      
      for (l=1;l<=elems->giveSize(); l++) {
        // get element level set gradient
        ielem = domain->giveElement(elems->at(l));
        inodes = ielem->giveNumberOfNodes();
        interface = (LevelSetPCSElementInterface*) 
          ielem->giveInterface(LevelSetPCSElementInterfaceType);
        
        if (interface) {
          
          interface->LS_PCS_computedN(dN);
          // assemble element vector with  level set values
          for (i=1; i<= inodes; i++) fi.at(i) = b.at(ielem->giveDofManagerNumber(i));
          // compute average gradient of level set
          v = interface->LS_PCS_computeVolume();
          for (j=1; j<=nsd; j++) {
            for (i=1; i<=inodes; i++) gfi.at(j)+=v*dN.at(i,j)*fi.at(i);
          }
          volume += v;
        }
      } // end loop over shared nodes
      gfi.times(1.0/volume); 
      gfi_norm = sqrt(dotProduct(gfi,gfi,nsd));
      if (gfi_norm > 0.5) {
        help = sgn(levelSetValues.at(inode))*(1.0-gfi_norm);
        levelSetValues.at(inode)= b.at(inode) - dt*help;
        cm = max (cm, dt*help);
      }
    }
    
    if (twostage) {
      this->pcs_stage1 (fs, w, atTime, PCS_levelSetRedistance);
      cm = 0.0;
      for (inode=1; inode<=ndofman; inode++) {
        if (w.at(inode) > 0.0) {
          //two stage integration
          // update
          levelSetValues.at(inode) = 0.5*(b.at(inode)+levelSetValues.at(inode))-
            0.5*dt*fs.at(inode)/w.at(inode);
          cm = max (cm, fabs (levelSetValues.at(inode)-b.at(inode)));
        }
      }
    }
    printf ("cm: %le\n", cm);
  } while ((cm > 1.e-4) && (++nite < 10));
  
}


void
LevelSetPCS::pcs_stage1 (FloatArray& fs, FloatArray& w, TimeStep* atTime, PCSEqType t)
{
  int i, j, l, inodes, nsd = 2;
  int ie, ndofman = domain->giveNumberOfDofManagers(), nelem   = domain->giveNumberOfElements();
  double alpha, dfi, help, sumkn, F, f, volume, gfi_norm;
  FloatMatrix dN;
  FloatArray gfi(nsd), fi(4), n(nsd), k(4), dfii(4);
  Element* ielem;
  LevelSetPCSElementInterface* interface;

  fs.resize(ndofman); w.resize(ndofman);
  fs.zero(); w.zero();

  // loop over elements
  for (ie=1; ie<=nelem; ie++) {
    ielem = domain->giveElement(ie);
    inodes = ielem->giveNumberOfNodes();
    interface = (LevelSetPCSElementInterface*) 
      ielem->giveInterface(LevelSetPCSElementInterfaceType);

    if (interface) {

      F = this->evalElemFContribution (t, ie, atTime);
      f = this->evalElemfContribution (t, ie, atTime);

      interface->LS_PCS_computedN(dN);
      volume = interface->LS_PCS_computeVolume();
      
      // assemble element vector with  level set values
      for (i=1; i<= inodes; i++) fi.at(i) = levelSetValues.at(ielem->giveDofManagerNumber(i));
      // compute gradient of level set
      for (j=1; j<=nsd; j++) {
        gfi.at(j) = 0.0;
        for (i=1; i<=inodes; i++) gfi.at(j)+=dN.at(i,j)*fi.at(i);
      }
      // compute ki
      for (i=1; i<=inodes; i++) {
        // eval size of gfi
        gfi_norm = sqrt(dotProduct(gfi, gfi, nsd));
        if (gfi_norm > 1.e-3) {
          // evaluate i-th normal (correcponding to side opposite to i-th vertex)
          for (j=1; j<=nsd; j++) n.at(j)=dN.at(i,j)*2.0*volume;
          k.at(i) = F*dotProduct(gfi, n, nsd) / (2.0*gfi_norm);
        } else k.at(i) = 0.0;
      }
      dfi = dotProduct (k,fi,inodes);
      for (i=1; i<=inodes; i++) {
        help = 0.0;sumkn = 0.0;
        for (l=1; l<=inodes; l++) {
          help += negbra(k.at(l))*(fi.at(i)-fi.at(l));
          sumkn+= negbra(k.at(l));
        }
        dfii.at(i)=macbra(k.at(i))*help/sumkn;
      }
      //compute alpha_i
      for (help=0.0, l=1; l<=inodes; l++) help+=max(0.0, dfii.at(l)/dfi);
      for (i=1; i<=inodes; i++) {
        if (fabs(help) > 0.0) {
          alpha=max(0.0,dfii.at(i)/dfi)/help;
          fs.at(ielem->giveDofManagerNumber(i))+=alpha*(dfi-f*volume);
          w.at(ielem->giveDofManagerNumber(i)) += alpha*volume;
        }
      }
    } else {
      OOFEM_ERROR2 ("LevelSetPCS::updatePosition: element %d does not implement LevelSetPCSElementInterfaceType",ie);
    }
  }// end loop over elements

}


double
LevelSetPCS::evalElemFContribution(PCSEqType t, int ie, TimeStep* atTime)
{
  LevelSetPCSElementInterface* interface = (LevelSetPCSElementInterface*) 
    domain->giveElement(ie)->giveInterface(LevelSetPCSElementInterfaceType);
  if (t == PCS_levelSetUpdate) 
    return interface->LS_PCS_computeF (this, atTime);
  else if (t == PCS_levelSetRedistance) 
    return interface->LS_PCS_computeS (this, atTime);
  return 0.0;
}

double
LevelSetPCS::evalElemfContribution(PCSEqType t, int ie, TimeStep* atTime)
{
  LevelSetPCSElementInterface* interface = (LevelSetPCSElementInterface*) 
    domain->giveElement(ie)->giveInterface(LevelSetPCSElementInterfaceType);
  if (t == PCS_levelSetUpdate) 
    return 0.0;
  else if (t == PCS_levelSetRedistance) 
    return interface->LS_PCS_computeS (this, atTime);
  return 0.0;
}
