// New class created by mj
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2010   Borek Patzak                                       



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

#include "gpexportmodule.h"
#include "gausspnt.h"
#include "element.h"
#include "integrationrule.h"

#include "timestep.h"
#include "engngm.h"
#include "strreader.h"
#include "node.h"
#include "mathfem.h"
#include "oofem_limits.h"
#include "structuralms.h"
#ifndef __MAKEDEPEND
#include <vector>
#endif


namespace oofem {


GPExportModule :: GPExportModule (EngngModel* e) : ExportModule(e)
{
}


GPExportModule::~GPExportModule ()
{
}


IRResultType
GPExportModule :: initializeFrom (InputRecord* ir)
{
 ExportModule::initializeFrom (ir);
 return IRRT_OK;
}


void    
GPExportModule::doOutput (TimeStep* tStep)
{

 if (!testTimeStepOutput(tStep)) return;

 FILE* stream = this->giveOutputStream(tStep);

 fprintf(stream, "# gp DataFile Version 1.2\n");
 fprintf(stream, "# Output for time %f\n", tStep->giveTime());
  
 Domain* d  = emodel->giveDomain(1);

  int ielem, j, nelem = d -> giveNumberOfElements();
  Element* elem;
  GaussPoint* gp;
  FloatArray gcoords;
  //FloatArray intvar1,intvar2;
  FloatArray intvar3,intvar4;
  double x, y, z;
  int ncoord;
  for (ielem = 1; ielem <= nelem; ielem++) {
    elem = d->giveElement(ielem);
    IntegrationRule* iRule = elem->giveDefaultIntegrationRulePtr();
    for (j=0 ; j < iRule->getNumberOfIntegrationPoints() ; j++) {
      gp = iRule->getIntegrationPoint(j) ;
      elem->computeGlobalCoordinates (gcoords, *(gp->giveCoordinates()));
      x = gcoords.at(1);
      ncoord = gcoords.giveSize();
      y = 0.;
      if (ncoord>1)
	y = gcoords.at(2);
      z = 0.;
      if (ncoord>2)
	z = gcoords.at(3);
      //elem->giveIPValue (intvar1, gp, IST_DissWorkDensity, tStep);
      //elem->giveIPValue (intvar2, gp, IST_StressWorkDensity, tStep);
      elem->giveIPValue (intvar3, gp, IST_DamageTensor, tStep);
      elem->giveIPValue (intvar4, gp, IST_MaxEquivalentStrainLevel, tStep);

  // integration weight (contributing area or volume)
  //    double weight = gp->giveElement()->computeVolumeAround(gp);
  // multiply by weight if you want the contribution to the total
  // here we work just with densities
  //    double stressWork = intvar2.at(1) * weight;
      //double stressWork = intvar2.at(1);
      //double dissWork = intvar1.at(1);
      //double freeEnergy = stressWork-dissWork;
      // write data only for Gauss points with nonzero dissipation 
      //if (dissWork > 1.e-6 * freeEnergy)
      {
	fprintf (stream, "%d %d %f %f %f ",elem->giveNumber(),j+1,x,y,z);
	//fprintf (stream, "%f %f %f ",dissWork,freeEnergy,stressWork);
	// for CST elements write also nodal coordinates
	/*
	int nnode = elem->giveNumberOfNodes();
	if (nnode==3){
	  for (int inod=1; inod<=3; inod++)
	    fprintf (stream, "%f %f ",elem->giveNode(inod)->giveCoordinate(1),elem->giveNode(inod)->giveCoordinate(2));
	}
	*/
	//StructuralMaterialStatus* stat = (StructuralMaterialStatus*) gp->giveMaterialStatus();
	//double eps = stat->giveStrainVector().at(1);
	//fprintf (stream, "%f %f %f ",eps,intvar3.at(1),intvar4.at(1));
	double damage = 0., kappa = 0.;
	if (intvar3.giveSize()>0)
	  damage = intvar3.at(1);
	if (intvar4.giveSize()>0)
	  kappa = intvar4.at(1);
	fprintf (stream, "%f %f ",damage,kappa);
	fprintf (stream, "\n");
	// we have written: element_number, gp_number, x_coordinate, y_coordinate, z_coordinate, damage, kappa
      }
    }
  }
 
 fclose (stream);
}

void
GPExportModule::initialize ()
{}


void
GPExportModule::terminate ()
{}


FILE* 
GPExportModule::giveOutputStream (TimeStep* tStep) 
{
 char baseFileName[MAX_FILENAME_LENGTH];
 char fileName[MAX_FILENAME_LENGTH];
 FILE* answer;

 emodel->giveOutputBaseFileName (baseFileName, MAX_FILENAME_LENGTH);
 sprintf (fileName, "%s.%d.gp", baseFileName, tStep->giveNumber());
 if ((answer = fopen (fileName,"w")) == NULL) {
   OOFEM_ERROR2 ("GPExportModule::giveOutputStream: failed to open file %s", fileName);
 }
 return answer;

}

} // namespace oofem
