/* $Header: /home/cvs/bp/oofem/oofemlib/src/sloanlevelstruct.h,v 1.5 2003/04/06 14:08:26 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



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

//   *****************************************************
//   *** CLASS SLOAN PROFILE OPTIMIZER LEVEL STRUCTURE ***
//   *****************************************************


/* Modified and optimized by: Borek Patzak */
/* Author: Milan Jirasek */


#ifndef sloanlevelstruct_h
#include "alist.h"
#include "intarray.h"

class SloanGraph;

/**
 Class representing level structure.
 This is partitioning of the nodes such that each node is assigned 
 to one of the levels in accordance with its distance from a specified root node.
*/
class SloanLevelStructure
{
private: 
  /// Reference to corresponding graph
  SloanGraph*   Graph;
 /// Root node of level structure
  int Root;
 /// End node of root structure
  int End;
 /// Data representation of structure: List of arrays, one array for each level.
  AList<IntArray> Structure;
 /// Depth of structure defined as number of levels
  int  Depth;
 /// Width of structure defined as max number of nodes in all levels.
  int  Width;

public:
 /// Constructor. Creates new level structure assignet to graph, with root being the root node.
 SloanLevelStructure(SloanGraph* graph, int root) : Structure(0)
  {Graph=graph; Root=root; End=0; Depth=Width=0; }
 /// Destructor
 ~SloanLevelStructure();

 /// Destroys all levels
  void  destroyLevels();
 /**
  Builds the level structure. The limitWidth parameter allows receiver build-up phase
  to be aborted during the assembly, when with of sone level is greater than given value. 
  Default value for limitWidth is -1 meaning no width limit.
  If assembly aborted, the destroyLevels() method is called.
  @return zero if assembly aborted due to width limit, positive value otherwise.
  */
  int   formYourself(int limitWidth = -1);

 /// Returns the depth of receiver
  int  giveDepth() {if (!Depth) computeDepth(); return Depth;}
 /// Returns the width of receiver
  int  giveWidth() {if (!Width) computeWidth(); return Width;}
 /// Returns the i-th level of receiver
  IntArray*  giveLevel(int num);
 /// Sets the end node of receiver
  void  setEnd(int end) {End=end;}
 /// Returns the end node of receiver
  int     giveEnd() {return End;}
 /// Return root node of receiver
  int     giveRoot() {return Root;}

private:
 /// Computes depth of receiver
  void  computeDepth();
 /// Computes the Width of receiver
  void  computeWidth();

};
#define sloanlevelstruct_h
#endif
