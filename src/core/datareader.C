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

#include "datareader.h"


namespace oofem{

InputRecord *DataReader::giveChildRecord( const std::shared_ptr<InputRecord> &ir, InputFieldType ift, const std::string &name, InputRecordType irType, bool optional )
{
    if ( ir->hasChild( ift, name, optional ) ) return &( this->giveInputRecord( irType, /*recordId*/ 1 ) );
    return nullptr;
};

DataReader::GroupRecords DataReader::giveGroupRecords( const std::shared_ptr<InputRecord> &ir, InputFieldType ift, const std::string &name, InputRecordType irType, bool optional )
{
    return GroupRecords( *this, name, irType, ir->giveGroupCount( ift, name, optional ) );
}

DataReader::GroupRecords DataReader::giveGroupRecords(const std::string& name, InputRecordType irType, int numRequired)
{
    int count=giveGroupCount(name);
    if(count>=0 && numRequired>=0 && count!=numRequired) OOFEM_ERROR("Mismatch in %s: %d records of type '%s' required, %d found.",giveReferenceName().c_str(),numRequired,name.c_str(),count);
    return GroupRecords(*this,name,irType,numRequired>=0?numRequired:count);
}

DataReader::GroupRecords::Iterator::Iterator( DataReader &dr_, const std::string &group_, InputRecordType irType_, int size_, int index_ ) :
    dr( dr_ ), group( group_ ), irType( irType_ ), size( size_ ), index( index_ )
{
    /* read the first line for begin(), unless there are no lines to read at all */
    if ( index == 0 && size_ > 0 ) {
        if ( !this->group.empty() ) {
            entered = true;
            dr.enterGroup( this->group );
        }
        irPtr = &dr.giveInputRecord( irType, /*recordId*/ 1 );
        if ( irPtr ) dr.enterRecord( irPtr );
    }
    #if 0
        // open 0-sized group to mark it as processed
        if(index==0 && size_==0 && !this->group.empty() && dr.hasGroup(this->group)){ dr.enterGroup(this->group); dr.leaveGroup(this->group); }
    #endif
};


DataReader::GroupRecords::Iterator &DataReader::GroupRecords::Iterator::operator++()
{
    if ( irPtr ) dr.leaveRecord( irPtr );
    index++;
    /* don't read past the last line, assign nullptr instead */
    if ( index >= size ) {
        irPtr = nullptr;
        if ( entered ) dr.leaveGroup( this->group );
    } else {
        irPtr = &dr.giveInputRecord( irType, /*recordId*/ index + 1 );
        if ( irPtr ) dr.enterRecord( irPtr );
    }
    return *this;
}


DataReader::GroupRecords::GroupRecords( DataReader &dr_, const std::string &group_, InputRecordType irType_, int size ) :
    dr( dr_ ), group( group_ ), irType( irType_ ), size_( size ){};



}; // namespace oofem
