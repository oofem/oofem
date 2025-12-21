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

#include "inputrecord.h"
#include "error.h"
#include "datareader.h"
#include <sstream>
#include <iomanip>
#include <set>


namespace oofem {

// maximum message length for error messages (otherwise ellipsized)
constexpr static int maxMsgLen=80;

/* MIT-licensed, from https://github.com/guilhermeagostinelli/levenshtein/blob/master/levenshtein.cpp */
int InputRecord::giveLevenshteinDist(const std::string& word1, const std::string& word2){
    int size1 = word1.size();
    int size2 = word2.size();
    // Verification matrix i.e. 2D array which will store the calculated distance.
    // int verif[size1 + 1][size2 + 1];
    std::vector<std::vector<int>> verif(size1+1);
    for(int i=0; i<=size1; i++){ verif[i].resize(size2+1); }


    // If one of the words has zero length, the distance is equal to the size of the other word.
    if (size1 == 0)
        return size2;
    if (size2 == 0)
        return size1;

    // Sets the first row and the first column of the verification matrix with the numerical order from 0 to the length of each word.
    for (int i = 0; i <= size1; i++)
        verif[i][0] = i;
    for (int j = 0; j <= size2; j++)
        verif[0][j] = j;

    // Verification step / matrix filling.
    for (int i = 1; i <= size1; i++) {
        for (int j = 1; j <= size2; j++) {
            // Sets the modification cost.
            // 0 means no modification (i.e. equal letters) and 1 means that a modification is needed (i.e. unequal letters).
            int cost = (word2[j - 1] == word1[i - 1]) ? 0 : 1;

            // Sets the current position of the matrix as the minimum value between a (deletion), b (insertion) and c (substitution).
            // a = the upper adjacent value plus 1: verif[i - 1][j] + 1
            // b = the left adjacent value plus 1: verif[i][j - 1] + 1
            // c = the upper left adjacent value plus the modification cost: verif[i - 1][j - 1] + cost
            verif[i][j] = std::min(
                std::min(verif[i - 1][j] + 1, verif[i][j - 1] + 1),
                              verif[i - 1][j - 1] + cost
            );
        }
    }
    // The last position of the matrix will contain the Levenshtein distance.
    return verif[size1][size2];
}

std::string InputRecord::error_msg_with_hints(const std::string& val, const std::map<int,std::vector<std::string>>& v2nn) {
     //std::string InputRecord::error_msg_with_hints(const std::string& val, const std::vector<std::string>& all_names){
    std::ostringstream oss, oss2;
    std::string minName; int minDist=1000;
    for(const auto& [num,names]: v2nn){
        oss2<<"  "<<std::setfill(' ')<<std::setw(3)<<num<<":";
        for(const auto& name: names){
            oss2<<" "<<name;
            if(!val.empty()){
                if(int ld=giveLevenshteinDist(val,name); ld<minDist){ minDist=ld; minName=name; }
            }
        }
        oss2<<"\n";
    }
    if(minDist<=4) oss<<" (did you mean '"<<minName<<"'?)";
    oss<<".\nPossible values:\n"<<oss2.str();
    return oss.str();
}


#ifdef _USE_TRACE_FIELDS
    bool InputRecord::TraceFields::active=false;
    std::ofstream InputRecord::TraceFields::out;
    void InputRecord::TraceFields::write(const std::string& s){
        if(!InputRecord::TraceFields::active) return;
        size_t hash=std::hash<std::string>{}(s);
        static std::set<size_t> written;
        if(written.count(hash)>0) return;
        written.insert(hash);
        InputRecord::TraceFields::out<<s<<std::endl;
    }

    void InputRecord::traceEnum(const std::string& name, const std::map<int,std::vector<std::string>>& val2names){
        if(!InputRecord::TraceFields::active) return;
        std::ostringstream rec;
        // write as tag;id;json
        rec<<"~Enum~;"<<name<<";";
        int i=0;
        for(const auto& kvv: val2names){
            rec<<(i++==0?'{':',')<<'"'<<kvv.first<<"\":";
            int j=0;
            for(const auto& v: kvv.second) rec<<(j++==0?'[':',')<<'"'<<v<<'"';
            rec<<']';
        }
        rec<<'}';
        InputRecord::TraceFields::write(rec.str());
    }
#endif

InputRecord :: InputRecord(DataReader* r){
    reader = r;
}

DataReader*
InputRecord :: giveReader() const {
    return reader;
}

std::shared_ptr<InputRecord>
InputRecord::ptr() {
    // we could just return this->shared_from_this, but it throws std::bad_weak_ptr (with no backtrace)
    // so we do essentially the same, but provide a nice message and abort immediately
    auto weak=this->weak_from_this();
    std::shared_ptr<InputRecord> ret=weak.lock();
    if(!ret) OOFEM_ERROR("shared_ptr<InputRecord>::ptr(): object lifetime expired (programming error)");
    return ret;
}


#ifdef _INPUTRECORD_OPTIONAL_OLD
void
InputRecord :: giveOptionalField(int &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(double &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(bool &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(std :: string &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(FloatArray &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(IntArray &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(FloatMatrix &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(std :: vector< std :: string > &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(Dictionary &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(std :: list< Range > &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}

void
InputRecord :: giveOptionalField(ScalarFunction &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch ( MissingKeywordInputException & ) { }
    }
}
#endif


InputException::InputException(const InputRecord& ir, std::string keyword, int number) : 
    record(ir.giveRecordAsString()), keyword(std::move(keyword)), number(number)
{ }


MissingKeywordInputException::MissingKeywordInputException(const InputRecord& ir, std::string kw, int n) :
    InputException(ir, std::move(kw), n)
{
    msg = ir.giveLocation()+": missing keyword \"" + keyword + "\"" \
          "\n  \"" + record.substr(0, maxMsgLen) + (record.size()>maxMsgLen?"...":"")+ "\"";
}


BadFormatInputException::BadFormatInputException(const InputRecord& ir, std::string kw, int n) :
    InputException(ir, std::move(kw), n)
{
    msg = ir.giveLocation()+": bad format for keyword \"" + keyword + "\"" \
          "\n   \"" + record.substr(0, maxMsgLen) + (record.size()>maxMsgLen?"...":"") + "\"";
}


ValueInputException::ValueInputException(const InputRecord& ir, std::string kw, const std::string &reason) :
    InputException(ir, std::move(kw), -1)
{
    msg = ir.giveLocation()+": value input error for keyword \"" + keyword + "\"" + \
          "\nReason: \"" + reason + "\"\n" + \
          "\n  \"" + record.substr(0, maxMsgLen) + (record.size()>maxMsgLen?"...":"") + "\"";
}


const char* MissingKeywordInputException::what() const noexcept
{ 
    return msg.c_str();
}


const char* BadFormatInputException::what() const noexcept
{
    return msg.c_str();
}

const char* ValueInputException::what() const noexcept
{
    return msg.c_str();
}


ComponentInputException::ComponentInputException(const std::string keyword, ComponentType ct, int number, const std::string &reason)
{
    std::string cname;
    if (ct==ctDofManager) {
        cname = "DofManager";
    } else if (ct ==ctElement) {
        cname = "Element";
    }
    msg = "Value input error for keyword \"" + keyword + "\"" + \
          "\nReason: \"" + reason + "\"\n" + \
          cname + "[" + std::to_string(number) + "]" + "\"";
}


ComponentInputException::ComponentInputException(ComponentType ct, int number, const std::string &reason)
{
    std::string cname;
    if (ct==ctDofManager) {
        cname = "DofManager";
    } else if (ct ==ctElement) {
        cname = "Element";
    }
    msg = "Value input error, reason: \"" + reason + "\"\n" + \
          cname + "[" + std::to_string(number) + "]" + "\"";
}

const char* ComponentInputException::what() const noexcept
{ 
    return msg.c_str();
}

} // end namespace oofem
