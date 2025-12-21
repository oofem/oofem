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

#ifndef paramkey_h
#define paramkey_h
 
 
#include <string>
#include <unordered_map>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <iostream>

namespace oofem {
/**
 * ParamKey class is used to create unique keys for parameters in a thread-safe manner.
 * It generates a unique index for each key and allows for debug name registration.
 */
class ParamKey {
public:
    using IndexType = size_t;

    explicit ParamKey(std::string name)
        : name_(std::move(name)), index_(generateUniqueIndex()) {
        registerDebugName(index_, name_);
    }

    IndexType getIndex() const { return index_; }
    const std::string& getName() const { return name_; }
    const char* getNameCStr() const { return name_.c_str(); }

    static std::string getDebugName(IndexType index) {
        std::shared_lock lock(debugMutex_);
        auto it = debugNames_.find(index);
        return it != debugNames_.end() ? it->second : "<unknown>";
    }

    friend std::ostream& operator<<(std::ostream& os, const ParamKey& key) {
        return os << "ParamKey{name=\"" << key.name_ << "\", index=" << key.index_ << "}";
    }

private:
    std::string name_;
    IndexType index_;

    static IndexType generateUniqueIndex() {
        static std::atomic<IndexType> counter{0};
        return counter.fetch_add(1, std::memory_order_relaxed);
    }

    static void registerDebugName(IndexType index, const std::string& name) {
        std::unique_lock lock(debugMutex_);
        debugNames_[index] = name;
    }

    static inline std::unordered_map<IndexType, std::string> debugNames_;
    static inline std::shared_mutex debugMutex_;
};

}   // end namespace oofem

#endif // paramkey_h