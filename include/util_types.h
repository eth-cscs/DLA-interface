//
// Distributed Linear Algebra Interface (DLAI)
//
// Copyright (c) 2018-2019, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#ifndef DLA_INTERFACE_UTIL_TYPES_H
#define DLA_INTERFACE_UTIL_TYPES_H

#include <map>
#include "error_message.h"
#include "types.h"

namespace dla_interface {
  namespace util {
    /// Returns the number of flop needed
    /// to perform 'add' additions and 'mul' multiplications of type T elements.
    template <class T>
    constexpr double nrOps(double add, double mul) {
      return TypeInfo<T>::ops_add * add + TypeInfo<T>::ops_mult * mul;
    }

    inline OpTrans getOpTrans(char trans) {
      switch (trans) {
        case 'n':
        case 'N':
          return NoTrans;
        case 't':
        case 'T':
          return Trans;
        case 'c':
        case 'C':
          return ConjTrans;
      }
      throw(std::invalid_argument(errorMessage("Wrong Trans char ", trans)));
    }

    inline UpLo getUpLo(char uplo) {
      switch (uplo) {
        case 'a':
        case 'A':
          return All;
        case 'u':
        case 'U':
          return Upper;
        case 'l':
        case 'L':
          return Lower;
      }
      throw(std::invalid_argument(errorMessage("Wrong UpLo char ", uplo)));
    }

    inline Diag getDiag(char diag) {
      switch (diag) {
        case 'u':
        case 'U':
          return Unit;
        case 'n':
        case 'N':
          return NonUnit;
      }
      throw(std::invalid_argument(errorMessage("Wrong Diag char ", diag)));
    }

    inline Ordering getOrdering(char ordering) {
      switch (ordering) {
        case 'r':
        case 'R':
          return RowMajor;
        case 'c':
        case 'C':
          return ColMajor;
      }
      throw(std::invalid_argument(errorMessage("Wrong Ordering char ", ordering)));
    }

    /// Return a string with the name of the solver.
    inline std::string getSolverString(SolverType solver) {
      return solver_names.at(solver);
    }

    /// Return a string with the name of the distribution.
    inline std::string getDistributionString(DistributionType dist) {
      return dist_names.at(dist);
    }

    namespace internal {
      template <class T, class U>
      inline std::map<T, U> get_inverse_map(const std::map<U, T>& map) {
        std::map<T, U> inv_map;
        for (auto& el : map) {
          inv_map[el.second] = el.first;
        }
        return inv_map;
      }
    }

    /// Return the solver which has name str.
    ///
    /// @throws if str do not correspond to any solver.
    inline SolverType getSolverType(std::string str) {
      static std::map<std::string, SolverType> inv_map = internal::get_inverse_map(solver_names);
      return inv_map.at(str);
    }

    /// Return the distribution which has name str.
    ///
    /// @throws if str do not correspond to any distribution.
    inline DistributionType getDistributionType(std::string str) {
      static std::map<std::string, DistributionType> inv_map = internal::get_inverse_map(dist_names);
      return inv_map.at(str);
    }

    inline std::ostream& operator<<(std::ostream& out, SolverType solver) {
      return out << getSolverString(solver);
    }

    inline std::ostream& operator<<(std::ostream& out, DistributionType dist) {
      return out << getDistributionString(dist);
    }

    inline std::ostream& operator<<(std::ostream& out, OpTrans trans) {
      return out << (char)trans;
    }

    inline std::ostream& operator<<(std::ostream& out, UpLo uplo) {
      return out << (char)uplo;
    }

    inline std::ostream& operator<<(std::ostream& out, Ordering ordering) {
      return out << (char)ordering;
    }

    template <class T>
    void debug_print(std::ostream& out, T obj) {
      out << obj;
    }

    template <class ElType>
    OpTrans RemoveConjOpTypeIfReal(OpTrans trans) {
      if (std::is_same<ElType, BaseType<ElType>>::value && trans == ConjTrans)
        return Trans;
      return trans;
    }
  }
}

using dla_interface::util::operator<<;

#endif  // DLA_INTERFACE_UTIL_TYPES_H
