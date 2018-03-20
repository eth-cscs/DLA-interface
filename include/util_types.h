#ifndef DLA_INTERFACE_UTIL_TYPES_H
#define DLA_INTERFACE_UTIL_TYPES_H

#include <map>
#include "error_message.h"
#include "types.h"

namespace dla_interface {
  namespace util {
    // Returns the number of flop needed
    // to perform 'add' additions and 'mul' multiplications of type T elements.
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
      return NoTrans;
    }

    inline UpLo getUpLo(char uplo) {
      switch (uplo) {
        case 'u':
        case 'U':
          return Upper;
        case 'l':
        case 'L':
          return Lower;
      }
      throw(std::invalid_argument(errorMessage("Wrong UpLo char ", uplo)));
      return Lower;
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
      return RowMajor;
    }

    // Return a string with the name of the solver.
    inline std::string getSolverString(SolverType solver) {
      return solverNames.at(solver);
    }

    namespace internal {
      inline std::map<std::string, SolverType> get_inverse_solver_map() {
        std::map<std::string, SolverType> inv_map;
        for (auto& el : solverNames) {
          inv_map[el.second] = el.first;
        }
        return inv_map;
      }
    }

    // Return the solver which has name str.
    // Throws if str do not correspond to any solver.
    inline SolverType getSolverType(std::string str) {
      static std::map<std::string, SolverType> inv_map = internal::get_inverse_solver_map();
      return inv_map.at(str);
    }

    inline std::ostream& operator<<(std::ostream& out, SolverType solver) {
      return out << getSolverString(solver);
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
  }
}

using dla_interface::util::operator<<;

#endif  // DLA_INTERFACE_UTIL_TYPES_H
