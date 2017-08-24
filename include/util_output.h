#ifndef DLA_INTERFACE_UTIL_OUTPUT_H
#define DLA_INTERFACE_UTIL_OUTPUT_H

#include <ostream>
#include "matrix_index.h"

namespace dla_interface {
  namespace util {

    template <class t1, class t2>
    std::ostream& operator<<(std::ostream& s, const std::pair<t1, t2>& value) {
      return s << "(" << value.first << ", " << value.second << ")";
    }

    inline std::ostream& operator<<(std::ostream& s, const Local2DIndex& index) {
      return s << "(" << index.row << ", " << index.col << ")";
    }
    inline std::ostream& operator<<(std::ostream& s, const Global2DIndex& index) {
      return s << "(" << index.row << ", " << index.col << ")";
    }
  }
}

using dla_interface::util::operator<<;

#endif  // DLA_INTERFACE_UTIL_OUTPUT_H
