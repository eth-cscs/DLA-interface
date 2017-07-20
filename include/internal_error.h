#ifndef DLA_INTERFACE_INTERNAL_ERROR_H
#define DLA_INTERFACE_INTERNAL_ERROR_H

#include <stdexcept>

namespace dla_interface {
  namespace error {
    class InternalError : public std::runtime_error {
      public:
      explicit InternalError(const std::string& what_arg) : std::runtime_error(what_arg) {}
      explicit InternalError(const char* what_arg) : std::runtime_error(what_arg) {}
    };
  }
}

#endif  // DLA_INTERFACE_INTERNAL_ERROR_H
