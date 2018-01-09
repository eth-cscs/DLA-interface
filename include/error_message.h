#ifndef DLA_INTERFACE_ERROR_MESSAGE_H
#define DLA_INTERFACE_ERROR_MESSAGE_H

#include <string>
#include <sstream>
#include "util_output.h"

namespace dla_interface {
  namespace util {
    // declaration of internal functions to compose the error message
    template <class... Args>
    std::string errorMessageInternal(const char* file, int line, const char* func,
                                     const Args&... args);
    void errorMessageCont(std::stringstream& s);
    template <class T, class... Args>
    void errorMessageCont(std::stringstream& s, const T& value, const Args&... args);

    // Implementations
    inline void errorMessageCont(std::stringstream& s) {}

    template <class T, class... Args>
    void errorMessageCont(std::stringstream& s, const T& value, const Args&... args) {
      s << value;
      errorMessageCont(s, args...);
    }

    template <class... Args>
    std::string errorMessageInternal(const char* file, int line, const char* func,
                                     const Args&... args) {
      std::stringstream s;
      s << "Calling function " << func << ": Error raised on " << file << ":" << line << ": ";
      errorMessageCont(s, args...);
      return s.str();
    }
  }
}

// Macro to generate the error message which add filename, line and function name.
// The message will have the following structure:
// "Calling function __func__: Error raised on __FILE__:__LINE__: <message>.
// where message is composed concatenating the rest of the arguments with operator<<.
#define errorMessageFunc(func, ...) \
  dla_interface::util::errorMessageInternal(__FILE__, __LINE__, func, __VA_ARGS__)
// Same as errorMessageFunc with func = __func__.
#define errorMessage(...) \
  dla_interface::util::errorMessageInternal(__FILE__, __LINE__, __func__, __VA_ARGS__)

#endif  // DLA_INTERFACE_ERROR_MESSAGE_H
