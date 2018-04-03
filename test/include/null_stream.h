#ifndef DLA_INTERFACE_TEST_INCLUDE_NULL_STREAM_H
#define DLA_INTERFACE_TEST_INCLUDE_NULL_STREAM_H

namespace testing {
  struct NullStream {
    template <class T>
    NullStream& operator<<(T&&) {
      return *this;
    }

    NullStream& operator<<(std::ostream&(std::ostream&)) {
      return *this;
    }
  };
}

#endif  // DLA_INTERFACE_TEST_INCLUDE_NULL_STREAM_H
