#ifndef DLA_INTERFACE_TEST_INCLUDE_TILE_MATRIX_TOOLS_H
#define DLA_INTERFACE_TEST_INCLUDE_TILE_MATRIX_TOOLS_H

#include <cstddef>
#include <utility>
#include "matrix_index.h"
#include "types.h"

namespace tile_matrix_tools {
  using dla_interface::SizeType;

  std::size_t arrayIndex(dla_interface::Local2DIndex ind, dla_interface::Local2DIndex base_ind,
                         std::pair<SizeType, SizeType> block_size, int leading_number_blocks, int ld);

}  // namespace tile_matrix_tools

#endif  // DLA_INTERFACE_TEST_INCLUDE_TILE_MATRIX_TOOLS_H
