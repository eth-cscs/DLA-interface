#include <cstddef>
#include <utility>
#include "matrix_index.h"
#include "tile_matrix_tools.h"

namespace tile_matrix_tools {

  std::size_t arrayIndex_0(dla_interface::Local2DIndex ind, std::pair<SizeType, SizeType> block_size,
                           int leading_nr_blocks, std::size_t ld) {
    int full_block_col = ind.col / block_size.second;
    int remaining_block_row = ind.row / block_size.first;
    int block_i = ind.row % block_size.first;
    int block_j = ind.col % block_size.second;

    return (ld * block_size.second) * (full_block_col * leading_nr_blocks + remaining_block_row) +
           block_i + ld * block_j;
  }

  std::size_t arrayIndex(dla_interface::Local2DIndex ind, dla_interface::Local2DIndex base_ind,
                         std::pair<SizeType, SizeType> block_size, int leading_number_blocks, int ld) {
    dla_interface::Local2DIndex full_index(ind.row + base_ind.row, ind.col + base_ind.col);
    return arrayIndex_0(full_index, block_size, leading_number_blocks, ld) -
           arrayIndex_0(base_ind, block_size, leading_number_blocks, ld);
  }

}  // namespace tile_matrix_tools
