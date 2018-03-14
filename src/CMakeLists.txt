
set(DLA_SRC
    communicator_grid.cpp
    communicator_manager.cpp
)

# C/Fortran interface enabled only if ScaLAPACK is available
if(DLA_HAVE_SCALAPACK)
list(APPEND DLA_SRC
     c_dla_interface.cpp
)
endif()
if(DLA_WITH_FORTRAN)
list(APPEND DLA_SRC
     f_dla_interface.F90
)
endif()

add_library(dla_interface ${DLA_SRC})

