
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

add_library(dla_interface ${DLA_SRC})

