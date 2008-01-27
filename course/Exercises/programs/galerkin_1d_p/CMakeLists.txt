# The name of our project is "galerkin_1d_p".  CMakeLists files in this project can
# refer to the root source directory of the project as ${GALERKIN_1D_P_SOURCE_DIR} and
# to the root binary directory of the project as ${GALERKIN_1D_P_BINARY_DIR}.
project (GALERKIN_1D_P)

add_executable (test_jacobi test_jacobi.cpp)
add_executable (test_legendreset test_legendreset.cpp)
add_executable (test_pset test_pset.cpp)
add_executable (test_fset test_fset.cpp)


