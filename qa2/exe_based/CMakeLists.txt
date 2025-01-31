cmake_minimum_required(VERSION 3.10)

SET(EXECUTABLES
    mc_qa2diff
    treeKF_qa2
)

# Project name
project(qa2)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)

# Find the ROOT package (you might need to specify the ROOT_DIR if it's not in the default location)
find_package(ROOT REQUIRED)

# Specify where ROOT headers are located
include_directories(${ROOT_INCLUDE_DIRS})

foreach(EXE ${EXECUTABLES})
    # Add the executable
    add_executable(${EXE} ${EXE}.cpp)
    # Link ROOT libraries to the executable
    target_link_libraries(${EXE} ${ROOT_LIBRARIES})
    install (TARGETS ${EXE} RUNTIME DESTINATION bin)
endforeach()
