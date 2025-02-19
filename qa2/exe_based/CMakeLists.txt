cmake_minimum_required(VERSION 3.10)

# Project name
project(qa2)

SET(EXECUTABLES
    mc_qa2diff
    treeKF_qa2
    mc_qa2
    treeKF_qa2diff
    pid_qa
    graph_qa3
    efficiency
)

SET(SOURCES
    Helper.cpp
    ShapeFitter.cpp
)

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
list(APPEND HEADERS "Shapes.hpp")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)

# Find the ROOT package (you might need to specify the ROOT_DIR if it's not in the default location)
find_package(ROOT REQUIRED)

# Specify where ROOT headers are located
include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR})

set(PCM_FILE_NAME libQa2)

ROOT_GENERATE_DICTIONARY( G__Qa2 ${HEADERS}
    LINKDEF
    Qa2LinkDef.h
    OPTIONS
    -I ${CMAKE_SOURCE_DIR}
)

add_library(Qa2 SHARED ${SOURCES} G__Qa2.cxx)
target_link_libraries(Qa2 PRIVATE ${ROOT_LIBRARIES})

install(TARGETS Qa2 EXPORT Qa2Targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

install(FILES
        ${HEADERS}
        DESTINATION
        include
        COMPONENT
        Devel)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        DESTINATION
        lib
        OPTIONAL)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}.rootmap"
        DESTINATION
        lib
        OPTIONAL)

foreach(EXE ${EXECUTABLES})
    # Add the executable
    add_executable(${EXE} tasks/${EXE}.cpp)
    # Link ROOT libraries to the executable
    target_link_libraries(${EXE} ${ROOT_LIBRARIES} Qa2)
    install (TARGETS ${EXE} RUNTIME DESTINATION bin)
endforeach()
