cmake_minimum_required(VERSION 3.1)
project(variational)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
#set(SRCS main.cpp MeshTransformation.cpp)
#set(HEADERS MeshTransformation.h)
#add_executable(${PROJECT_NAME}_bin ${SRCS} ${HEADERS})

file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})

target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw)
