cmake_minimum_required(VERSION 3.1)

project(2DVectorFieldDecom)

find_package(OpenGL)

find_package(Qt5Widgets)
find_package(Qt5OpenGL)
set(CMAKE_AUTOMOC ON)

find_package(CGAL)

include(${CGAL_USE_FILE})

add_subdirectory(glm)

include_directories(".")
include_directories("eigen")

file(GLOB Headers "src/*.h")
file(GLOB Sources "src/*.cpp")
file(GLOB Resources "resources/*.png")

add_executable( 2DVectorFieldDecom ${Headers} ${Sources} ${Resources})

qt5_use_modules( 2DVectorFieldDecom Widgets)
qt5_use_modules( 2DVectorFieldDecom OpenGL)

target_link_libraries( 2DVectorFieldDecom ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
target_link_libraries( 2DVectorFieldDecom ${QT_LIBRARIES} ${OPENGL_LIBRARIES})