cmake_minimum_required(VERSION 3.25)
project(Day2)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} src/main.cpp src/Image.cpp src/Ray.h src/Sphere.cpp src/Sphere.h src/Film.h src/Camera.h src/Camera.cpp src/Body.h src/Renderer.cpp src/Renderer.h src/Material.h src/Color.cpp
        src/Box.h
        src/Box.cpp
        src/Triangle.h
        src/Triangle.cpp)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

find_package(OpenCV CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC opencv_ml opencv_dnn opencv_core opencv_flann opencv_imgcodecs opencv_highgui)

#for OpenMP
find_package(OpenMP)
if(OPENMP_CXX_FOUND)
    # add flags for OpenMP
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
