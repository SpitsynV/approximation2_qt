cmake_minimum_required(VERSION 3.16)
project(approx2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5 REQUIRED COMPONENTS Widgets)

add_executable(approx2d
    main.cpp
    func.cpp      func.h
    task3.cpp     task3.h
    task31.cpp    task31.h
    approximator2d.cpp  approximator2d.h
    plotwidget2d.cpp    plotwidget2d.h
    error.cpp   error.h
    task5.cpp
)

target_link_libraries(approx2d Qt5::Widgets)