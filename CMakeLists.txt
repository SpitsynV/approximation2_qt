cmake_minimum_required(VERSION 3.16)
project(approx2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5 REQUIRED COMPONENTS Widgets)

find_package(Threads REQUIRED)

add_executable(approx2d
    main.cpp

    func.cpp       func.h
    task3.cpp      task3.h
    task31.cpp     task31.h
    task5.cpp      task5.h    
    task6.cpp      task6.h

    approximator2d.cpp  approximator2d.h
    plotwidget2d.cpp    plotwidget2d.h
    error.cpp           error.h

    # ── Новые файлы ──────────────────────────────────────────────
    shared_input.h
    method_result.h
    method_runner.cpp   method_runner.h
    method_selector.cpp method_selector.h
    parallel_runner.cpp parallel_runner.h
)

target_link_libraries(approx2d
    Qt5::Widgets
    Threads::Threads
)