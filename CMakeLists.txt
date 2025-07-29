cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(sFFG4MC LANGUAGES CXX)

find_package(Geant4 CONFIG REQUIRED ui_all vis_all )
find_package(ROOT CONFIG REQUIRED COMPONENTS RIO Tree Core )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Prefer standard C++
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cc")

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cc ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        ${Geant4_LIBRARIES}
        ROOT::RIO
        ROOT::Core
        ROOT::Tree
)

