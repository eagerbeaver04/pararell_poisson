cmake_minimum_required(VERSION 3.10)

project(Project C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED TRUE)

add_executable(main_executable src/main.c)

target_compile_options(main_executable PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wformat=2
    -Wconversion
    -O2
    -g
)
