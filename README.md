# Compile-time pathtracer for C++20

## Description
This project was inspired by [Jason Turners Homework Assignments](https://www.youtube.com/watch?v=cpdjQiRxEJ8).
The goal is to create a pathtracer that generates an image at compile-time and writes it to a file at runtime.

## Current problems
GCC and Clang are not capable of evaluating **HUGE** amounts of constexpr code.

## TODO
- Implement constexpr string
