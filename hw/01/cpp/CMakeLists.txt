cmake_minimum_required(VERSION 3.17)
project(geo1004.hw01)

add_executable(hw01 
  src/main.cpp
)
set_target_properties(hw01 PROPERTIES 
  CXX_STANDARD 17
)

# Make the executable a default target to build & run in Visual Studio
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT hw01)

