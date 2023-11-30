# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "CMakeFiles\\basic_plus_ImageSeg_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\basic_plus_ImageSeg_autogen.dir\\ParseCache.txt"
  "basic_plus_ImageSeg_autogen"
  )
endif()
