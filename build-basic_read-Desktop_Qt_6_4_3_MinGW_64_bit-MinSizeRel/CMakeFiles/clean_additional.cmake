# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "MinSizeRel")
  file(REMOVE_RECURSE
  "CMakeFiles\\basic_read_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\basic_read_autogen.dir\\ParseCache.txt"
  "basic_read_autogen"
  )
endif()
