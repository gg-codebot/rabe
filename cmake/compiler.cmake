message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -ffpe-summary=invalid,zero,overflow")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(STATUS "Setting Intel flags.")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -mkl -heap-arrays -stand f08")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
  message(STATUS "Setting LLVM flang flags.")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp")
  if(APPLE)
    # overwrite the flang incompatible MacOS default flag -dynamiclib
    set(CMAKE_SHARED_MODULE_CREATE_Fortran_FLAGS "-shared")
  endif()
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(MY_DEBUG_FLAG_LIST
  " -O0 -cpp -g"
 " -ffpe-summary=invalid,zero,overflow,underflow "
 " -fcheck=all "
 " -fbacktrace "
 " -finit-real=nan"
 " -Wall "
 " -Wextra "
 " -Wconversion "
 " -fimplicit-none "
 " -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label "
 " -Wno-do-subscript "
 " -Wshadow "
 " -Wno-external-argument-mismatch "  # GCC 15 bug
  )
  string(REPLACE ";" "" MY_DEBUG_FLAG  ${MY_DEBUG_FLAG_LIST})
  set(CMAKE_Fortran_FLAGS_DEBUG "${MY_DEBUG_FLAG}")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "-cpp -mkl -g -warn all -stand f08 ")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
  set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -cpp -g")
endif()
