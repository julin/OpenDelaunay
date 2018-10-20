macro(find_all_libraries VARNAME LISTNAME PATH SUFFIX)
  set(${VARNAME})
  list(LENGTH ${LISTNAME} NUM_LIST)
  foreach(LIB ${${LISTNAME}})
    if("${PATH}" STREQUAL "")
      find_library(FOUND_LIB ${LIB} PATH_SUFFIXES ${SUFFIX})
    else()
      find_library(FOUND_LIB ${LIB} PATHS ${PATH} NO_DEFAULT_PATH)
    endif()
    if(FOUND_LIB)
      list(APPEND ${VARNAME} ${FOUND_LIB})
    endif()
    unset(FOUND_LIB CACHE)
  endforeach(LIB)
  list(LENGTH ${VARNAME} NUM_FOUND_LIBRARIES)
  if(NUM_FOUND_LIBRARIES LESS NUM_LIST)
    set(${VARNAME})
  endif(NUM_FOUND_LIBRARIES LESS NUM_LIST)
endmacro()



macro(add_boost)
  if(MSVC)
    add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX -D_CRT_SECURE_NO_DEPRECATE 
                    -D_SCL_SECURE_NO_DEPRECATE)
    link_directories("${BOOST_LIB_DIR}")
    foreach(VAR
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
          CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
          CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
      if(${VAR} MATCHES "/RTC1")
          string(REGEX REPLACE "/RTC1" "" ${VAR} "${${VAR}}")
      endif()
    endforeach() 
    string (REPLACE "/DWIN32" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  else()
    add_definitions(-DLINUX)
    
    link_directories("${BOOST_LIB_DIR}")
    set(BOOST_LIBS_REQUIRED
        boost_filesystem boost_thread boost_system boost_timer boost_chrono boost_python)
  
    find_all_libraries(BOOST_LIBS BOOST_LIBS_REQUIRED "${BOOST_LIB_DIR}" "")
  
    message(STATUS "BOOST_LIBS ${BOOST_LIBS}")
    if(UNIX)
      list(APPEND EXTERNAL_LIBRARIES ${BOOST_LIBS} pthread dl )
      if(ENABLE_VALGRIND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -ggdb")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -ggdb")
      endif()
    else()
      list(APPEND EXTERNAL_LIBRARIES ${BOOST_LIBS} )
    endif()
  endif()
endmacro()


macro(set_static_lib)
  if(MSVC)
    add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
    set_target_properties(${MODULE_NAME} PROPERTIES LINK_FLAGS_RELEASE "/nodefaultlib:LIBCMT")
    foreach(VAR
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
          CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
          CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
      if(${VAR} MATCHES "/RTC1")
          string(REGEX REPLACE "/RTC1" "" ${VAR} "${${VAR}}")
      endif()
    endforeach()  
    #message(STATUS "CMAKE_CXX_FLAGS_DEBUG=" ${CMAKE_CXX_FLAGS_DEBUG}) 
  elseif(UNIX)
    if(APPLE)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
	  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    else()
      add_definitions(-DLINUX)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
	  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
      if(ENABLE_VALGRIND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -ggdb")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -ggdb")
      endif()
    endif()
  endif()
endmacro()

macro(configure_runtime)
endmacro()

if(ENABLE_OPENMP)
  include(FindOpenMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-D_USE_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

