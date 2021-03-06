
if(MSVC)
  if(${MSVC_VERSION} STREQUAL "1800")
    set(ENABLE_INTESIM 1)
    add_definitions(-D_MR_INTESIM)
    set(USE_OCC710 1)
    set(MR_MSVC_VERSION 2013)
  else()
    set(ENABLE_INTESIM 0)
    set(MR_MSVC_VERSION 2017)
  endif()
endif()


if(ENABLE_INTESIM)
  add_definitions(-D_MR_FUNCTION_LEVEL=0)
else()
  add_definitions(-D_MR_FUNCTION_LEVEL=1)
endif()

if(MSVC)
  set(OUTPUT_NAME "output")
  set(DIR_3RD_PARTY "$ENV{DIR_3RDPARTY}")
  if(ENABLE_INTESIM)
    set(BOOST_DIR "${DIR_3RD_PARTY}/intesim/boost/include")
    if(USE_OCC710)
      add_definitions(-D_MR_FUNCTION_LEVEL=1)
      set(OCC_DIR "${DIR_3RD_PARTY}/occ.sdk/7.2/vs2013")
      set(OCC_LIB_DIR_DEBUG "${OCC_DIR}/win64/vc12/libd")
      set(OCC_LIB_DIR_OPT "${OCC_DIR}/win64/vc12/lib")
      add_definitions( -D_MR_OCC_VERSION=710)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2013.win64)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2013.win64)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2013.win64)
    else()
      set(OCC_DIR "${DIR_3RD_PARTY}/intesim/OpenCASCADE")
      set(OCC_LIB_DIR_DEBUG "${OCC_DIR}/libd")
      set(OCC_LIB_DIR_OPT "${OCC_DIR}/lib")
      add_definitions( -D_MR_OCC_VERSION=660)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2013.win64.intesim)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2013.win64.intesim)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2013.win64.intesim)
    endif()
    
    set(QT_DIR "${DIR_3RD_PARTY}/qt4.8.7/vs2013")
    set(OUTPUT_NAME "output2013")
    set(VTK_VER "-7.1")
    set(VTK_DIR "${DIR_3RD_PARTY}/VTK-7.1.1.SDK/vs2013")
    set(VTK_DIR_INC "${VTK_DIR}/include/vtk-7.1")
    set(VTK_LIB_DIR_OPT ${VTK_DIR}/lib)
    set(VTK_LIB_DIR_DEBUG "${VTK_DIR}/libd")
    set(BOOST_LIB_DIR "${DIR_3RD_PARTY}/intesim/boost/lib")
    
  elseif(${MR_MSVC_VERSION} STREQUAL "2017" )
    set(VTK_VER "-7.1")
    set(VTK_DIR "${DIR_3RD_PARTY}/VTK-7.1.1.SDK/vs2017")
    set(BOOST_DIR "${DIR_3RD_PARTY}/boost_1_63_0")
    set(OCC_DIR "${DIR_3RD_PARTY}/occ.sdk/7.2/vs2017")

    set(VTK_DIR_INC "${VTK_DIR}/include/vtk-7.1")
    set(VTK_LIB_DIR_OPT ${VTK_DIR}/lib)
    set(VTK_LIB_DIR_DEBUG ${VTK_DIR}/libd)

    add_definitions( -D_MR_OCC_VERSION=710)

    set(OCC_LIB_DIR_DEBUG "${OCC_DIR}/win64/vc14/libd")
    set(OCC_LIB_DIR_OPT "${OCC_DIR}/win64/vc14/lib")
    set(BOOST_LIB_DIR "${BOOST_DIR}/stage/lib")
    set(ASM_DIR "${DIR_3RD_PARTY}/asm/223.10.2_M140/win64")
    set(ASM_LIB_DIR_DEBUG "${ASM_DIR}/lib/NT_DLLD_A_140-64") 
    set(ASM_LIB_DIR_OPT "${ASM_DIR}/lib/NT_DLL_A_140-64")

    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2017.win64)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2017.win64)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/2017.win64)
  else()
    set(DIR_3RD_PARTY "$ENV{DIR_3RDPARTY}")
    set(CRYPTCPP_LIB_DIR "${DIR_3RD_PARTY}/cryptcpp/win64/2012")
    set(VTK_VER "-7.0")
    set(VTK_DIR "${DIR_3RD_PARTY}/VTK-7.0.0")
    set(BOOST_DIR "${DIR_3RD_PARTY}/boost_1_57_0")
    set(OCC_DIR "${DIR_3RD_PARTY}/opencascade-6.9.0")
    add_definitions( -D_MR_OCC_VERSION=690)

    set(OCC_LIB_DIR_DEBUG "${OCC_DIR}/bin/win64/libd")
    set(OCC_LIB_DIR_OPT "${OCC_DIR}/bin/win64/lib")
    set(BOOST_LIB_DIR "${BOOST_DIR}/stage/win64")
    set(ASM_DIR "${DIR_3RD_PARTY}/asm")
    set(ASM_LIB_DIR_DEBUG "${ASM_DIR}/lib/NT_DLLD110-64") 
    set(ASM_LIB_DIR_OPT "${ASM_DIR}/lib/NT_DLL110-64")

    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/win64)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/win64)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/win64/win64)
  endif()
elseif(MINGW)
  set(DIR_3RD_PARTY "C:/MinGW64")
  add_definitions(-D_MINGW)
  set(VTK_VER "-8.0")
  set(VTK_DIR "${DIR_3RD_PARTY}/include/vtk-8.0")
  set(VTK_OUTPUT_DIR "${VTK_DIR}")

  set(OCC_DIR "${DIR_3RD_PARTY}/include/opencascade")
  add_definitions( -D_MR_OCC_VERSION=710)
  set(OCC_LIB_DIR "${DIR_3RD_PARTY}/lib")
  set(BOOST_DIR "${DIR_3RD_PARTY}/include/boost")
  set(BOOST_LIB_DIR "${DIR_3RD_PARTY}/lib")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/mingw64)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/mingw64)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/mingw64)
elseif(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
  add_definitions(-D_DARWIN)
  set(VTK_VER "-7.1")
  set(BOOST_DIR "${DIR_3RD_PARTY}/boost_1_63_0")
  set(OCC_DIR "${DIR_3RD_PARTY}/occ/7.1.0/output")
  add_definitions( -D_MR_OCC_VERSION=710)

  set(BOOST_LIB_DIR "${BOOST_DIR}/stage/lib")
  set(OCC_LIB_DIR "${OCC_DIR}/mac64/clang/lib")
  set(VTK_DIR "${DIR_3RD_PARTY}/VTK-7.1.0")
  set(VTK_OUTPUT_DIR "${VTK_DIR}/output")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/macx64)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/macx64)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/macx64)

else()
  add_definitions(-D_LINUX)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/linux64)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/linux64)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR_OUTPUT_ROOT}/linux64)
  if(ENABLE_DOCKER_BUILD)
    set(DIR_3RD_PARTY "/app/meshright/3rd")
    set(VTK_VER "-7.1")
    set(VTK_DIR "${DIR_3RD_PARTY}/VTK-7.1.1")
    set(VTK_DIR_INC "${VTK_DIR}/inc")
    set(VTK_OUTPUT_DIR "${VTK_DIR}")
    set(BOOST_DIR "${DIR_3RD_PARTY}/boost_1_63_0")
    set(OCC_DIR "${DIR_3RD_PARTY}/opencascade-7.1.0")
    add_definitions( -D_MR_OCC_VERSION=710)
    set(BOOST_LIB_DIR "${BOOST_DIR}/stage/lib")
    set(OCC_LIB_DIR "${OCC_DIR}/bin/linux64")
  else()
    set(DIR_3RD_PARTY "$ENV{DIR_3RDPARTY}")
    set(VTK_VER "-7.0")
    set(VTK_DIR "/home/julin/VTK-7.0.0")
    set(VTK_OUTPUT_DIR "${VTK_DIR}/output/linux")
    set(BOOST_DIR "${DIR_3RD_PARTY}/boost_1_57_0")
    set(OCC_DIR "${DIR_3RD_PARTY}/opencascade-6.9.0")
    add_definitions( -D_MR_OCC_VERSION=690)
    set(BOOST_LIB_DIR "${BOOST_DIR}/stage/linux64/lib")
    set(OCC_LIB_DIR "${OCC_DIR}/bin/linux64")
  endif()
endif()

set(OPENNURBS_DIR ${DIR_3P}/opennurbs)
set(TETGEN_DIR ${DIR_3P}/tetgen)
set(OPENMESH_DIR ${DIR_3P}/OpenMesh/OpenMesh-6.3)