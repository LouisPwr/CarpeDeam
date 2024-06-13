# Install script for directory: /vol/cloud/louis/apps/CarpeDeam15.3/lib/mmseqs

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/vol/cloud/louis/apps/CarpeDeam15.3/build")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/vol/cloud/louis/miniconda3/envs/all/bin/x86_64-conda-linux-gnu-objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/vol/cloud/louis/apps/CarpeDeam15.3/build/lib/mmseqs/lib/microtar/cmake_install.cmake")
  include("/vol/cloud/louis/apps/CarpeDeam15.3/build/lib/mmseqs/lib/cacode/cmake_install.cmake")
  include("/vol/cloud/louis/apps/CarpeDeam15.3/build/lib/mmseqs/lib/alp/cmake_install.cmake")
  include("/vol/cloud/louis/apps/CarpeDeam15.3/build/lib/mmseqs/lib/ksw2/cmake_install.cmake")
  include("/vol/cloud/louis/apps/CarpeDeam15.3/build/lib/mmseqs/data/cmake_install.cmake")
  include("/vol/cloud/louis/apps/CarpeDeam15.3/build/lib/mmseqs/src/cmake_install.cmake")

endif()

