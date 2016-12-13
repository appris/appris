# Install script for directory: /local/jmrodriguez/appris/code/opt/hh-suite/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/local/jmrodriguez/appris/code/opt/hh-suite")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhblits_omp")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhblits"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhmake"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhfilter"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhsearch"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhalign"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhconsensus"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/a3m_extract"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/a3m_database_reduce"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/a3m_database_extract"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/a3m_database_filter"
    "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/cstranslate"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/local/jmrodriguez/appris/code/opt/hh-suite/build/src/cs/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

