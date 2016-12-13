# Install script for directory: /local/jmrodriguez/appris/code/opt/hh-suite/scripts

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts" TYPE PROGRAM PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ FILES
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/Align.pm"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/HHPaths.pm"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts" TYPE PROGRAM FILES
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/addss.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/create_profile_from_hhm.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/hhmakemodel.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/hhsuitedb.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/mergeali.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/multithread.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/pdbfilter.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/renumberpdb.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/create_profile_from_hmmer.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/pdb2fasta.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/reformat.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/splitfasta.pl"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/check_a3m.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/ffindex.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/a3m.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/get_a3m_size.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/pdbfilter.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/cif2fasta.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/hhmakemodel.py"
    "/local/jmrodriguez/appris/code/opt/hh-suite/scripts/hh_reader.py"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

