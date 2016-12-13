# CMake generated Testfile for 
# Source directory: /local/jmrodriguez/appris/code/opt/hh-suite/test
# Build directory: /local/jmrodriguez/appris/code/opt/hh-suite/build/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(hhblits "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhblits" "-h")
SET_TESTS_PROPERTIES(hhblits PROPERTIES  ENVIRONMENT "HHLIB=/local/jmrodriguez/appris/code/opt/hh-suite")
ADD_TEST(hhsearch "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhsearch" "-h")
SET_TESTS_PROPERTIES(hhsearch PROPERTIES  ENVIRONMENT "HHLIB=/local/jmrodriguez/appris/code/opt/hh-suite")
ADD_TEST(hhalign "/local/jmrodriguez/appris/code/opt/hh-suite/build/bin/hhalign" "-h")
SET_TESTS_PROPERTIES(hhalign PROPERTIES  ENVIRONMENT "HHLIB=/local/jmrodriguez/appris/code/opt/hh-suite")
