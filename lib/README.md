Perl APPRIS version 2.1
=======================

This directory contains the Perl library of APPRIS codebase.

Installation
------------

To install *all* the general Perl APPRIS libraries on your system, 
type the following:
```
   perl Makefile.PL
   make
   make test
   make install
```

Dependencies
------------

There are several dependencies used by APPRIS codebase
that are not cleanly installed during the make install phase 
of installation. 

If you are having any problems automatically getting the dependencies 
to install properly, please either try command line cpan or download
each dependency independently and install them as per their installation
instructions.

Get HTML documentation - Deprecated
----------------------

Get the HTML documentation by means pdoc script:
```
$ softs/pdoc-1.1> perl scripts/perlmod2www.pl \
  -source /home/jmrodriguez/projects/Encode/release_7/lib/Perl/APPRIS/ \
  -target /home/jmrodriguez/projects/Encode/release_7/lib/Perl/docs/ \
  -css data/perl.css
```

