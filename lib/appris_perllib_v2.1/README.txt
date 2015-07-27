Perl APPRIS version 1.00

======================

This directory contains the Perl library of APPRIS codebase.

INSTALLATION

To install *all* the general Perl APPRIS libraries on your system, 
type the following:

   perl Makefile.PL
   #perl Makefile.PL INSTALL_BASE=~/tmp
   make
   make test
   make install

DEPENDENCIES

There are several dependencies used by APPRIS codebase
that are not cleanly installed during the make install phase 
of installation. 

If you are having any problems automatically getting the dependencies 
to install properly, please either try command line cpan or download
each dependency independently and install them as per their installation
instructions.


COPYRIGHT AND LICENCE

Copyright (C) 2003 BioMoby Developers Group (www.biomoby.org)

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 



GET HTML DOCUMENTATION

Get the HTML documentation by means pdoc script:
	softs/pdoc-1.1> perl scripts/perlmod2www.pl -source /home/jmrodriguez/projects/Encode/release_7/lib/Perl/APPRIS/ -target /home/jmrodriguez/projects/Encode/release_7/lib/Perl/docs/ -css data/perl.css

