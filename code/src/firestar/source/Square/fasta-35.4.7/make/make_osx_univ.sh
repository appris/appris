#!/bin/csh

set bin = ../bin
if (! -d ../bin ) mkdir ../bin

cd ../src
if (! -d ../bin/ppc) mkdir ../bin/ppc
rm *.o
make -f ../make/Makefile.os_x all
make -f ../make/Makefile.os_x uinstall

if (! -d ../bin/i386) mkdir ../bin/i386
rm *.o
make -f ../make/Makefile.os_x86 all
make -f ../make/Makefile.os_x86 uinstall

cd ../bin
foreach n ( ppc/* )
set f=$n:t
lipo -create ppc/$f i386/$f -output $f
echo "Universal $f built"
end
rm -rf ppc/ i386/
echo "Done!"

