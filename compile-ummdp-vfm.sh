#!/bin/sh

cp source/plug_ummdp_abaqus.f source/tmp.f

cat source/ummdp*.f >> source/tmp.f

mv source/tmp.f compiled/ummmdp_vfm.f

cd compiled

rm ummdp_vfm.so >nul

python -m numpy.f2py ummdp_vfm.f -m ummdp_vfm -h ummdp_vfm.pyf --overwrite-signature

python SignatureFile.py

python -m numpy.f2py -c ummdp_vfm2.pyf ummdp_vfm.f --fcompiler=intelem --opt="-heap-arrays -fast -QxHost" 1> ummdp_vfm.log

mv ummdp_vfm.*.so ummdp_vfm.so >nul