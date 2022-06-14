#!/bin/sh

cp source/plug_ummdp_abaqus.f source/tmp.f

cat source/ummdp*.f >> source/tmp.f

mv source/tmp.f compiled/ummmdp_vfm.f

cd compiled

rm ummdp_vfm.pyd >nul

python -m numpy.f2py ummdp_vfm.f -m ummdp_vfm -h ummdp_vfm.pyf --overwrite-signature