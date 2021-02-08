#!/bin/sh
cp source/plug_ummdp_abaqus.f source/tmp.f
cat source/ummdp*.f >> source/tmp.f
mv source/tmp.f compiled/UMMDp.f
cp compiled/UMMDp.f example/UMMDp.f