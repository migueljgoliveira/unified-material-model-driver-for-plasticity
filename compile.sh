#!/bin/sh
cp _sourcecode/plug_ummdp_abaqus.f _sourcecode/tmp.f
cat _sourcecode/ummdp*.f >> _sourcecode/tmp.f
mv _sourcecode/tmp.f _compiled/UMMDp.f