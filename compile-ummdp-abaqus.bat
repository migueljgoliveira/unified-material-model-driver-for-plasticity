@ECHO OFF
cls

copy "source\plug_ummdp_abaqus.f" "source\tmp.f"
type "source\ummdp*.f" >> "source\tmp.f"
move "source\tmp.f" "build\ummdp_abaqus.f"

:: compile fortran umat
cd build

call abaqus make library=ummdp_abaqus

@REM move ummdp_abaqus-std.obj ummdp_abaqus.obj

cd ..

copy "build\ummdp_abaqus.f" "release\ummdp_abaqus.f"
copy "build\ummdp_abaqus-std.obj" "release\ummdp_abaqus.obj"
copy "build\ummdp_abaqus.f" "example\ummdp_abaqus.f"