@ECHO OFF
copy "source\plug_ummdp_abaqus.f" "source\tmp.f"
type "source\ummdp*.f" >> "source\tmp.f"
move "source\tmp.f" "compiled\ummdp_abaqus.f"
copy "compiled\ummdp_abaqus.f" "example\ummdp_abaqus.f"
:: compile fortran umat
cd compiled
call abaqus make library=ummdp_abaqus
move ummdp_abaqus-std.obj ummdp_abaqus.obj
del standardU.dll
cd ..