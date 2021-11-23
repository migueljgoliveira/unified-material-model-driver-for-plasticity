@ECHO OFF
copy "source\plug_ummdp_abaqus.f" "source\tmp.f"
type "source\ummdp*.f" >> "source\tmp.f"
move "source\tmp.f" "compiled\ummdp_abaqus.f"
copy "compiled\ummdp_abaqus.f" "example\ummdp_abaqus.f"