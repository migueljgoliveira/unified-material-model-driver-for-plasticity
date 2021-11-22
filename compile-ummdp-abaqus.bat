@ECHO OFF
copy "source\plug_ummdp_abaqus.f" "source\tmp.f"
type "source\ummdp*.f" >> "source\tmp.f"
move "source\tmp.f" "compiled\UMMDp-abaqus.f"
copy "compiled\UMMDp-abaqus.f" "example\UMMDp-abaqus.f"