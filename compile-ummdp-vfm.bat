@ECHO OFF
cls

copy "source\plug_ummdp_vfm.f" "source\tmp.f"
type "source\ummdp.f" >> "source\tmp.f"
type "source\ummdp_isotropic.f" >> "source\tmp.f"
type "source\ummdp_kinematic.f" >> "source\tmp.f"
type "source\ummdp_print.f" >> "source\tmp.f"
type "source\ummdp_utl.f" >> "source\tmp.f"
type "source\ummdp_yfunc.f" >> "source\tmp.f"
type "source\ummdp_yfunc_bbc2005.f" >> "source\tmp.f"
type "source\ummdp_yfunc_bbc2008.f" >> "source\tmp.f"
type "source\ummdp_yfunc_cpb2006.f" >> "source\tmp.f"
type "source\ummdp_yfunc_gotoh.f" >> "source\tmp.f"
type "source\ummdp_yfunc_hill1948.f" >> "source\tmp.f"
type "source\ummdp_yfunc_hill1990.f" >> "source\tmp.f"
type "source\ummdp_yfunc_hu2005.f" >> "source\tmp.f"
type "source\ummdp_yfunc_karafillisboyce.f" >> "source\tmp.f"
type "source\ummdp_yfunc_mises.f" >> "source\tmp.f"
type "source\ummdp_yfunc_vegter.f" >> "source\tmp.f"
type "source\ummdp_yfunc_yld2000.f" >> "source\tmp.f"
type "source\ummdp_yfunc_yld2004.f" >> "source\tmp.f"
type "source\ummdp_yfunc_yld89.f" >> "source\tmp.f"
type "source\ummdp_yfunc_yoshida2011.f" >> "source\tmp.f"
move "source\tmp.f" "compiled\UMMDp-vfm.f"

cd compiled 

del UMMDp-vfm.pyd >nul

python -m numpy.f2py -c -m UMAT UMMDp-vfm.f --fcompiler=intelvem --opt=/heap-arrays:0 1> UMMDp-vfm.log 2>&1

rename UMAT.*.pyd UMMDp-vfm.pyd >nul

cd ..
