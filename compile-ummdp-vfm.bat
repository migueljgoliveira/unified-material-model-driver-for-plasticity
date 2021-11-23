@ECHO OFF
cls

cd source

copy plug_ummdp_vfm.f tmp.f
type ummdp.f >> tmp.f
type ummdp_isotropic.f >> tmp.f
type ummdp_kinematic.f >> tmp.f
type ummdp_print.f >> tmp.f
type ummdp_utl.f >> tmp.f
type ummdp_yfunc.f >> tmp.f
type ummdp_yfunc_bbc2005.f >> tmp.f
type ummdp_yfunc_bbc2008.f >> tmp.f
type ummdp_yfunc_cpb2006.f >> tmp.f
type ummdp_yfunc_gotoh.f >> tmp.f
type ummdp_yfunc_hill1948.f >> tmp.f
type ummdp_yfunc_hill1990.f >> tmp.f
type ummdp_yfunc_hu2005.f >> tmp.f
type ummdp_yfunc_karafillisboyce.f >> tmp.f
type ummdp_yfunc_mises.f >> tmp.f
type ummdp_yfunc_vegter.f >> tmp.f
type ummdp_yfunc_yld2000.f >> tmp.f
type ummdp_yfunc_yld2004.f >> tmp.f
type ummdp_yfunc_yld89.f >> tmp.f"
type ummdp_yfunc_yoshida2011.f >> tmp.f

cd .. 

move "source\tmp.f" "compiled\ummdp_vfm.f"

cd compiled 

del ummdp_vfm.pyd >nul

python -m numpy.f2py -c -m umat ummdp_vfm.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_vfm.log 2>&1
::python -m numpy.f2py -c -m UMAT ummdp_vfm.f

rename umat.*.pyd ummdp_vfm.pyd >nul

cd ..
