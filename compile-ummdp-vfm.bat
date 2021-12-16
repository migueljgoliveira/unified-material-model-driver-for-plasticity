@ECHO OFF
cls

cd source

copy plug_ummdp_vfm.f tmp.f
type ummdp.f >> tmp.f
type ummdp_isotropic.f >> tmp.f
type ummdp_kinematic.f >> tmp.f
type ummdp_print.f >> tmp.f
type ummdp_utility.f >> tmp.f
type ummdp_yield.f >> tmp.f
type ummdp_yield_bbc2005.f >> tmp.f
type ummdp_yield_bbc2008.f >> tmp.f
type ummdp_yield_cpb2006.f >> tmp.f
type ummdp_yield_gotoh.f >> tmp.f
type ummdp_yield_hill1948.f >> tmp.f
type ummdp_yield_hill1990.f >> tmp.f
type ummdp_yield_hu2005.f >> tmp.f
type ummdp_yield_karafillisboyce.f >> tmp.f
type ummdp_yield_mises.f >> tmp.f
type ummdp_yield_vegter.f >> tmp.f
type ummdp_yield_yld2000.f >> tmp.f
type ummdp_yield_yld2004.f >> tmp.f
type ummdp_yield_yld89.f >> tmp.f"
type ummdp_yield_yoshida2011.f >> tmp.f

cd .. 

move "source\tmp.f" "compiled\ummdp_vfm.f"

cd compiled 

del ummdp_vfm.pyd >nul

python -m numpy.f2py -c -m ummdp_vfm ummdp_vfm.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_vfm.log
rem 

rename ummdp_vfm.*.pyd ummdp_vfm.pyd >nul

cd ..
