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
type ummdp_yield_yld89.f >> tmp.f
type ummdp_yield_yoshida2011.f >> tmp.f

cd ..

move "source\tmp.f" "build\ummdp_vfm.f"

cd build 

del ummdp_vfm.pyd >nul

:: ifort optimization flags
:: /heap-arrays temporary arrays of minimum size n (in kilobytes) are allocated in heap memory rather than on the stack
:: /Ofast = /O3 /Qprec-div- /fp:fast=2
:: /QxHost : generate instructions for the highest instruction set and processor available on the compilation host machine
:: python -m numpy.f2py -c -m ummdp_vfm ummdp_vfm.f --fcompiler=intelvem --opt="/heap-arrays /Ofast /QxHost" 1> ummdp_vfm.log
python -m numpy.f2py ummdp_vfm.f -m ummdp_vfm -h ummdp_vfm.pyf --overwrite-signature

call python SignatureFile.py ummdp_vfm

python -m numpy.f2py -c ummdp_vfm2.pyf ummdp_vfm.f --fcompiler=intelvem --opt="/heap-arrays /Ofast /QxHost" 1> ummdp_vfm.log

rename ummdp_vfm.*.pyd ummdp_vfm.pyd >nul

cd ..

copy "build\ummdp_vfm.pyd" "release\ummdp_vfm.pyd"
copy "build\ummdp_vfm.f" "release\ummdp_vfm.f"