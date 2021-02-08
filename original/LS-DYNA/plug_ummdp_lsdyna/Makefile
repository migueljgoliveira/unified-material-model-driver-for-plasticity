#
# Makefile for building lsdyna
# on 'Xeon64 System                     '
# for 'Linux 2.6.18                      '
# using 'Intel Fortran Compiler 13.1 SSE2  '
#
lsdyna: init_dyn21.o plug_ummdp_lsdyna_other.o dyn21b.o couple2other_user.o dynrfn_user.o plug_ummdp_lsdyna_umat_utan.o ummdp.o ummdp_curve.o ummdp_kinematic.o ummdp_matconst.o ummdp_setprop.o ummdp_utl.o ummdp_yfunc_bbc2005.o ummdp_yfunc_bbc2008.o ummdp_yfunc_cazacu2006.o ummdp_yfunc_control.o ummdp_yfunc_gotoh.o ummdp_yfunc_hill48.o ummdp_yfunc_hill90.o ummdp_yfunc_hu2005.o ummdp_yfunc_karafillisboyce.o ummdp_yfunc_mises.o ummdp_yfunc_vegter.o ummdp_yfunc_yld89.o ummdp_yfunc_yld2000.o ummdp_yfunc_yld2004.o ummdp_yfunc_yoshida2011.o
	ifort -o lsdyna init_once.o init_keywrd.o init_dyn21.o dynm.o plug_ummdp_lsdyna_other.o dyn21b.o couple2other_user.o dynrfn_user.o plug_ummdp_lsdyna_umat_utan.o ummdp.o ummdp_curve.o ummdp_kinematic.o ummdp_matconst.o ummdp_setprop.o ummdp_utl.o ummdp_yfunc_bbc2005.o ummdp_yfunc_bbc2008.o ummdp_yfunc_cazacu2006.o ummdp_yfunc_control.o ummdp_yfunc_gotoh.o ummdp_yfunc_hill48.o ummdp_yfunc_hill90.o ummdp_yfunc_hu2005.o ummdp_yfunc_karafillisboyce.o ummdp_yfunc_mises.o ummdp_yfunc_vegter.o ummdp_yfunc_yld89.o ummdp_yfunc_yld2000.o ummdp_yfunc_yld2004.o ummdp_yfunc_yoshida2011.o adummy_graph.o orderByMetis.o adummy_msc.o mscapi.o libdyna.a libbcsext4.a liblsda.a liblssecurity.a  liblcpack.a libspooles.a libcparse.a libmf2.a liblsm.a liblscrypt.a libresurf.a libsfg.a libmetis.a  libfemster_wrap.a libfemster_wrap2d.a libfemster_wrap1d.a libfemster.a libfemster2d.a libfemster1d.a libpfem.a libmetis.a libblas.a liblapack.a libfftw3.a  intel64_131_libansysdp.a  libdyna.a libblas.a liblapack.a libmetis.a  -i-static -L/usr/X11R6/lib64 -lX11 -openmp  -lrt -lstdc++ 
#
init_dyn21.o: init_dyn21.f nhisparm.inc 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. init_dyn21.f
plug_ummdp_lsdyna_other.o: plug_ummdp_lsdyna_other.f nhisparm.inc 
	ifort -c -fpe0 -traceback -g -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O0 -I. plug_ummdp_lsdyna_other.f
dyn21b.o: dyn21b.f nhisparm.inc 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. dyn21b.f
couple2other_user.o: couple2other_user.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. couple2other_user.f
dynrfn_user.o: dynrfn_user.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. dynrfn_user.f
plug_ummdp_lsdyna_umat_utan.o: plug_ummdp_lsdyna_umat_utan.f 
	ifort -c -fpe0 -traceback -g -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O0 -I. plug_ummdp_lsdyna_umat_utan.f
ummdp.o: ummdp.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp.f
ummdp_curve.o: ummdp_curve.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_curve.f
ummdp_kinematic.o: ummdp_kinematic.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_kinematic.f
ummdp_matconst.o: ummdp_matconst.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_matconst.f
ummdp_setprop.o: ummdp_setprop.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_setprop.f
ummdp_utl.o: ummdp_utl.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_utl.f
ummdp_yfunc_bbc2005.o: ummdp_yfunc_bbc2005.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_bbc2005.f
ummdp_yfunc_bbc2008.o: ummdp_yfunc_bbc2008.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_bbc2008.f
ummdp_yfunc_cazacu2006.o: ummdp_yfunc_cazacu2006.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_cazacu2006.f
ummdp_yfunc_control.o: ummdp_yfunc_control.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_control.f
ummdp_yfunc_gotoh.o: ummdp_yfunc_gotoh.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_gotoh.f
ummdp_yfunc_hill48.o: ummdp_yfunc_hill48.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_hill48.f
ummdp_yfunc_hill90.o: ummdp_yfunc_hill90.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_hill90.f
ummdp_yfunc_hu2005.o: ummdp_yfunc_hu2005.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_hu2005.f
ummdp_yfunc_karafillisboyce.o: ummdp_yfunc_karafillisboyce.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_karafillisboyce.f
ummdp_yfunc_mises.o: ummdp_yfunc_mises.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_mises.f
ummdp_yfunc_vegter.o: ummdp_yfunc_vegter.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_vegter.f
ummdp_yfunc_yld89.o: ummdp_yfunc_yld89.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_yld89.f
ummdp_yfunc_yld2000.o: ummdp_yfunc_yld2000.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_yld2000.f
ummdp_yfunc_yld2004.o: ummdp_yfunc_yld2004.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_yld2004.f
ummdp_yfunc_yoshida2011.o: ummdp_yfunc_yoshida2011.f 
	ifort -c -safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -axSSE2 -align array16byte -nopad -fpp2 -openmp -i8 -r8  -DINTEL -DAdd_ -O2 -I. ummdp_yfunc_yoshida2011.f
