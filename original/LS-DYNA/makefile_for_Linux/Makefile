#
# Makefile for building lsdyna
# on 'Xeon64 System                     '
# for 'Linux 2.6.32                      '
# using 'Intel Fortran Compiler 13.1 SSE2  '
#
# For double precision may need system BLAS lib.
# Please add following string in the link line.
# 
#
#
FC=ifort
#
FF=-safe_cray_ptr -assume byterecl,buffered_io,protect_parens -warn nousage -zero -ftz -fp-model precise -mP2OPT_hpo_dist_factor=21 -diag-disable 10212,10010 -fimf-arch-consistency=true -nopad -openmp -i8 -r8  -DINTEL -DAdd_ -xSSE2 -align array16byte -fPIC -O2
#
LD=ifort -sox
#
LF= -i-static -L/usr/X11R6/lib64 -lX11 -openmp -lstdc++ -lrt -lstdc++ -lpthread -lstdc++ -lpthread
#
OBJS=init_once.o dynm.o dyn21.o dyn21b.o init_dyn21.o couple2other_user.o dynrfn_user.o adummy_graph.o orderByMetis.o    adummy_msc.o mscapi.o
#
UMMDp_OBJS=plug_ummdp_lsdyna_R931_umat_utan.o ummdp.o ummdp_curve.o ummdp_kinematic.o ummdp_matconst.o ummdp_setprop.o ummdp_utl.o ummdp_yfunc_bbc2005.o ummdp_yfunc_bbc2008.o ummdp_yfunc_cazacu2006.o ummdp_yfunc_control.o ummdp_yfunc_gotoh.o ummdp_yfunc_hill48.o ummdp_yfunc_hill90.o ummdp_yfunc_hu2005.o ummdp_yfunc_karafillisboyce.o ummdp_yfunc_mises.o ummdp_yfunc_vegter.o ummdp_yfunc_yld89.o ummdp_yfunc_yld2000.o ummdp_yfunc_yld2004.o ummdp_yfunc_yoshida2011.o
#
LIBS=libdyna.a libbcsext4.a liblsda.a liblssecurity.a  liblcpack.a libspooles.a libcparse.a libmf2.a liblsm.a liblscrypt.a libresurf.a libsfg.a  libmetis.a libim_rotor_dynamics.a libarpack.a   liblso.a  libem.a libfemster_wrap.a libfemster_wrap2d.a libfemster_wrap1d.a libfemster.a libfemster2d.a libfemster1d.a libpfem.a libcese.a libchemistry.a libparticles.a libeosfl.a libpfem.a libmoving.a libblas.a liblapack.a libfftw3.a  intel64_131_libansysdp.a  libdyna.a libblas.a liblapack.a libmetis.a libspooles.a
#
UMMDp_R931: $(OBJS) $(UMMDp_OBJS)
	$(LD) -o UMMDp_R931 $(OBJS) $(UMMDp_OBJS) $(LIBS) $(LF)
#
init_dyn21.o: init_dyn21.f nhisparm.inc 
	$(FC) -c $(FF) -I. init_dyn21.f
dyn21.o: dyn21.f nhisparm.inc
	$(FC) -c $(FF) -I. dyn21.f
dyn21b.o: dyn21b.f nhisparm.inc
	$(FC) -c $(FF) -I. dyn21b.f
couple2other_user.o: couple2other_user.f 
	$(FC) -c $(FF) -I. couple2other_user.f
dynrfn_user.o: dynrfn_user.f 
	$(FC) -c $(FF) -I. dynrfn_user.f
plug_ummdp_lsdyna_R931_umat_utan.o: plug_ummdp_lsdyna_R931_umat_utan.f 
	$(FC) -c $(FF) -I. plug_ummdp_lsdyna_R931_umat_utan.f
ummdp.o: ummdp.f 
	$(FC) -c $(FF) -I. ummdp.f
ummdp_curve.o: ummdp_curve.f 
	$(FC) -c $(FF) -I. ummdp_curve.f
ummdp_kinematic.o: ummdp_kinematic.f 
	$(FC) -c $(FF) -I. ummdp_kinematic.f
ummdp_matconst.o: ummdp_matconst.f 
	$(FC) -c $(FF) -I. ummdp_matconst.f
ummdp_setprop.o: ummdp_setprop.f 
	$(FC) -c $(FF) -I. ummdp_setprop.f
ummdp_utl.o: ummdp_utl.f 
	$(FC) -c $(FF) -I. ummdp_utl.f
ummdp_yfunc_bbc2005.o: ummdp_yfunc_bbc2005.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_bbc2005.f
ummdp_yfunc_bbc2008.o: ummdp_yfunc_bbc2008.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_bbc2008.f
ummdp_yfunc_cazacu2006.o: ummdp_yfunc_cazacu2006.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_cazacu2006.f
ummdp_yfunc_control.o: ummdp_yfunc_control.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_control.f
ummdp_yfunc_gotoh.o: ummdp_yfunc_gotoh.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_gotoh.f
ummdp_yfunc_hill48.o: ummdp_yfunc_hill48.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_hill48.f
ummdp_yfunc_hill90.o: ummdp_yfunc_hill90.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_hill90.f
ummdp_yfunc_hu2005.o: ummdp_yfunc_hu2005.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_hu2005.f
ummdp_yfunc_karafillisboyce.o: ummdp_yfunc_karafillisboyce.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_karafillisboyce.f
ummdp_yfunc_mises.o: ummdp_yfunc_mises.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_mises.f
ummdp_yfunc_vegter.o: ummdp_yfunc_vegter.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_vegter.f
ummdp_yfunc_yld89.o: ummdp_yfunc_yld89.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_yld89.f
ummdp_yfunc_yld2000.o: ummdp_yfunc_yld2000.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_yld2000.f
ummdp_yfunc_yld2004.o: ummdp_yfunc_yld2004.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_yld2004.f
ummdp_yfunc_yoshida2011.o: ummdp_yfunc_yoshida2011.f 
	$(FC) -c $(FF) -I. ummdp_yfunc_yoshida2011.f
