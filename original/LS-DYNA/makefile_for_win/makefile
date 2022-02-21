#
# Makefile for building lsdyna.exe
# on 'Windows X64
# SMPD
FC=ifort
FLAGS = -c -traceback  -4Yportlib -fp:precise -assume:byterecl,buffered_io -QxW -w -WB \
-fpp -DPCWIN -DINTEL -DIA64 -DAMD64 -DWINX64 -DNET_SECURITY -DXEON64 -DLSTCODE -DLONGFORMAT 

TARGET=UMMDp_R931.exe  
SMPS = -DOPENMP -Qopenmp -DDOUBLEXR 
SMPD = -DOPENMP -Qopenmp -DAUTODOUBLE -4R8 -4I8 

FFLAGS = $(FLAGS) $(SMPD)
LFLAGS = -F:1610612736 -link -force -nodefaultlib:msvcrt.lib -nodefaultlib:vcomp.lib 
OBJS= init_dyn21.obj dyn21.obj dyn21b.obj couple2other_user.obj
OBJS_UMMDp=  plug_ummdp_lsdyna_R931_umat_utan.obj\
             ummdp.obj ummdp_utl.obj\
             ummdp_curve.obj\
             ummdp_kinematic.obj\
             ummdp_setprop.obj\
             ummdp_yfunc_control.obj\
             ummdp_yfunc_bbc2005.obj\
             ummdp_yfunc_bbc2008.obj\
             ummdp_yfunc_cazacu2006.obj\
             ummdp_yfunc_gotoh.obj\
             ummdp_yfunc_hill48.obj\
             ummdp_yfunc_hill90.obj\
             ummdp_yfunc_karafillisboyce.obj\
             ummdp_yfunc_mises.obj\
             ummdp_yfunc_vegter.obj\
             ummdp_yfunc_yld2000.obj\
             ummdp_yfunc_yld2004.obj\
             ummdp_yfunc_yld89.obj\
             ummdp_yfunc_yoshida2011.obj\
             ummdp_yfunc_hu2005.obj
LIBS= libdyna.lib libansys.lib shell32.lib user32.lib comctl32.lib comdlg32.lib Iphlpapi.lib \
gdi32.lib mkl_intel_ilp64.lib mkl_sequential.lib mkl_core.lib

$(TARGET): $(OBJS) $(OBJS_UMMDp)
	$(FC) -w -q -o $(TARGET) $(OBJS) $(OBJS_UMMDp) $(LIBS) $(LFLAGS)
$(OBJS) :nhisparm.inc
.F.obj:  
 	$(FC) $(FFLAGS) $*.F
$(OBJS_UMMDp) :nhisparm.inc
.F.obj:  
 	$(FC) $(FFLAGS) $*.F

clean:
	-if exist *.obj erase *.obj
#	-if exist del $(OBJS)
	-if exist $(TARGET) erase $(TARGET)

