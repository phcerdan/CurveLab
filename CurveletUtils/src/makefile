
MATLABDIR = /Applications/MATLAB_R2008a/bin
MEX = $(MATLABDIR)/mex
CC = gcc

CFLAGS = -O3 -msse2 -funroll-loops
MEXFLAGS = 

CURVELABDIR = /Users/geback/Code/Curvelets/CurveLab-2.1.1
CURVELABSRC = $(CURVELABDIR)/fdct_wrapping_cpp/src 

FFTWDIR = /Users/geback/Code/fftw-2.1.5

LIBNAME = CurveletUtils
FULLLIBNAME = lib$(LIBNAME).a

LIBS = -lfdct_wrapping -lfftw -l$(LIBNAME)
LIBDIRS = -L$(CURVELABSRC) -L$(FFTWDIR)/fftw/.libs -L.
INCLUDES = -I$(CURVELABSRC) -I$(FFTWDIR)/fftw

#MEXEXT = mexmaci
MEXEXT = $(shell $(MATLABDIR)/mexext)

TARGETS = $(FULLLIBNAME) ../crvlt_extractdirs.$(MEXEXT) ../crvlt_keeplargest.$(MEXEXT) ../crvlt_getmagnitude.$(MEXEXT)  
OBJS = cu_basic_utils.o cu_extractdirs.o cu_thresholding.o

all: $(TARGETS)
	@

../%.$(MEXEXT): %.cpp $(FULLLIBNAME) 
	@$(MEX) $(MEXFLAGS) -outdir .. $< $(INCLUDES) $(LIBDIRS) $(LIBS)

$(FULLLIBNAME): $(OBJS)
	@ar -curs $(FULLLIBNAME) $(OBJS) 
	@ranlib $(FULLLIBNAME)

%.o: %.cpp CurveletUtils.h
	@$(CC) $(CFLAGS) -c $< $(INCLUDES)

clean: 
	rm -f $(TARGETS) $(OBJS)

