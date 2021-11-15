OPTIMIZATION = -O3
FFLAGS = $(HEADERS:%=-I%) -fbounds-check -mcmodel=medium -fopenmp -fPIC
LIBS = -lm -lc -lstdc++

FHEAD = BSE/const_bse.h BSE/zdata.h

FC = gfortran $(OPTIMIZATION)
CC = gcc $(OPTIMIZATION)

$(BSE)/%.o: $(BSE)/%.f -c $<


FSOURCE = \
	BSE/comenv.f BSE/corerd.f BSE/deltat.f BSE/dgcore.f BSE/evolv1.f BSE/evolv2.f \
	BSE/gntage.f BSE/hrdiag.f BSE/instar.f BSE/kick.f BSE/mix.f BSE/mlwind.f \
	BSE/mrenv.f BSE/ran3.f BSE/rl.f BSE/star.f BSE/zcnsts.f BSE/zfuncs.f input.f

FOBJ = $(FSOURCE:.f=$(BSE).o)


CSOURCE = \
	inipar/dictionary.c inipar/iniparser.c inipar/iniparser_interface.c \
	main.c

default: all

all: clean mcluster

mcluster: $(FOBJ) $(FHEAD)
	$(CC) $(FOBJ) $(CSOURCE) $(LIBS) $(FFLAGS) -o mcluster

clean:
	rm --f BSE/*.o mcluster
	rm -f *.dat *.out dat.10
