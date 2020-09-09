GCC = R CMD SHLIB

SOURCE = collapsed_lpcm.c
TARGET = collpcm.so
FLAGS = 
LIBS = 
MYLIBS = nrutil.c cat.c GaussianMixtureModel.c NetworkLib.c Polya_gamma_rng.c

all:
	$(GCC) -o $(TARGET) $(FLAGS) $(MYLIBS) $(SOURCE) $(LIBS)
	

clean:
	rm *.o *.so
	

