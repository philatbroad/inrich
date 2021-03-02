CXX = g++
CFLAGS = -O2
OUTPUT = inrich
SRC = main.cpp inrichmain.cpp crandom.cpp helper.cpp loaders.cpp inputdata.cpp common.cpp interval.cpp gene.cpp snp.cpp pairintervaldist.cpp allres.cpp pathwayres.cpp generes.cpp intervalres.cpp interval2gene.cpp log.cpp aligatormain.cpp inrichmultimain.cpp
HDR = inrichmain.h crandom.h helper.h loaders.h inputdata.h common.h interval.h gene.h snp.h pairintervaldist.h allres.h pathwayres.h generes.h intervalres.h interval2gene.h log.h aligatormain.h inrichmultimain.h
OBJ = main.o inrichmain.o crandom.o helper.o loaders.o inputdata.o common.o interval.o allres.o pathwayres.o generes.o intervalres.o aligatormain.o inrichmultimain.o

all : $(OUTPUT)

$(OUTPUT) :
	$(CXX) $(CFLAGS) -o $(OUTPUT) $(OBJ) -I. 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CFLAGS) -c $*.cpp
.SUFFIXES : .cpp .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	rm -f *.o *~ $(OUTPUT)
