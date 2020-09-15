CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type
UTILITIES_DIR = ./thirdUtils
BIO_DIR = ./bioUtils
BAM_DIR = ./thirdUtils/BamTools
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/include/api \
           -I$(BIO_DIR)
BIO_LIBS = -L$(UTILITIES_DIR)/BamTools/lib/ -lbamtools \
             -L$(BIO_DIR)/ -lbiotools

all:
	cd $(BAM_DIR); make api; make
	cd $(BIO_DIR); make
	make endSeeker

clean:
	cd $(BAM_DIR); make clean_api
	cd $(BIO_DIR); make clean
	rm -f *.o

endSeeker: endSeeker.o endSeekerMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o endSeeker endSeekerMain.o endSeeker.o $(BIO_LIBS) $(LIBS)

endSeeker.o: endSeeker.cpp endSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c endSeeker.cpp
  
endSeekerMain.o: endSeekerMain.cpp endSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c endSeekerMain.cpp
