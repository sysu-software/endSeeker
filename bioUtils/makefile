UTILITIES_DIR = ../thirdUtils

INCLUDES = -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/include/api

SOURCES = samFile.cpp bedFile.cpp bioUtils.cpp fasta.cpp faiFile.cpp

OBJECTS= $(SOURCES:.cpp=.o)

biolib.a: $(OBJECTS) $(SOURCES)
	ar rcus libbiotools.a $(OBJECTS)

$(OBJECTS): $(SOURCES)
	@$(CXX) -c $(*D)/$(*F).cpp $(INCLUDES) -o $(*F).o

clean:
	rm -f *.o libbiotools.a
