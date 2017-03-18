CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -std=c++11 -fopenmp

BASE_PATH = files/impl/BasePath.o
CREATE_FOLDER = files/impl/CreateFolder.o
GENERATOR = generator/Generator.o

SORTER = sorter/Sorter.o

OBJS =		Yvenss.o $(BASE_PATH) $(CREATE_FOLDER) $(GENERATOR) $(SORTER)

LIBS = -fopenmp

TARGET =	yvenss

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
