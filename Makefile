CXXFLAGS =	-O3 -g -Wall -fmessage-length=0 -fopenmp

OBJS =		verletTestArray.o

LIBS =

TARGET =	verletTestArray

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(CXXFLAGS) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
