CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		NMR_2015revamp.o

LIBS =

TARGET =	NMR_2015revamp

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
