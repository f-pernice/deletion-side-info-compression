
TARGET = sim
TARGET2 = upper

CXXBASE = g++
CXX = $(CXXBASE) -std=c++17
CXXFLAGS = -g -ggdb -O -Wall -Werror

CPPFLAGS =
LIBS =

OBJS = sim.o
OBJS2 = upper.o

HEADERS =

all: $(TARGET) $(TARGET2)

$(OBJS): $(HEADERS)

$(TARGET): $(OBJS)
	$(CXX) -o exec_sim $(OBJS) $(LIBS)

$(TARGET2): $(OBJS2)
	$(CXX) -o exec_upp $(OBJS2) $(LIBS)


clean:
	rm -f exec_sim $(TARGET) $(LIB) $(OBJS) $(LIBOBJS) *~ .*~ _test_data*
	rm -f exec_upp $(TARGET2) $(LIB) $(OBJS2) $(LIBOBJS) *~ .*~ _test_data*

.PHONY: all clean starter
