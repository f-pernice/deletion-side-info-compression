
TARGET = sim
TARGET2 = upper
TARGET3 = mitz-lb

CXXBASE = g++
CXX = $(CXXBASE) -std=c++17
CXXFLAGS = -g -ggdb -O -Wall -Werror

CPPFLAGS =
LIBS =

OBJS = sim.o
OBJS2 = upper.o
OBJS3 = mitz-lb.o

HEADERS =

all: $(TARGET) $(TARGET2) $(TARGET3)

$(OBJS): $(HEADERS)

$(TARGET): $(OBJS)
	$(CXX) -o exec_sim $(OBJS) $(LIBS)

$(TARGET2): $(OBJS2)
	$(CXX) -o exec_upp $(OBJS2) $(LIBS)

$(TARGET3): $(OBJS3)
	$(CXX) -o exec_mitz $(OBJS3) $(LIBS)

clean:
	rm -f exec_sim $(TARGET) $(LIB) $(OBJS) $(LIBOBJS) *~ .*~ _test_data*
	rm -f exec_upp $(TARGET2) $(LIB) $(OBJS2) $(LIBOBJS) *~ .*~ _test_data*
	rm -f exec_mitz $(TARGET3) $(LIB) $(OBJS3) $(LIBOBJS) *~ .*~ _test_data*

.PHONY: all clean starter
