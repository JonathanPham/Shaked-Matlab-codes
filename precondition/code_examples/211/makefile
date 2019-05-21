CXX := g++
CXXFLAGS := -O3 -Wall -Wextra -Wconversion -std=c++11
TARGET := main
OBJS := main.o CGSolver.o COO2CSR.o matvecops.o heat.o sparse.o
INCS := CGSolver.hpp COO2CSR.hpp matvecops.hpp heat.hpp spare.hpp
$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS)
	# this is a make pattern rule
# %.o: %.cpp $(INCS)
	 # $(CXX) -c -o $@ $< $(CXXFLAGS)
.PHONY: clean
clean:
	 $(RM) $(OBJS) $(TARGET)
	
