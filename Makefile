CXX=g++
%.o: seven-op.h
seven-op: seven-op.o seven-op-common.o
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LDFLAGS)
