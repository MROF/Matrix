CXX = g++
CXXFLAGS = -g -Wall

all:

nw-sse3: nw_sse3.cc
	$(CXX) $(CXXFLAGS) -msse3 $+ -o $@

nw-avx: nw_avx.cc
	$(CXX) $(CXXFLAGS) -mavx $+ -o $@

nw-avx2: nw_avx2.cc
	$(CXX) $(CXXFLAGS) -mavx2 $+ -o $@

clean:
	rm nw-sse3 nw-avx nw-avx2






