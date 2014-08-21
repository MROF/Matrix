CXX = g++
CXXFLAGS = -g -Wall

all: matrix-naive matrix-sse matrix-avx

matrix-naive: final_matrix_multiplication.cc
	$(CXX) $(CXXFLAGS) $+ -o $@

matrix-sse: final_matrix_multiplication.cc
	$(CXX) $(CXXFLAGS) -msse3 -D__SSE $+ -o $@

matrix-avx: final_matrix_multiplication.cc
	$(CXX) $(CXXFLAGS) -mavx -D__AVX $+ -o $@

clean:
	rm matrix-naive matrix-sse matrix-avx
