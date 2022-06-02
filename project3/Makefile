all: ga

ga: ga.cpp
	g++ -std=c++14 -o ga -O3 ga.cpp

run: ga
	./ga < maxcut.in > maxcut.out

clean:
	rm ga maxcut.in maxcut.out
