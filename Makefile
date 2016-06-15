all: CImg.h filters.cpp
	g++ filters.cpp -o filters -lX11 -lpthread

clean:
	rm filters
