CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise1.exe : exercise1.o random.o statistics.o
	$(CC) random.o statistics.o exercise1.o -o exercise1.exe -O2 -larmadillo

exercise1.o : exercise1.cpp
	$(CC) -c exercise1.cpp -o exercise1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

statistics.o : statistics.cpp statistics.h
	$(CC) -c statistics.cpp -o statistics.o $(CFLAGS)

clean :
	rm *.o *.exe
