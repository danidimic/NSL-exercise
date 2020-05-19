CC = g++
CFLAGS = -Wall -O3 --std=c++11

variationalMC.exe : variationalMC.o random.o statistics.o
	$(CC) random.o statistics.o variationalMC.o -o variationalMC.exe -O2 -larmadillo

variationalMC.o : variationalMC.cpp
	$(CC) -c variationalMC.cpp -o variationalMC.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

statistics.o : statistics.cpp statistics.h
	$(CC) -c statistics.cpp -o statistics.o $(CFLAGS)

clean :
	rm *.o *.exe
