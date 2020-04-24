CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: exercise1.exe exercise2.exe exercise3.exe

exercise1.exe : exercise1.o random.o statistics.o
	$(CC) random.o statistics.o exercise1.o -o exercise1.exe

exercise1.o : exercise1.cpp
	$(CC) -c exercise1.cpp -o exercise1.o $(CFLAGS)

exercise2.exe : exercise2.o random.o statistics.o
	$(CC) random.o statistics.o exercise2.o -o exercise2.exe

exercise2.o : exercise2.cpp
	$(CC) -c exercise2.cpp -o exercise2.o $(CFLAGS)

exercise3.exe : exercise3.o random.o statistics.o
	$(CC) random.o statistics.o exercise3.o -o exercise3.exe

exercise3.o : exercise3.cpp
	$(CC) -c exercise3.cpp -o exercise3.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

statistics.o : statistics.cpp statistics.h
	$(CC) -c statistics.cpp -o statistics.o $(CFLAGS)

clean :
	rm *.o *.exe
