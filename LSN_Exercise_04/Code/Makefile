CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise_1.exe : exercise_1.o random.o statistics.o
	$(CC) random.o statistics.o exercise_1.o -o exercise_1.exe

exercise_1.o : exercise_1.cpp
	$(CC) -c exercise_1.cpp -o exercise_1.o $(CFLAGS)

exercise_2.exe : exercise_2.o random.o statistics.o
	$(CC) random.o statistics.o exercise_2.o -o exercise_2.exe

exercise_2.o : exercise_2.cpp
	$(CC) -c exercise_2.cpp -o exercise_2.o $(CFLAGS)

exercise_3.exe : exercise_3.o random.o statistics.o
	$(CC) random.o statistics.o exercise_3.o -o exercise_3.exe

exercise_3.o : exercise_3.cpp
	$(CC) -c exercise_3.cpp -o exercise_3.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

statistics.o : statistics.cpp statistics.h
	$(CC) -c statistics.cpp -o statistics.o $(CFLAGS)

clean :
	rm *.o *.exe
