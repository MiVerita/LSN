CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o System.o
	$(CC) random.o System.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
system.o: System.cpp System.h
	$(CC) -c System.cpp -o System.o $(CFLAGS)	
clean :
	rm *.o main.exe seed.out
	