CC = g++
CFLAGS = -Wall -O3 --std=c++11

main : main.o
	$(CC) main.o -o main
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main *.txt
