CC = g++
CFLAGS = -Wall -O3 --std=c++11

main : main.o random.o point.o functions.o
	$(CC) random.o point.o main.o functions.o -o main
rmax : rmax.o random.o point.o functions.o
	$(CC) random.o point.o rmax.o functions.o -o rmax
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
rmax.o : rmax.cpp
	$(CC) -c rmax.cpp -o rmax.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
point.o: point.cpp point.h
	$(CC) -c point.cpp -o point.o $(CFLAGS)
functions.o: functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.txt
