CC	:= mpic++   # use cc on okeanos
CFLAGS	:= -O3 -c --std=c++2a -fopenmp 
LFLAGS	:= -O3
# Add new targets below:
zadanie2	:= zadanie2.exe



zadanie2: zadanie2.cpp
	${CC} -Wall -fopenmp -O3 -c --std=c++2a --std=c++2a zadanie2.cpp data_source.cpp
 
#%.exe : %.o
#	$(CC) $(LFLAGS) -o $@ $<


#%.o : %.c
#	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o *.out *.err $(ALL)

