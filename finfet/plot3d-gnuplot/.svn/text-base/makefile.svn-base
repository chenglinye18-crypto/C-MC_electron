CC=gcc
CFLAGS+= -g -DDEBUG

gen_data: gen_data.o utils.o
	$(CC) $(CFLAGS) $^ -o $@ 

gen_data.o: gen_data.c 
	$(CC) $(CFLAGS) -c $^ -o $@ 
utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c $< -o $@ 
clean:
	rm -fr *.o gen_data
