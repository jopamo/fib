NAME=fib
LOWERNAME=${NAME}
OBJS = ${NAME}.o
CFLAGS= -O3 -Wall -Werror
CC = gcc
LIBS = -lm

${NAME}: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -o ${NAME} ${LIBS}

chCount.o: chCount.c chCount.h chConvert.h
	${CC} ${CFLAGS} -c chCount.c

chConvert.o: chConvert.c chConvert.h
	${CC} ${CFLAGS} -c chConvert.c

clean:
	rm -f *.o ${NAME}
