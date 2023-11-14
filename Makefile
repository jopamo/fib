NAME=fib
LOWERNAME=${NAME}
OBJS = ${NAME}.o
CFLAGS= -O3 -Wall -Werror
CC = gcc
LIBS = -lm -lgmp -lmpfr

${NAME}: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -o ${NAME} ${LIBS}

clean:
	rm -f *.o ${NAME}
