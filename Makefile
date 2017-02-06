blas ?= /opt/openblas

all:
	gcc -o contactMatrix -O2 src/contactMatrix.c src/KRnorm.c -I$(blas)/include -L$(blas)/lib -lopenblas -lm