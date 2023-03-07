CC = gcc
LD = gcc
CFLAGS =  
LDFLAGS = -shared -llapacke -llapack -lblas -lm
OBJFILES = scg.o cg.o chol.o
SHARED = iter
SOURCE = scg.c cg.c chol.c

all: $(SHARED)

$(SHARED): $(OBJFILES)
	$(LD) -o $(SHARED) $(OBJFILES) $(LDFLAGS)

$(OBJFILES): $(SOURCE)
	$(CC) $(CFLAGS) $(INC) -c $(SOURCE) 

clean:
	rm -f $(OBJFILES) $(SHARED) *.jpg *.pdf
