# Define compiler CC, flags CFLAGS, linked libraries LIBS, source files SRCS, resulting objects OBJS and executable MAIN
CC = clang
CFLAGS = -Wall
LIBS = -lm  # links math.h
SRCS = dostuff.c read_trajec.c chemistry.c calc/msd.c calc/rdf.c calc/oacf.c kissFFT/kiss_fft.c kissFFT/kiss_fftr.c
OBJS = $(SRCS:.c=.o)  # .o-files generated from .c-files in SRCS
MAIN = dostuff.out

# Default action is to compile executable MAIN
all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

# Generate .o-files from corresponding .c-files where needed
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

# Delete all .o-files named in OBJS with "make clean"
clean:
	$(RM) $(OBJS)
