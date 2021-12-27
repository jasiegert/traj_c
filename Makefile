# Reference: https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
# Define compiler CC, flags CFLAGS, linked libraries LIBS, source files SRCS, resulting objects OBJS and executable MAIN
CC = gcc
CFLAGS = -Wall -O3 -march=native
LIBS = -lm  # links math.h
SRCDIR = src/
SRCFILES = trajec_io/read_trajec.c trajec_io/chemistry.c calc/mathtools.c calc/msd.c calc/rdf.c calc/oacf.c kissFFT/kiss_fft.c kissFFT/kiss_fftr.c
MAINSRC = dostuff.c     # special treatment, because it doesn't have a header file

HEADERS = $(SRCFILES:%.c=$(SRCDIR)%.h)
INCLUDE = $(SRCDIR:%=-I%)
SRCS = $(SRCFILES:%=$(SRCDIR)%) $(MAINSRC:%=$(SRCDIR)%)
OBJS = $(SRCS:%.c=%.o)  # .o-files generated from .c-files in SRCS
MAIN = dostuff.out

# Default action is to compile executable MAIN
all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(INCLUDE) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

# Generate .o-files from corresponding .c-files where needed
.c.o:
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

# Delete all .o-files named in OBJS with "make clean"
clean:
	$(RM) $(OBJS)

remove:
	$(RM) $(OBJS) $(MAIN)

manual: manual/html/index.html

manual/html/index.html: $(MAIN) $(SRCS) $(HEADERS)
	doxygen Doxyfile
