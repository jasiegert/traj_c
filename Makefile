# Reference: https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
# Define compiler CC, flags CFLAGS, linked libraries LIBS, source files SRCS, resulting objects OBJS and executable MAIN
CC = gcc
CFLAGS = -Wextra -Wall -O3 -march=native
LIBS = -lm  # links math.h
SRCDIR = src/
SRCFILES = dostuff.c trajec_io/read_trajec.c trajec_io/chemistry.c calc/mathtools.c calc/msd.c calc/rdf.c calc/oacf.c kissFFT/kiss_fft.c kissFFT/kiss_fftr.c

HEADERS = $(SRCFILES:%.c=$(SRCDIR)%.h)
INCLUDE = $(SRCDIR:%=-I%)
SRCS = $(SRCFILES:%=$(SRCDIR)%)
OBJS = $(SRCS:%.c=%.o)  # .o-files generated from .c-files in SRCS
MAIN = traj_analyzer

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

manual: docs/html/index.html

docs/html/index.html: $(SRCS) $(HEADERS) .Doxyfile src/doxygen_modules.h README.md docs/THEORY.md
	doxygen .Doxyfile

openmanual:
	see docs/html/index.html
