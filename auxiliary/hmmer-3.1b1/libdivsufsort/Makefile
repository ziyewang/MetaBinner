# Makefile for fm-index

# VPATH and shell configuration
top_srcdir = ..
srcdir     = .

SHELL      = /bin/sh

# sources
OBJS				= divsufsort.o 
TARGET				= libdivsufsort.a
MAKEFILE			= Makefile

# options
CC        = gcc -std=gnu99
CFLAGS    = -I/usr/usc/gnu/mpfr/3.1.2/include -L/usr/usc/gnu/mpfr/3.1.2/lib -I/usr/usc/gnu/mpc/1.0.1/include -L/usr/usc/gnu/mpc/1.0.1/lib -I/usr/usc/gnu/gmp/5.1.2/include -L/usr/usc/gnu/gmp/5.1.2/lib -I/usr/usc/gnu/cloog/0.18.0/include -L/usr/usc/gnu/cloog/0.18.0/lib -pthread -fPIC
CPPFLAGS  = -I/usr/usc/gnu/mpfr/3.1.2/include -L/usr/usc/gnu/mpfr/3.1.2/lib -I/usr/usc/gnu/mpc/1.0.1/include -L/usr/usc/gnu/mpc/1.0.1/lib -I/usr/usc/gnu/gmp/5.1.2/include -L/usr/usc/gnu/gmp/5.1.2/lib -I/usr/usc/gnu/cloog/0.18.0/include -L/usr/usc/gnu/cloog/0.18.0/lib
MPILIBS   = 
AR        = /usr/bin/ar rc
RANLIB    = ranlib
INSTALL   = /usr/bin/install -c

#CC					= gcc 
#CFLAGS  	 		= -g -O3 -fomit-frame-pointer -fstrict-aliasing 
#CPPFLAGS			= -Wall 
#AR       = /usr/bin/ar rc
#RANLIB   = ranlib


# beautification magic stolen from git (added within hmmer source)
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_AR      = @echo '    ' AR $@;
endif


# targets
all: $(TARGET) 
$(TARGET): $(OBJS)

.c.o:
	${QUIET_CC}${CC} -I. ${CFLAGS} ${CPPFLAGS} -o $@ -c $<


libdivsufsort.a: $(OBJS)
	${QUIET_AR}${AR} libdivsufsort.a $(OBJS)
	@${RANLIB} libdivsufsort.a
	@chmod 644 libdivsufsort.a


clean:
	$(RM) $(TARGET) $(OBJS)

distclean:
	$(RM) $(TARGET) $(OBJS) $(MAKEFILE) divsufsort.h

# dependencies
$(OBJS): divsufsort.h $(MAKEFILE)
