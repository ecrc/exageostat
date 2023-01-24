# Common binaries
CC             = gcc
LD             = gcc

#user flags -std=c++11
CCFLAGS   := -O3  -w -Ofast -Wall   `pkg-config plasma --cflags`   
FCFLAGS   := -O3  -w -Ofast
LDFLAGS   := -O3  -w -Ofast `pkg-config plasma --libs`  -lmkl_gf_lp64 -lmpi -lmagma -lstarpu-1.2 -lstarpumpi-1.2   -lchameleon  -lchameleon_starpu  /home/abdullsm/software/spack/opt/spack/linux-x86_64/gcc-4.8/chameleon-0.9.1-he4lewz2dezrda33megkd6yy3rfzfsjn/lib/libcoreblas.so  /home/abdullsm/software/nlopt_install/lib/libnlopt.a /home/abdullsm/software/gsl/lib/libgsl.so   -Wl,--no-as-needed -lpthread -lm -ldl  
# With GNU compilers
#LDFLAGS   := -O3 -L$(PLASMADIR)/lib -lplasma   -lcoreblas -lplasma -lquark  -llapacke -L/opt/share/INTEL/v15.app/composer_xe_2015.0.090/mkl/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--no-as-needed -lpthread -lm -ldl


EXTRA_CCFLAGS = -I$(PLASMADIR)/include -I.  


# Optimisation flags
CXXFLAGS := $(CCFLAGS)

# DEBUG flags
#CCFLAGS   += -g

# PIC flags
#CCFLAGS   += -fPIC 

MAIN = MLE.c
#SOURCE_FILES = $(filter-out $(MAIN), $(wildcard *.c))
OBJ_FILES = $(patsubst %.c,%.o,$(SOURCE_FILES))

#####################
# default make rule #
#####################

TARGET = MLE

.PHONY: all
all: $(TARGET)

$(TARGET): MLE.o
	$(LD) $(FCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS)

%.o: %.c
	$(CC) -c $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o
