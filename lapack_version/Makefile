
# Common binaries
CC             = gcc
LD             = gcc

#user flags
CCFLAGS   := -O3 -w -std=c++11 -fopenmp -Ofast -Wall  `pkg-config plasma --cflags`
FCFLAGS   := -O3 -w -Ofast
LDFLAGS   := -O3   -w -Ofast -L/home/abdullsm/codes/plasma-20.9.20/build/install_dir/lib -lnlopt  -lgsl  -lmkl_gf_lp64  -L/home/abdullsm/develop/hicma/chameleon/build/install_dir/lib  -llapacke  -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--no-as-needed -lmkl_core -lgomp  -lpthread -lm -ldl

# Extra user flags
EXTRA_CCFLAGS = -I/home/abdullsm/develop/hicma/chameleon/build/install_dir/include/coreblas  -I/home/abdullsm/codes/plasma-20.9.20/build/install_dir/include

# Optimisation flags
CXXFLAGS := $(CCFLAGS)

# DEBUG flags
#CCFLAGS   += -g

# PIC flags
#CCFLAGS   += -fPIC 

MAIN = LAPACKE_MLE.c
#SOURCE_FILES = $(filter-out $(MAIN), $(wildcard *.c))
OBJ_FILES = $(patsubst %.c,%.o,$(SOURCE_FILES))

#####################
# default make rule #
#####################

TARGET = LAPACKE_MLE

.PHONY: all
all: $(TARGET)

$(TARGET): LAPACKE_MLE.o
	$(LD) $(FCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS)

%.o: %.c
	$(CC) -c $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o
