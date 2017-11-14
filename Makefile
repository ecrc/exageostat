GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

# Common binaries
CC             := $(CC)
LD             := $(CC)
OPTS           := -c -fPIC


flagsc := $(shell pkg-config nlopt gsl libstarpu chameleon  --cflags)
flagsl := $(shell pkg-config nlopt gsl libstarpu chameleon  --libs)

#user flags
CCFLAGS   := -O3  -w -Ofast -Wall $(flagsc) -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS   := -O3  -w -Ofast -lstarpu-1.2   -lchameleon  -lchameleon_starpu -lcoreblas -lstdc++  -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl $(flagsl)
# Extra user flags
EXTRA_CCFLAGS = -I$(PLASMADIR)/include -I./src/include/ -I./exageostat_exact/core/include/ -I./exageostat_exact/runtime/starpu/include/ -I./misc/include/ -I ./exageostat_exact/src/include/ -I ./r-wrappers/include -I/home/abdullsm/codes/chameleon/coreblas/include/coreblas
# Optimisation flags
#CXXFLAGS := $(CCFLAGS)

# DEBUG flags
#CCFLAGS   += -g

# PIC flags
#CCFLAGS   += -fPIC

MAIN =  ./examples/zgen_mle_test.c
MAIN2 =  ./examples/r_zgen_mle_test.c
SOURCE_FILES = $(filter-out $(MAIN), $(wildcard ./src/compute/*.c ./exageostat_exact/core/compute/*.c ./exageostat_exact/runtime/starpu/codelets/*.c ./misc/compute/*.c    ./exageostat_exact/src/compute/*.c ./r-wrappers/compute/*.c))

OBJ_FILES = $(patsubst %.c,%.o,$(SOURCE_FILES))

OBJECTS = $(patsubst %.o)
exageostat.so:    $(OBJ_FILES)
	R CMD SHLIB -o exageostat.so $(OBJ_FILES)  -Ofast -lstarpu-1.2   -lchameleon  -lchameleon_starpu -lcoreblas -lstdc++  -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl $(flagsl)

#.c.o; $(CC) $(OPTS) -c $<

#####################
# default make rule #
#####################

ARGET = ./examples/zgen_mle_test
TARGET2 =  ./examples/r_zgen_mle_test

.PHONY: all
all: $(TARGET)
all: $(TARGET2)


$(TARGET): $(MAIN) $(SOURCE_FILES)
	$(LD) $(OPTS) $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS)
%.o: %.c
	$(CC) $(OPTS) -c $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@

$(TARGET2): $(MAIN2) $(SOURCE_FILES)
	$(LD) $(OPTS) $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS)
%.o: %.c
	$(CC) $(OPTS) -c $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -rf *.o *.so
	rm -f $(TARGET)
	find . -name "*.o" -type f -delete
