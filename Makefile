# Common binaries
CC             := $(CC)
LD             := $(CC)
OPTS           := -c -fPIC -DEXAGEOSTAT_USE_HICMA


flagsc := $(shell pkg-config nlopt gsl libstarpu chameleon hicma  starsh --cflags)
flagsl := $(shell pkg-config nlopt gsl libstarpu chameleon hicma  starsh --libs)
flagss := $(shell pkg-config  gsl libstarpu chameleon netcdf  --libs)
#user flags
CCFLAGS   := -O3  -w -Ofast -Wall $(flagsc) -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS   := -O3  -w -Ofast -lstarpu-1.2   -lchameleon  -lchameleon_starpu -lhicma -lcoreblas -lstdc++  -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl $(flagsl)
# Extra user flags
EXTRA_CCFLAGS = -I./include/ -I./src/include/ -I./exageostat_exact/core/include/ -I./exageostat_exact/runtime/starpu/include/ -I./misc/include/ -I ./exageostat_exact/src/include/ -I ./r-wrappers/include  -I./exageostat_approx/runtime/starpu/include/ -I./exageostat_approx/src/include/
# Optimisation flags
#CXXFLAGS := $(CCFLAGS)

# DEBUG flags
#CCFLAGS   += -g

# PIC flags
#CCFLAGS   += -fPIC

MAIN =  ./examples/zgen_mle_testr.c
SOURCE_FILES = $(filter-out $(MAIN), $(wildcard ./src/compute/*.c ./exageostat_exact/core/compute/*.c ./exageostat_exact/src/compute/*.c ./exageostat_exact/runtime/starpu/codelets/*.c  ./misc/compute/flat_file.c  ./misc/compute/MLE_misc.c    ./r-wrappers/compute/*.c ./exageostat_approx/runtime/starpu/codelets/*.c ./exageostat_approx/src/compute/*.c))

OBJ_FILES = $(patsubst %.c,%.o,$(SOURCE_FILES))

OBJECTS = $(patsubst %.o)
exageostat.so:    $(OBJ_FILES)
	R CMD SHLIB -o exageostat.so $(OBJ_FILES)  -Ofast -lchameleon_starpu -lcoreblas -lstdc++  -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl $(flagsl)

#.c.o; $(CC) $(OPTS) -c $<

#####################
# default make rule #
#####################

TARGET = ./examples/zgen_mle_testr

.PHONY: all
all: $(TARGET)


$(TARGET): $(MAIN) $(SOURCE_FILES)
	$(LD) $(OPTS) $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS)
%.o: %.c
	$(CC) $(OPTS) -c $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -rf *.o *.so
	rm -f $(TARGET)
	find . -name "*.o" -type f -delete
