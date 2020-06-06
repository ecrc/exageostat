.PHONY: all clean
all:
	(cd .. && make VERBOSE=1 && cp ./lib*.so ./src/exageostat.so)

