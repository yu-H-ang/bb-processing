## DEPENDENCIES ##
HDF5_DIR = /opt/hdf5
CGNS_DIR = /home-4/yzhan175@jhu.edu/deps/cgns
CUSP_DIR = /home-4/yzhan175@jhu.edu/deps/cusplibrary

## COMPILERS
GCC = gcc

SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include -I $(CUSP_DIR)
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

SRCC =	main.c cgns_reader.c
EXTRA = Makefile main.h	cgns_reader.h

# compile normally:
# make -g -O2 to optimize
all: COPT += -O2
all: main

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: CUDAOPT += -DDEBUG -g
debug: main

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

main: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -lstdc++ -o main

clean:
	rm -f $(SRC_DIR)/*.o main
