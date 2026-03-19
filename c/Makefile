# Franck-Condon I-V Simulation — C+GSL Implementation
# Build: make
# Run:   ./c_benchmark quick
#        ./c_benchmark default

CC      ?= gcc
CFLAGS  ?= -std=c99 -O2 -flto -march=native -Wall -Wextra -Wno-unused-parameter

# GSL flags via gsl-config (works on macOS Homebrew and Linux)
GSL_CFLAGS := $(shell gsl-config --cflags 2>/dev/null)
GSL_LIBS   := $(shell gsl-config --libs 2>/dev/null)

# Fallback if gsl-config not found
ifeq ($(GSL_CFLAGS),)
  GSL_CFLAGS := -I/opt/homebrew/include -I/usr/local/include
  GSL_LIBS   := -L/opt/homebrew/lib -L/usr/local/lib -lgsl -lgslcblas
endif

# OpenMP detection
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
  LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(LIBOMP_PREFIX),)
    OMP_CFLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    OMP_LIBS   := -L$(LIBOMP_PREFIX)/lib -lomp
  endif
else
  OMP_CFLAGS := -fopenmp
  OMP_LIBS   := -fopenmp
endif

CFLAGS  += $(OMP_CFLAGS) $(GSL_CFLAGS)
LDFLAGS += $(OMP_LIBS) $(GSL_LIBS) -lm

# Source files
SRC_DIR = src
SRCS = $(SRC_DIR)/laguerre.c \
       $(SRC_DIR)/fc_matrix.c \
       $(SRC_DIR)/fermi_bose.c \
       $(SRC_DIR)/digamma.c \
       $(SRC_DIR)/regularized.c \
       $(SRC_DIR)/cotunneling.c \
       $(SRC_DIR)/rate.c \
       $(SRC_DIR)/matrix.c \
       $(SRC_DIR)/solver.c \
       $(SRC_DIR)/current.c \
       $(SRC_DIR)/simulate.c \
       $(SRC_DIR)/json_io.c \
       $(SRC_DIR)/plotting.c \
       cJSON.c

OBJS = $(SRCS:.c=.o)
TARGET = c_benchmark

.PHONY: all clean

all: $(TARGET)

$(TARGET): main.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

main.o: main.c $(SRC_DIR)/simulate.h $(SRC_DIR)/json_io.h $(SRC_DIR)/plotting.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

cJSON.o: cJSON.c cJSON.h
	$(CC) $(CFLAGS) -Wno-unused-function -c -o $@ $<

clean:
	rm -f $(TARGET) main.o $(OBJS) cJSON.o
	rm -f $(SRC_DIR)/*.o
