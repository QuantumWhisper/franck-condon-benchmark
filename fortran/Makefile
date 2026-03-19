# Franck-Condon I-V Simulation — Fortran Implementation
# Build: make
# Run:   ./fortran_benchmark quick
#        ./fortran_benchmark default

FC      := gfortran
FFLAGS  := -O3 -march=native -flto -funroll-loops -std=f2008 -fall-intrinsics -fopenmp
LDFLAGS := -llapack -lblas

# macOS Accelerate framework (fallback if LAPACK not found separately)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LDFLAGS = -framework Accelerate
endif

SRC_DIR = src
MODDIR  = mod

# Source files in compilation order (module dependencies)
SRCS = $(SRC_DIR)/constants.f90 \
       $(SRC_DIR)/laguerre.f90 \
       $(SRC_DIR)/fermi_bose.f90 \
       $(SRC_DIR)/digamma.f90 \
       $(SRC_DIR)/fc_matrix.f90 \
       $(SRC_DIR)/solver.f90 \
       $(SRC_DIR)/regularized.f90 \
       $(SRC_DIR)/cotunneling.f90 \
       $(SRC_DIR)/rate.f90 \
       $(SRC_DIR)/matrix.f90 \
       $(SRC_DIR)/current.f90 \
       $(SRC_DIR)/simulate.f90 \
       $(SRC_DIR)/json_io.f90 \
       $(SRC_DIR)/plotting.f90

OBJS = $(SRCS:.f90=.o)
MAIN_OBJ = main.o
TARGET = fortran_benchmark

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Pattern rule: compile .f90 to .o with module output to mod/
$(SRC_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c -o $@ $<

main.o: main.f90 $(OBJS)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c -o $@ $<

clean:
	rm -f $(TARGET) $(MAIN_OBJ) $(OBJS)
	rm -rf $(MODDIR)
