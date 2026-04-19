# ============================================================
#  Makefile for MCQsim25
#  Works on: macOS (Homebrew), Linux, and HPC clusters
#
#  Usage:
#    make          -> auto-detect platform and build
#    make mac      -> build explicitly for macOS (Homebrew)
#    make linux    -> build explicitly for Linux
#    make hpc      -> build for HPC (edit flags below if needed)
#    make clean    -> remove binary and build artefacts
#
#  HPC usage (SLURM example):
#    module load gcc openmpi openblas
#    make hpc
#    mpirun -np 32 ./MCQsim25
# ============================================================

TARGET = MCQsim25

SRCS = MCQsim25_Main.c \
       MCQsim25_SupportFuncs.c \
       MCQsim25_LoadParaInput.c \
       MCQsim25_LoadRghStghBnd.c \
       MCQsim25_AssignRefValues.c \
       MCQsim25_KhGeneration.c \
       MCQsim25_TectLoading.c \
       MCQsim_StrainHS.c \
       MCQsim_StrainFS.c

CFLAGS  = -Wall -O3
LDFLAGS = -lm

UNAME := $(shell uname -s)
ARCH  := $(shell uname -m)

ifeq ($(UNAME), Darwin)
HOMEBREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /opt/homebrew)
OPENBLAS_PREFIX := $(shell brew --prefix openblas 2>/dev/null || echo $(HOMEBREW_PREFIX)/opt/openblas)
CC       = mpicc
CFLAGS  += -I$(OPENBLAS_PREFIX)/include
LDFLAGS += -L$(OPENBLAS_PREFIX)/lib -lopenblas -lmpich
CLEAN_EXTRA = rm -rf $(TARGET).dSYM
else
CC       = mpicc
LDFLAGS += -lopenblas -lmpich
CLEAN_EXTRA = true
endif

.PHONY: all mac linux hpc clean

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SRCS)
	$(CLEAN_EXTRA)
	@echo "Built: $(TARGET) [$(UNAME)/$(ARCH)]"
	@echo "Run:   mpirun -np <N> ./$(TARGET)"

mac:
	$(eval OPENBLAS_PREFIX := $(shell brew --prefix openblas 2>/dev/null || echo /opt/homebrew/opt/openblas))
	mpicc -I$(OPENBLAS_PREFIX)/include \
	      -L$(OPENBLAS_PREFIX)/lib \
	      -lopenblas -lmpich -lm \
	      -Wall -O3 \
	      -o $(TARGET) $(SRCS)
	@rm -rf $(TARGET).dSYM
	@echo "Built: $(TARGET) [macOS]"

linux:
	mpicc -lopenblas -lmpich -lm \
	      -Wall -O3 \
	      -o $(TARGET) $(SRCS)
	@echo "Built: $(TARGET) [Linux]"

hpc:
	mpicc -lopenblas -lmpich -lm \
	      -Wall -O3 \
	      -o $(TARGET) $(SRCS)
	@echo "Built: $(TARGET) [HPC]"
	@echo "Run:   mpirun -np <N> ./$(TARGET)"

clean:
	rm -f $(TARGET) *.o
	@rm -rf $(TARGET).dSYM
	@echo "Cleaned."
