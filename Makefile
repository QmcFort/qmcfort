# Entry ========================================================================

# Select toolchain: make TOOLCHAIN=ifort-mkl-intelmpi qmcfort ==================
TOOLCHAIN         ?= ifort-mkl-intelmpi


# Relevant info ================================================================
QMCFORT_ROOT      := $(abspath .)
MK_DIR            := makefiles
BUILD_DIR         ?= build/$(TOOLCHAIN)
BIN_DIR           := $(BUILD_DIR)/bin
TOOLS_DIR         := $(abspath tools)


# Fail early if toolchain file missing =========================================
TOOLCHAIN_MK      := $(MK_DIR)/toolchain/$(TOOLCHAIN).mk
ifeq ($(wildcard $(TOOLCHAIN_MK)),)
$(error Unknown TOOLCHAIN '$(TOOLCHAIN)'. Available: $(notdir $(basename $(wildcard $(MK_DIR)/toolchain/*.mk))))
endif


# Build info ===================================================================
GIT_VERSION       := "'$(shell git log -1 --pretty=format:%h)'"
COMPILER_VERSION  := $(shell $(FC) --version | head -n 1)
MPI_VERSION       := $(shell mpirun --version | head -n 1)
DATE              := $(shell git log -1 --format="%cd")


# include sub-makefiles ========================================================
include $(MK_DIR)/help.mk  	        				# 1) how to use
include $(TOOLCHAIN_MK)         					# 2) compiler & libs
include Sources.mk 									# 3) source and test files
include $(MK_DIR)/qmcfort.mk 						# 4) qmcfort build rules
include $(MK_DIR)/python_venv.mk					# 5) python virtual environment setup
include $(MK_DIR)/unit_tests.mk						# 6) unit tests build rules
include $(MK_DIR)/regression_tests.mk 				# 7) regression tests rules
include $(MK_DIR)/deps.mk                           # 8) dependencies rules
include $(MK_DIR)/basis_sets.mk						# 9) Gaussian basis sets 


# Default goal
.DEFAULT_GOAL := qmcfort

.PHONY: distclean veryclean

distclean:
	rm -rf $(BUILD_DIR)

veryclean:
	rm -rf build