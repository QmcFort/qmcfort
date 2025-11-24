# QmcFort source ===============================================================
SRC_DIR           := src
SRC_MAIN          := $(SRC_DIR)/main.f90


# QmcFort build tree ===========================================================
SRC_BUILD_DIR     := $(BUILD_DIR)/main
SRC_OBJ_DIR       := $(SRC_BUILD_DIR)/obj
SRC_MOD_DIR       := $(SRC_BUILD_DIR)/mod
SRC_BIN           := $(BIN_DIR)/qmcfort


# Define object files using source files SRC ===================================
SRC_OBJ           := $(patsubst $(SRC_DIR)/%.f90,$(SRC_OBJ_DIR)/%.o,$(SRC))


# Fortran compiler flags =======================================================
FFLAGS            += $(OPENMP) -O3
#FFLAGS           += $(OPENMP) $(DEBUG)


# Fortran preprocessor flags ===================================================
FPP               += -DMPI
FPP               += -DTOOLS_DIR="'$(TOOLS_DIR)'"
FPP               += -DVERSION=$(GIT_VERSION)


.PHONY: qmcfort clean


# Default target ===============================================================
qmcfort: $(SRC_BIN)
	@mkdir -p bin
	cp $(SRC_BIN) bin/$(notdir $(SRC_BIN))


# Clean build dir ==============================================================
clean:
	rm -rf $(SRC_BUILD_DIR)


# Compile object files (<fname>.f90 --> <fname>.o) =============================
$(SRC_OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(SRC_OBJ_DIR) $(SRC_MOD_DIR)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $(FPP) $(IMOD)$(SRC_MOD_DIR) $(MODFLAG)$(SRC_MOD_DIR) -c $< -o $@
	@echo ""


# Link object files into a final binary ========================================
$(SRC_BIN): $(SRC_OBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $(SRC_OBJ) $(LDFLAGS) $(LIBS)
	@echo ""


# Create necessary directories =================================================
$(SRC_OBJ_DIR) $(SRC_MOD_DIR) $(BIN_DIR) $(SRC_BUILD_DIR):
	@mkdir -p $@