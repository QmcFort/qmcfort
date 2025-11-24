# Unit tests sources ===========================================================
TEST_DIR          := tests/unit
TEST_MAIN         := $(TEST_DIR)/test_main.f90


# Unit tests build tree ========================================================
TEST_BUILD_DIR    := $(BUILD_DIR)/tests/unit
TEST_OBJ_DIR      := $(TEST_BUILD_DIR)/obj
TEST_MOD_DIR      := $(TEST_BUILD_DIR)/mod
TEST_BIN          := $(BIN_DIR)/qmcfort_tests


# Collect sources and objects for unit tests ==================================
SRCX              := $(filter-out $(SRC_MAIN),$(SRC))
SRCX_OBJ          := $(patsubst $(SRC_DIR)/%.f90,$(SRC_OBJ_DIR)/%.o,$(SRCX))
TEST_OBJ          := $(patsubst $(TEST_DIR)/%.f90,$(TEST_OBJ_DIR)/%.o,$(TEST_SRC))
TEST_ALL_SRC      := $(SRCX) $(TEST_SRC)
TEST_ALL_OBJ      := $(SRCX_OBJ) $(TEST_OBJ)


.PHONY: ut ut-clean

# Execute unit tests ===========================================================
ut: $(TEST_BIN)
	@mkdir -p bin
	cp $(TEST_BIN) bin/$(notdir $(TEST_BIN))
	./$(TEST_BIN)


# Clean unit tests ============================================================
ut-clean:
	rm -rf $(TEST_BUILD_DIR)


# Compile object files (<fname>.f90 --> <fname>.o) =============================
$(TEST_OBJ_DIR)/%.o: $(TEST_DIR)/%.f90 | $(TEST_OBJ_DIR) $(TEST_MOD_DIR)
	@mkdir -p $(dir $@)
	@echo ""
	$(FC) $(FFLAGS) $(FPP) $(IMOD)$(SRC_MOD_DIR) $(IMOD)$(TEST_MOD_DIR) $(MODFLAG)$(TEST_MOD_DIR) -c $< -o $@


# Link object files into a final binary ========================================
$(TEST_BIN): $(TEST_ALL_OBJ) | $(BIN_DIR)
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(TEST_ALL_OBJ) $(LDFLAGS) $(LIBS)


# Create necessary directories =================================================
$(TEST_OBJ_DIR) $(TEST_MOD_DIR) $(TEST_BUILD_DIR):
	@mkdir -p $@