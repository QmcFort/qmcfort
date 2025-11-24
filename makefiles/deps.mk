# Control dependencies =========================================================
DEPS              ?= 0
DEPS_ON           := 1
DEPS_FILE         := $(BUILD_DIR)/deps.mk
MAKEDEPS_SCRIPT   := $(MK_DIR)/tools/makedeps.py


# Only include and require deps when DEPS defined ==============================
ifeq ($(DEPS),$(DEPS_ON))
  -include $(DEPS_FILE)
  qmcfort: deps
  ut: deps
endif


.PHONY: deps deps-clean


deps: $(DEPS_FILE) deps-clean
	@:


$(DEPS_FILE):
	python3 $(MAKEDEPS_SCRIPT) $(SRC_DIR) $(SRC_OBJ_DIR) $(TEST_DIR) $(TEST_OBJ_DIR)
	@mkdir -p $(BUILD_DIR)
	@mv $(notdir $(DEPS_FILE)) $(DEPS_FILE)


deps-clean:
	rm -f $(DEPS_FILE)