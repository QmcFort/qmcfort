# Default directories for reg_tests ============================================
REG_DIR := tests/regression
BIN     ?= bin/qmcfort
OUT_DIR ?= $(REG_DIR)/results
CASES   ?=


.PHONY: rt rt-update rt-clean rt-list

rt: rt-clean-report
	@bash -Eeuo pipefail -c '\
	source "$(abspath $(REG_DIR))/src/regression_test_driver.sh"; \
	run_regression_tests \
		"$(abspath $(BIN))" \
		"$(abspath $(OUT_DIR))" \
		"$(CASES)"'


rt-update:
	@bash -Eeuo pipefail -c '\
	source "$(abspath $(REG_DIR))/src/regression_test_driver.sh"; \
	update_regression_tests \
		"$(abspath $(OUT_DIR))" \
		"$(CASES)" \'


rt-clean: rt-clean-report
	rm -rf $(OUT_DIR)


rt-clean-report:
	rm -rf $(REG_DIR)/report.log


rt-list:
	@echo "Available regression tests (see tests/regression/cases/<case>):"
	@echo "==============================================================="
	@find "$(REG_DIR)/cases" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort