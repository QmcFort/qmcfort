# Gaussian basis sets variables ================================================
BS_DIR            ?= data/basis_sets


.PHONY: basis_sets


basis_sets: venv-bse
	$(VENV_PYBIN) $(TOOLS_DIR)/update_basis_sets.py $(BS_DIR)