# Basic Python variables =======================================================
PY                := python3
PYBIN             := python3
PIP               := pip3
VENV_DIR          ?= .venv
VENV_PYBIN        := $(VENV_DIR)/bin/python
VENV_PIP          := $(VENV_DIR)/bin/pip


# QmcFortPy variables ==========================================================
QMCFORTPY_DIR     := qmcfortpy


# Separate stamp for different tools ===========================================
BSE_OK            := $(VENV_DIR)/.venv-bse.ok
BSE_REQ           := $(TOOLS_DIR)/requirements-bse.txt
QMCFORTPY_OK      := $(VENV_DIR)/.venv-qmcfortpy.ok
QCMFORTPY_REQ     := $(TOOLS_DIR)/requirements-qmcfortpy.txt


.PHONY: venv venv-clean


venv: $(VENV_DIR)/.ok
$(VENV_DIR)/.ok: $(TOOLS_DIR)/requirements-base.txt
	$(PYBIN) -m venv $(VENV_DIR)
	$(VENV_PIP) install --upgrade pip
	$(VENV_PIP) install -r $<
	touch $@


venv-clean:
	rm -rf $(VENV_DIR)


.PHONY: venv-all venv-bse venv-qmcfortpy


venv-all: venv-bse venv-qmcfortpy
	@:


venv-bse: $(BSE_OK)
$(BSE_OK): $(BSE_REQ) $(VENV_DIR)/.ok
	$(VENV_PIP) install -r $(BSE_REQ)
	touch $(BSE_OK)


venv-qmcfortpy: $(QMCFORTPY_OK)
$(QMCFORTPY_OK): $(QCMFORTPY_REQ) $(VENV_DIR)/.ok
	$(VENV_PIP) install -r $(QCMFORTPY_REQ)
	$(VENV_PIP) install -e $(QMCFORTPY_DIR)/.
	touch $(QMCFORTPY_OK)