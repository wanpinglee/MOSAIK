# =============================
# MOSAIK Common Source Makefile
# (c) 2008 Michael Stromberg
# =============================

# define our build subdirectories
SUBDIRS = AssemblyFormats DataStructures ExternalReadFormats MosaikReadFormat PairwiseAlignment Utilities

all:
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

.PHONY: all
