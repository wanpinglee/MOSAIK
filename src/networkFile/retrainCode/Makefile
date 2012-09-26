# ==========================
# MOSAIK Makefile
# (c) 2012 Wan-Ping Lee
# ==========================

# define our source subdirectories
SUBDIRS = attachXC trainNetwork

all:
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

.PHONY: all

clean:
	@for dir in $(SUBDIRS); do \
		$(MAKE) clean --no-print-directory -C $$dir; \
		echo ""; \
	done

.PHONY: clean
