PC := $(shell hostname -s)
-include Makefile.$(PC)

OJ_VERIFY ?= oj-verify
VERIFY_JOBS ?= 2

.PHONY: verify

verify:
	$(OJ_VERIFY) run -j $(VERIFY_JOBS) $(TEST)