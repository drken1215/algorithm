PC := $(shell hostname -s)
-include Makefile.$(PC)

OJ_VERIFY ?= oj-verify
VERIFY_CXX ?= g++
VERIFY_JOBS ?= 2

.PHONY: configure verify add docs

configure:
	python3 .verify-helper/tools/setup_includes.py --cxx "$(VERIFY_CXX)"

verify: configure
	$(OJ_VERIFY) run -j $(VERIFY_JOBS) $(TEST)

add:
	python3 .verify-helper/tools/add_library.py --lib "$(LIB)"

docs: configure
	python3 .verify-helper/tools/update_readme_links.py
	cp README.md .verify-helper/docs/index.md
	python3 .verify-helper/tools/docs_with_includes.py \
		--cxx "$(VERIFY_CXX)" -j $(VERIFY_JOBS)
