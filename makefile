PC := $(shell hostname -s)
-include Makefile.$(PC)

OJ_VERIFY ?= oj-verify
VERIFY_CXX ?= g++
VERIFY_JOBS ?= 2

.PHONY: configure readme docs verify add add_atcoder

configure:
	python3 .verify-helper/tools/setup_includes.py --cxx "$(VERIFY_CXX)"

readme: configure
	python3 .verify-helper/tools/update_readme_links.py

docs: readme
	cp README.md .verify-helper/docs/index.md
	python3 .verify-helper/tools/docs_with_includes.py --cxx "$(VERIFY_CXX)" -j $(VERIFY_JOBS)

verify: configure
	$(OJ_VERIFY) run -j $(VERIFY_JOBS) $(TEST)
	python3 .verify-helper/tools/sync_local_verification.py export

add:
	python3 .verify-helper/tools/add_library.py --lib "$(LIB)"

add_atcoder:
	@if [ -z "$(TITLE)" ]; then \
		echo "Error: TITLE is required."; \
		echo "Usage: make add_atcoder TITLE=EDPC_Z_monge_shortest_path_by_cht"; \
		exit 1; \
	fi
	@mkdir -p tests/AtCoder
	@if [ -e "tests/AtCoder/$(TITLE).cpp" ] || \
	    [ -e "tests/AtCoder/$(TITLE).md" ]; then \
		echo "Error: File already exists."; \
		exit 1; \
	fi
	@touch "tests/AtCoder/$(TITLE).cpp"
	@printf '%s\n' \
		'---' \
		'title: ' \
		'documentation_of: $(TITLE).cpp' \
		'---' \
		'' \
		'## 問題概要' \
		'' \
		'## 制約' \
		'' \
		> "tests/AtCoder/$(TITLE).md"
	@echo "Created: tests/AtCoder/$(TITLE).cpp"
	@echo "Created: tests/AtCoder/$(TITLE).md"

