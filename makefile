PC := $(shell hostname -s)
-include Makefile.$(PC)

OJ_VERIFY ?= oj-verify
VERIFY_CXX ?= g++
VERIFY_JOBS ?= 2

.PHONY: configure readme docs verify add add_atcoder add_yukicoder add_aoj add_codeforces add_library_checker

add_atcoder: TEST_DIR = tests/AtCoder
add_atcoder: TEST_EXT = .cpp
add_atcoder: add_test

add_yukicoder: TEST_DIR = tests/yukicoder
add_yukicoder: TEST_EXT = .test.cpp
add_yukicoder: add_test

add_aoj: TEST_DIR = tests/AOJ
add_aoj: TEST_EXT = .test.cpp
add_aoj: add_test

add_codeforces: TEST_DIR = tests/Codeforces
add_codeforces: TEST_EXT = .cpp
add_codeforces: add_test

add_library_checker: TEST_DIR = tests/LibraryChecker
add_library_checker: TEST_EXT = .test.cpp
add_library_checker: add_test

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

add_test:
	@if [ -z "$(TITLE)" ]; then \
		echo "Error: TITLE is required."; \
		echo "Usage: make add_atcoder TITLE=ABC_250_G_dual_slope_trick"; \
		exit 1; \
	fi
	@mkdir -p "$(TEST_DIR)"
	@cpp="$(TEST_DIR)/$(TITLE)$(TEST_EXT)"; \
	md="$(TEST_DIR)/$(TITLE)$(patsubst %.cpp,%.md,$(TEST_EXT))"; \
	if [ -e "$$cpp" ] || [ -e "$$md" ]; then \
		echo "Error: File already exists."; \
		exit 1; \
	fi; \
	touch "$$cpp"; \
	printf '%s\n' \
		'---' \
		'title: ' \
		'documentation_of: $(TITLE)$(TEST_EXT)' \
		'---' \
		'' \
		'## 問題概要' \
		'' \
		'## 制約' \
		'' \
		> "$$md"; \
	echo "Created: $$cpp"; \
	echo "Created: $$md"
