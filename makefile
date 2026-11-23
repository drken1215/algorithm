PC := $(shell hostname -s)
-include Makefile.$(PC)

OJ_VERIFY ?= oj-verify
VERIFY_CXX ?= g++
VERIFY_JOBS ?= 2

.PHONY: configure readme docs verify add add_atcoder add_yukicoder add_aoj add_codeforces add_library_checker add_mytest add_test

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

add_mytest: TEST_DIR = tests/MyTests
add_mytest: TEST_EXT = .test.cpp
add_mytest: TEST_TYPE = mytest
add_mytest: add_test

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
		echo "Usage: make add_atcoder TITLE=EDPC_Z_monge_shortest_path_by_cht"; \
		exit 1; \
	fi
	@if [ -z "$(TEST_DIR)" ] || [ -z "$(TEST_EXT)" ]; then \
		echo "Error: TEST_DIR and TEST_EXT are required."; \
		exit 1; \
	fi
	@mkdir -p "$(TEST_DIR)"
	@cpp="$(TEST_DIR)/$(TITLE)$(TEST_EXT)"; \
	md="$(TEST_DIR)/$(TITLE)$(patsubst %.cpp,%.md,$(TEST_EXT))"; \
	if [ -e "$$cpp" ] || [ -e "$$md" ]; then \
		echo "Error: File already exists."; \
		exit 1; \
	fi; \
	if [ "$(TEST_TYPE)" = "mytest" ]; then \
		printf '%s\n' \
			'#define PROBLEM "https://judge.yosupo.jp/problem/aplusb"' \
			'' \
			'#include <bits/stdc++.h>' \
			'using namespace std;' \
			'' \
			'#include "template.hpp"' \
			'' \
			'void my_test() {' \
			'    ' \
			'}' \
			'' \
			'void librarychecker() {' \
			'    int A, B;' \
			'    cin >> A >> B;' \
			'    cout << A + B << endl;' \
			'}' \
			'' \
			'int main() {' \
			'    my_test();' \
			'    librarychecker();' \
			'}' \
			> "$$cpp"; \
	else \
		printf '%s\n' \
			'#define PROBLEM ""' \
			'' \
			'#include <bits/stdc++.h>' \
			'using namespace std;' \
			'' \
			'int main() {' \
			'    ' \
			'}' \
			> "$$cpp"; \
	fi; \
	if [ "$(TEST_TYPE)" = "mytest" ]; then \
		printf '%s\n' \
			'---' \
			'title: ' \
			'documentation_of: $(TITLE)$(TEST_EXT)' \
			'---' \
			> "$$md"; \
	else \
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
			'## 解法メモ' \
			'' \
			> "$$md"; \
	fi; \
	echo "Created: $$cpp"; \
	echo "Created: $$md"
