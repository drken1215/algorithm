
#!/usr/bin/env bash
set -euo pipefail

# Move to the repository root
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Environment settings
PC="$(hostname -s)"
VENV="$HOME/.venvs/oj-helper-py312"
TOOLS=".verify-helper/tools"
REQUIREMENTS=".verify-helper/requirements.txt"

echo "=== Setup environment: $PC ==="

# Check macOS
if [[ "$(uname -s)" != "Darwin" ]]; then
    echo "Error: This setup script currently supports macOS only."
    exit 1
fi

# Check Homebrew
if ! command -v brew >/dev/null 2>&1; then
    echo "Error: Homebrew is not installed."
    echo "Please visit https://brew.sh/"
    exit 1
fi

BREW_PREFIX="$(brew --prefix)"
CXX="$BREW_PREFIX/bin/g++-16"
PYTHON="$BREW_PREFIX/bin/python3.12"

# Install GCC
echo "=== GCC ==="

if [[ ! -x "$CXX" ]]; then
    brew install gcc
fi

if [[ ! -x "$CXX" ]]; then
    echo "Error: g++-16 was not found."
    exit 1
fi

"$CXX" --version

# Install Python
echo "=== Python ==="

if [[ ! -x "$PYTHON" ]]; then
    brew install python@3.12
fi

"$PYTHON" --version

# Create Python virtual environment
echo "=== Virtual environment ==="

if [[ ! -e "$VENV" ]]; then
    mkdir -p "$HOME/.venvs"
    "$PYTHON" -m venv "$VENV"
fi

"$VENV/bin/python3" -c \
    'import sys; assert sys.version_info[:2] == (3, 12)'

# Install packages
echo "=== Install oj-helper ==="

"$VENV/bin/python3" -m pip install \
    -r "$REQUIREMENTS"

# Apply patches
echo "=== Apply patches ==="

"$VENV/bin/python3" \
    "$TOOLS/patch_oj_helper.py"

# Generate machine-specific makefile
echo "=== Generate makefile ==="

MAKEFILE="makefile.$PC"

if [[ -e "$MAKEFILE" ]]; then
    echo "Already exists: $MAKEFILE (skipped)"
else
    cat > "$MAKEFILE" <<EOF
OJ_HELPER := \$(HOME)/.venvs/oj-helper-py312

OJ_BUNDLE := \$(OJ_HELPER)/bin/python3 .verify-helper/tools/bundle_quiet.py
BUNDLE_CXX := $CXX

OJ_VERIFY := \$(OJ_HELPER)/bin/oj-verify
VERIFY_CXX := $CXX

export CXX := \$(VERIFY_CXX)
export PATH := \$(OJ_HELPER)/bin:\$(PATH)
EOF

    echo "Generated: $MAKEFILE"
fi

# Generate VS Code settings
echo "=== VS Code settings ==="

"$VENV/bin/python3" \
    "$TOOLS/setup_includes.py" \
    --vscode --cxx "$CXX"

# Final checks
echo "=== Environment check ==="

"$VENV/bin/python3" --version
"$VENV/bin/python3" -c \
    'import ssl; print(ssl.OPENSSL_VERSION)'
"$VENV/bin/oj" --version

echo
echo "=== Setup completed ==="
echo
echo 'Add the following lines to ~/.zshrc if needed:'
echo
echo 'typeset -U path PATH'
echo 'export PATH="$HOME/.venvs/oj-helper-py312/bin:$PATH"'
echo 'ulimit -Ss "$(ulimit -Hs)"'
echo
echo 'Set YUKICODER_TOKEN separately.'
echo
echo 'Next: make, make verify, make docs'
