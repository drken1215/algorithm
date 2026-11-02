
#!/usr/bin/env python3

import argparse
import json
import os
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ライブラリ探索の対象外とするディレクトリ
EXCLUDED = {
    ".git",
    ".verify-helper",
    ".github",
    "test",
    "tests",
    "docs",
    "tools",
    "include",
    "__pycache__",
}


def find_include_dirs():
    """リポジトリ内の .hpp を持つディレクトリを列挙する。"""
    dirs = {ROOT}
    filenames = {}

    for current, children, files in os.walk(ROOT):
        children[:] = sorted(
            d for d in children
            if d not in EXCLUDED and not d.startswith(".")
        )

        directory = Path(current)
        for filename in files:
            if not filename.endswith(".hpp"):
                continue

            path = directory / filename
            dirs.add(directory)
            filenames.setdefault(filename, []).append(path)

    # 短い include 名の衝突を検出
    duplicates = {
        name: paths
        for name, paths in filenames.items()
        if len(paths) > 1
    }

    if duplicates:
        print("Duplicate header names detected:", file=sys.stderr)
        for name, paths in sorted(duplicates.items()):
            print("  " + name, file=sys.stderr)
            for path in paths:
                print("    " + str(path.relative_to(ROOT)),
                      file=sys.stderr)
        raise SystemExit(1)

    # リポジトリルートを先頭にし、残りは名前順
    return [ROOT] + sorted(dirs - {ROOT})


def write_config(dirs, cxx):
    """oj-helper 用の config.toml を生成する。"""
    flags = ["-std=c++20", "-O2"]

    for directory in dirs:
        relative = directory.relative_to(ROOT)
        name = "." if relative == Path(".") else relative.as_posix()
        flags.append("-I" + name)

    config = ROOT / ".verify-helper" / "config.toml"
    config.parent.mkdir(parents=True, exist_ok=True)

    contents = (
        "[[languages.cpp.environments]]\n"
        f"CXX = {json.dumps(cxx)}\n"
        f"CXXFLAGS = {json.dumps(flags)}\n"
    )
    config.write_text(contents, encoding="utf-8")

    print(f"Generated: {config.relative_to(ROOT)}")
    print(f"Compiler: {cxx}")
    print(f"Include directories: {len(dirs)}")


def print_bundle_flags(dirs):
    """コンテスト用 Makefile で使う -I オプションを出力する。"""
    print(" ".join(
        "-I " + shlex.quote(str(directory))
        for directory in dirs
    ))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cxx",
        default="g++",
        help="Compiler used by oj-verify",
    )
    parser.add_argument(
        "--bundle-flags",
        action="store_true",
        help="Print include flags for oj-bundle",
    )
    args = parser.parse_args()

    dirs = find_include_dirs()

    if args.bundle_flags:
        print_bundle_flags(dirs)
    else:
        write_config(dirs, args.cxx)


if __name__ == "__main__":
    main()
