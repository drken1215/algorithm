
#!/usr/bin/env python3

import argparse
import json
import os
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

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


def write_vscode_config(dirs, cxx):
    """VS Code IntelliSense 用の設定を生成する。"""
    import platform
    import shutil

    # コンパイラの絶対パスを取得
    compiler = shutil.which(cxx)
    if compiler is None:
        raise RuntimeError(f"Compiler not found: {cxx}")

    # include パスを生成
    include_paths = []

    for directory in dirs:
        relative = directory.relative_to(ROOT)

        if relative == Path("."):
            include_paths.append("${workspaceFolder}")
        else:
            include_paths.append(
                "${workspaceFolder}/" + relative.as_posix()
            )

    # OS と CPU アーキテクチャを判定
    system = platform.system()
    machine = platform.machine()

    if system == "Darwin":
        os_name = "macos"
    elif system == "Linux":
        os_name = "linux"
    else:
        raise RuntimeError(f"Unsupported OS: {system}")

    if machine in ("arm64", "aarch64"):
        arch = "arm64"
    elif machine in ("x86_64", "AMD64"):
        arch = "x64"
    else:
        raise RuntimeError(f"Unsupported architecture: {machine}")

    # VS Code 設定ファイル
    path = ROOT / ".vscode" / "c_cpp_properties.json"
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.exists():
        config = json.loads(path.read_text(encoding="utf-8"))
    else:
        config = {
            "configurations": [],
            "version": 4,
        }

    # 管理対象の設定
    name = "Mac" if system == "Darwin" else "Linux"
    settings = {
        "name": name,
        "includePath": include_paths,
        "compilerPath": compiler,
        "cStandard": "c17",
        "cppStandard": "c++20",
        "intelliSenseMode": f"{os_name}-gcc-{arch}",
    }

    # 既存の設定を更新する
    configurations = config.setdefault("configurations", [])

    for entry in configurations:
        if entry.get("name") == name:
            entry.update(settings)
            break
    else:
        configurations.append(settings)

    config["version"] = 4

    path.write_text(
        json.dumps(config, indent=4, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"Generated: {path.relative_to(ROOT)}")
    print(f"Compiler: {compiler}")
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
    parser.add_argument(
        "--vscode",
        action="store_true",
        help="Generate VS Code IntelliSense configuration",
    )
    args = parser.parse_args()

    dirs = find_include_dirs()

    if args.bundle_flags:
        print_bundle_flags(dirs)
    elif args.vscode:
        write_vscode_config(dirs, args.cxx)
    else:
        write_config(dirs, args.cxx)


if __name__ == "__main__":
    main()
