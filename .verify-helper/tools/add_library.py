
#!/usr/bin/env python3

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def add_library(lib):
    relative = Path(lib)

    # 不正なパスを防ぐ
    if (
        not lib
        or relative.is_absolute()
        or ".." in relative.parts
        or relative.suffix
    ):
        raise SystemExit(
            "Error: LIB must be a relative path without extension."
        )

    name = relative.name
    hpp = ROOT / (str(relative) + ".hpp")
    md = ROOT / (str(relative) + ".md")

    # 既存のファイルは上書きしない
    for path in (hpp, md):
        if path.exists():
            raise SystemExit(f"Error: already exists: {path}")

    # ディレクトリを作成
    hpp.parent.mkdir(parents=True, exist_ok=True)

    hpp_content = (
        "#pragma once\n"
        "#include <bits/stdc++.h>\n"
        "#include <cassert>\n"
        "using namespace std;\n"
    )

    md_content = (
        "---\n"
        f"title: {name}\n"
        f"documentation_of: ./{name}.hpp\n"
        "---\n"
        "\n"
        "## 概要\n"
        "\n"
        "## 仕様\n"
        "\n"
        "## 制約\n"
        "\n"
        "## 問題例\n"
        "\n"
        "## 参考文献\n"
        "\n"
    )

    hpp.write_text(hpp_content, encoding="utf-8")
    md.write_text(md_content, encoding="utf-8")

    print(f"Created: {hpp.relative_to(ROOT)}")
    print(f"Created: {md.relative_to(ROOT)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lib",
        required=True,
        help="Library path, e.g. Others/io_for_contests",
    )
    args = parser.parse_args()
    add_library(args.lib)


if __name__ == "__main__":
    main()
