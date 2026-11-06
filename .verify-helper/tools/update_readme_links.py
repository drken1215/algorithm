
#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"

PAGES_BASE = "https://drken1215.github.io/algorithm/"
GITHUB_PREFIX = "https://github.com/drken1215/algorithm/blob/"

EXCLUDED = {
    ".git",
    ".github",
    ".verify-helper",
    "__pycache__",
}


def find_sources():
    """利用可能な .hpp と tests/*.cpp を索引化する。"""
    index = defaultdict(list)

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue

        relative = path.relative_to(ROOT)

        if any(part in EXCLUDED or part.startswith(".")
               for part in relative.parts):
            continue

        is_header = path.suffix == ".hpp"
        is_test = (
            path.suffix == ".cpp"
            and len(relative.parts) >= 2
            and relative.parts[0] == "tests"
        )

        if is_header or is_test:
            index[path.name].append(relative.as_posix())

    return index


def parse_source_url(url):
    """
    URL のドメインやディレクトリ構成に関係なく、
    最後のファイル名を使って検索する。
    """
    parts = urlsplit(url)

    path = unquote(parts.path)

    # oj-helper の .html ページにも対応
    if path.endswith(".hpp.html") or path.endswith(".cpp.html"):
        path = path[:-5]

    # .hpp / .cpp 以外のリンクは対象外
    if not path.endswith((".hpp", ".cpp")):
        return None

    return path


def choose_path(old_path, index):
    """完全一致を優先し、なければファイル名で探索する。"""
    old_path = old_path.lstrip("/")
    basename = Path(old_path).name

    if not basename.endswith((".hpp", ".cpp")):
        return None, None

    candidates = index.get(basename, [])

    # 元のパスが現在も存在する場合
    if old_path in candidates:
        return old_path, None

    # ファイル名が一意なら、そのファイルを使用
    if len(candidates) == 1:
        return candidates[0], None

    # 同名ファイルが複数ある場合は変更しない
    if len(candidates) > 1:
        return None, f"Ambiguous: {old_path} -> {candidates}"

    return None, None


def rewrite_url(url, index):
    old_path = parse_source_url(url)

    if old_path is None:
        return url, None

    new_path, warning = choose_path(old_path, index)

    if new_path is None:
        return url, warning

    parts = urlsplit(url)

    # .hpp は常に oj-helper の HTML ページへリンクする
    if new_path.endswith(".hpp"):
        html_suffix = ".html"
    else:
        # .cpp は従来どおり元の URL 形式を維持する
        html_suffix = ".html" if parts.path.endswith(".html") else ""

    new_url = PAGES_BASE + quote(new_path, safe="/()_.-")
    new_url += html_suffix

    if parts.query:
        new_url += "?" + parts.query

    if parts.fragment:
        new_url += "#" + parts.fragment

    return new_url, None


def update_links(content, index):
    """
    Markdown の [表示名](URL) を処理する。
    URL 内の括弧にも対応する。
    """
    result = []
    changes = []
    warnings = []

    cursor = 0

    for match in re.finditer(r"\]\(", content):
        start = match.end()

        if start < cursor:
            continue

        depth = 1
        pos = start

        while pos < len(content) and depth:
            if content[pos] == "\\":
                pos += 2
                continue

            if content[pos] == "(":
                depth += 1
            elif content[pos] == ")":
                depth -= 1

            pos += 1

        if depth:
            continue

        # URL のみを取得。タイトル付きリンクは対象外
        url = content[start:pos - 1]

        if any(c.isspace() for c in url):
            continue

        new_url, warning = rewrite_url(url, index)

        if warning:
            warnings.append(warning)

        if new_url != url:
            result.append(content[cursor:start])
            result.append(new_url)
            cursor = pos - 1
            changes.append((url, new_url))

    result.append(content[cursor:])

    return "".join(result), changes, warnings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="Preview changes without modifying README.md",
    )
    args = parser.parse_args()

    index = find_sources()
    original = README.read_text(encoding="utf-8")

    updated, changes, warnings = update_links(original, index)

    for old, new in changes:
        print(f"UPDATE: {old}")
        print(f"     -> {new}")

    for warning in sorted(set(warnings)):
        print(f"WARNING: {warning}")

    if changes and not args.check:
        README.write_text(updated, encoding="utf-8")

    print(f"Updated links: {len(changes)}")

    if args.check:
        print("Check only: README.md was not modified.")


if __name__ == "__main__":
    main()
