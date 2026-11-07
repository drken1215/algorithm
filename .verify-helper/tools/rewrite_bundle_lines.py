
#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from urllib.parse import quote

PAGES_BASE = "https://drken1215.github.io/algorithm/"

LINE_PATTERN = re.compile(
    r'^(\s*#line\s+\d+\s+)"([^"\n]+)"',
    re.MULTILINE,
)


def rewrite(content, library_root):
    library_root = library_root.resolve()

    def replace(match):
        prefix, filename = match.groups()
        path = Path(filename)

        if not path.is_absolute():
            return match.group(0)

        try:
            relative = path.resolve().relative_to(library_root)
        except ValueError:
            return match.group(0)

        # ライブラリのヘッダファイルだけを変換
        if relative.suffix != ".hpp":
            return match.group(0)

        url = PAGES_BASE + quote(
            relative.as_posix(), safe="/()_.-"
        ) + ".html"

        return f'{prefix}"{url}"'

    return LINE_PATTERN.sub(replace, content)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path)
    parser.add_argument("--library", type=Path, required=True)
    args = parser.parse_args()

    original = args.file.read_text(encoding="utf-8")
    updated = rewrite(original, args.library)

    args.file.write_text(updated, encoding="utf-8")

    count = sum(
        1 for a, b in zip(
            original.splitlines(),
            updated.splitlines(),
        ) if a != b
    )
    print(f"Rewritten #line directives: {count}")


if __name__ == "__main__":
    main()
