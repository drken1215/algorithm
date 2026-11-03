#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

from setup_includes import ROOT, find_include_dirs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", type=int, default=2)
    parser.add_argument("--cxx", default="g++")
    args = parser.parse_args()

    os.chdir(ROOT)

    # oj-helper を import する前に設定する
    os.environ["CXX"] = args.cxx

    # 重要: ここで初めて oj-helper を import
    from onlinejudge_verify.languages.cplusplus import CPlusPlusLanguage
    from onlinejudge_verify.documentation.main import main as generate_docs

    include_dirs = find_include_dirs()

    print("Compiler:", args.cxx, flush=True)
    print("Include directories:", flush=True)
    for directory in include_dirs:
        print(" ", directory, flush=True)

    original_bundle = CPlusPlusLanguage.bundle

    def bundle_with_includes(self, path, *, basedir=Path.cwd(), options):
        options = dict(options)

        original_paths = options.get("include_paths", [])

        options["include_paths"] = list(
            dict.fromkeys(
                [Path(p) for p in original_paths] + include_dirs
            )
        )

        return original_bundle(
            self,
            path,
            basedir=basedir,
            options=options,
        )

    CPlusPlusLanguage.bundle = bundle_with_includes

    try:
        generate_docs(jobs=args.jobs)
    finally:
        CPlusPlusLanguage.bundle = original_bundle


if __name__ == "__main__":
    main()
