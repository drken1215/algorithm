import logging

logging.basicConfig(level=logging.WARNING)

from onlinejudge_bundle.main import main

if __name__ == "__main__":
    raise SystemExit(main())
