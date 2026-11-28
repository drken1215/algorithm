from pathlib import Path
import importlib.util


def patch(module, old, new):
    spec = importlib.util.find_spec(module)

    if spec is None or spec.origin is None:
        raise RuntimeError(f"Module not found: {module}")

    path = Path(spec.origin)
    source = path.read_text(encoding="utf-8")

    if new in source:
        print("Already patched:", path)

    elif source.count(old) == 1:
        path.write_text(
            source.replace(old, new, 1),
            encoding="utf-8"
        )
        print("Successfully patched:", path)

    else:
        raise RuntimeError(
            "Patch target not found or ambiguous: " + str(path)
        )


# Suppress pragma once warnings
patch(
    "onlinejudge_verify.languages.cplusplus_bundle",
    "'-fpreprocessed', '-dD', '-E', '-'",
    "'-fpreprocessed', '-dD', '-Wno-pragma-once-outside-header', '-E', '-'"
)


# Avoid unnecessary stack size changes
patch(
    "onlinejudge_verify.verify",
    "_, hard = resource.getrlimit(resource.RLIMIT_STACK)  # type: ignore\n"
    "        resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))",
    "soft, hard = resource.getrlimit(resource.RLIMIT_STACK)  # type: ignore\n"
    "        if soft != hard:\n"
    "            resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))"
)
