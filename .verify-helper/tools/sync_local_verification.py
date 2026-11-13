
#!/usr/bin/env python3

import hashlib
import json
import pathlib
import sys

from onlinejudge_verify import marker
from onlinejudge_verify.languages import list as languages

ROOT = pathlib.Path.cwd().resolve()
LOCAL = ROOT / ".verify-helper/timestamps.local.json"
REMOTE = ROOT / ".verify-helper/timestamps.remote.json"
MANIFEST = ROOT / ".verify-helper/verified-local.json"


def test_files():
    return sorted(pathlib.Path("tests").rglob("*.test.cpp"))


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependencies(test):
    language = languages.get(test)
    if language is None:
        raise RuntimeError(f"Unsupported language: {test}")

    result = {}
    for dep in language.list_dependencies(test, basedir=ROOT):
        absolute = dep.resolve(strict=True)
        try:
            relative = absolute.relative_to(ROOT)
        except ValueError:
            raise RuntimeError(
                f"Dependency outside repository: {absolute}"
            )
        result[relative.as_posix()] = file_hash(absolute)

    # テスト本体も必ず記録する
    result[test.as_posix()] = file_hash(ROOT / test)
    return dict(sorted(result.items()))


def export_local():
    if not LOCAL.exists():
        raise RuntimeError(
            "Local verification history not found. "
            "Run make verify first."
        )

    local_marker = marker.VerificationMarker(
        json_path=LOCAL,
        use_git_timestamp=False,
        jobs=None,
    )

    data = {}
    missing = []

    for test in test_files():
        name = test.as_posix()
        if not local_marker.is_verified(test):
            missing.append(name)
            continue

        data[name] = dependencies(test)

    if missing:
        raise RuntimeError(
            "Local verification is missing or outdated:\n"
            + "\n".join(missing)
            + "\nRun make verify again."
        )

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Exported {len(data)} local verifications.")


def import_remote():
    if not MANIFEST.exists():
        raise RuntimeError(
            "Missing verified-local.json. "
            "Run make verify locally and commit "
            ".verify-helper/verified-local.json."
        )

    saved = json.loads(MANIFEST.read_text(encoding="utf-8"))

    # Actions 側の検証履歴はローカル結果のみから再構築
    remote_data = {}
    missing = []

    for test in test_files():
        name = test.as_posix()
        expected = saved.get(name)

        if expected is None:
            missing.append(name)
            continue

        try:
            actual = dependencies(test)
        except Exception as e:
            print(f"Dependency error in {name}: {e}", file=sys.stderr)
            missing.append(name)
            continue

        if actual != expected:
            missing.append(name)
            continue

        # GitHub Actions 側のコミット日時を使用
        timestamp = marker.get_last_commit_time_to_verify(test)
        remote_data[name] = timestamp.strftime(
            "%Y-%m-%d %H:%M:%S %z"
        )

    if missing:
        raise RuntimeError(
            "Local verification is missing or outdated:\n"
            + "\n".join(missing)
            + "\nRe-run make verify locally and commit "
              ".verify-helper/verified-local.json."
        )

    REMOTE.parent.mkdir(parents=True, exist_ok=True)
    REMOTE.write_text(
        json.dumps(remote_data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    # GitHub Actions と同じ判定方式で全件を確認
    remote_marker = marker.VerificationMarker(
        json_path=REMOTE,
        use_git_timestamp=True,
        jobs=None,
    )

    unverified = [
        test.as_posix()
        for test in test_files()
        if not remote_marker.is_verified(test)
    ]

    if unverified:
        raise RuntimeError(
            "Remote verification status is not verified:\n"
            + "\n".join(unverified)
        )

    print(f"Imported {len(remote_data)} local verifications.")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in (
        "export", "import"
    ):
        print("Usage: sync_local_verification.py export|import")
        sys.exit(2)

    try:
        if sys.argv[1] == "export":
            export_local()
        else:
            import_remote()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
