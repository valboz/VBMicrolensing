#!/usr/bin/env python3

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path


def fetch_sdist(package: str, version: str) -> tuple[str, str]:
    url = f"https://pypi.org/pypi/{package}/{version}/json"
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)

    sdists = [u for u in data.get("urls", []) if u.get("packagetype") == "sdist"]
    if not sdists:
        raise SystemExit(f"No sdist found on PyPI for version {version}")

    sdist = sdists[0]
    return sdist["url"], sdist["digests"]["sha256"]


def patch_recipe(
    recipe_path: Path, version: str, sha256: str, *, drop_stdlib_c: bool = False
) -> None:
    text = recipe_path.read_text()

    text, n1 = re.subn(
        r'(\{%\s*set\s+version\s*=\s*")[^"]+("\s*%\})',
        rf"\g<1>{version}\2",
        text,
        count=1,
    )
    text, n2 = re.subn(
        r"(^\s*sha256:\s*)[0-9a-f]{64}(\s*$)",
        rf"\g<1>{sha256}\2",
        text,
        count=1,
        flags=re.MULTILINE,
    )

    if n1 != 1 or n2 != 1:
        raise SystemExit("Failed to patch conda/meta.yaml version/sha256")

    if drop_stdlib_c:
        text, n3 = re.subn(
            r"(?m)^\s*-\s+\{\{\s*stdlib\('c'\)\s*\}\}\s*$\n?",
            "",
            text,
            count=1,
        )
        if n3:
            sys.stderr.write(
                "Smoke-test patch: removed {{ stdlib('c') }} requirement for standalone CI\n"
            )

    recipe_path.write_text(text)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch PyPI sdist metadata and patch a conda recipe version/sha256."
    )
    parser.add_argument("version", help="Version to fetch from PyPI")
    parser.add_argument(
        "--package", default="vbmicrolensing", help="PyPI package name (default: vbmicrolensing)"
    )
    parser.add_argument(
        "--recipe", default="conda/meta.yaml", help="Path to recipe file (default: conda/meta.yaml)"
    )
    parser.add_argument(
        "--drop-stdlib-c",
        action="store_true",
        help="Remove {{ stdlib('c') }} from build requirements (for standalone smoke CI only)",
    )
    args = parser.parse_args()

    source_url, source_sha256 = fetch_sdist(args.package, args.version)
    patch_recipe(
        Path(args.recipe),
        args.version,
        source_sha256,
        drop_stdlib_c=args.drop_stdlib_c,
    )

    sys.stdout.write(source_url + "\n")
    sys.stdout.write(source_sha256 + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
