#!/usr/bin/env python3

import os
import platform
import sys
import time
import urllib.request
from pathlib import Path


VIRTUALENV_VERSION = "20.27.1"
URL = (
    "https://raw.githubusercontent.com/pypa/get-virtualenv/"
    f"{VIRTUALENV_VERSION}/public/virtualenv.pyz"
)


def cache_dir() -> Path:
    system = platform.system()
    if system == "Windows":
        localappdata = os.environ.get("LOCALAPPDATA")
        if not localappdata:
            raise RuntimeError("LOCALAPPDATA is not set")
        return Path(localappdata) / "pypa" / "cibuildwheel" / "Cache"
    if system == "Darwin":
        return Path.home() / "Library" / "Caches" / "cibuildwheel"
    return Path.home() / ".cache" / "cibuildwheel"


def download_with_retry(url: str, dest: Path, attempts: int = 5) -> None:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                dest.write_bytes(resp.read())
            return
        except Exception as exc:
            last_error = exc
            if attempt == attempts:
                break
            print(
                f"Download failed (attempt {attempt}/{attempts}): {exc}. Retrying...",
                file=sys.stderr,
            )
            time.sleep(2 * attempt)
    raise RuntimeError(f"Failed to download {url}") from last_error


def main() -> int:
    cdir = cache_dir()
    cdir.mkdir(parents=True, exist_ok=True)
    dest = cdir / f"virtualenv-{VIRTUALENV_VERSION}.pyz"

    if dest.exists() and dest.stat().st_size > 0:
        print(f"Using existing cibuildwheel virtualenv cache: {dest}")
        return 0

    print(f"Prefetching cibuildwheel virtualenv helper to {dest}")
    download_with_retry(URL, dest)
    print(f"Saved {dest} ({dest.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
