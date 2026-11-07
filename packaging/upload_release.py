"""Upload release assets with per-file retry + state verification.

Network drops mid-upload are routine on long transfers to GitHub; gh exits
non-zero and the asset is left absent/partial, so each file is retried until
the release shows it 'uploaded'.
"""
import subprocess
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")

FILES = [
    r"D:\zcode_rea\pka_replica\packaging\dist\pka_replica_v1.0_lite.zip",
    r"D:\zcode_rea\pka_replica\packaging\dist\pka_replica_v1.0_full.zip",
    r"D:\zcode_rea\pka_replica\packaging\dist\pka_replica_v1.0_ai-addon.zip",
    r"D:\zcode_rea\pka_replica\packaging\dist\SHA256SUMS.txt",
]


def assets():
    out = subprocess.run(
        ["gh", "api", "repos/CyMing-MX/pka-ensemble/releases",
         "--jq", ".[0].assets[] | \"\\(.name)|\\(.state)\""],
        capture_output=True, text=True, timeout=60)
    return dict(line.split("|") for line in out.stdout.splitlines() if "|")


def uploaded(name):
    return assets().get(name) == "uploaded"


def main():
    for path in FILES:
        name = path.rsplit("\\", 1)[-1]
        for attempt in range(1, 5):
            if uploaded(name):
                print(f"[skip] {name} already uploaded", flush=True)
                break
            print(f"[upload] {name} attempt {attempt}", flush=True)
            r = subprocess.run(
                ["gh", "release", "upload", "v1.0", path, "--clobber",
                 "--repo", "CyMing-MX/pka-ensemble"],
                capture_output=True, text=True, timeout=3600)
            if r.returncode == 0 and uploaded(name):
                print(f"[ok] {name}", flush=True)
                break
            print(f"  failed: {(r.stderr or r.stdout).strip()[-200:]}",
                  flush=True)
            time.sleep(10 * attempt)
        else:
            print(f"[FAIL] {name} after 4 attempts", flush=True)
            sys.exit(1)
    print("\nFINAL STATE:", flush=True)
    for k, v in assets().items():
        print(f"  {k}  {v}", flush=True)


if __name__ == "__main__":
    main()
