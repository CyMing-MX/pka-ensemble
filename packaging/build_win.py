"""Build the portable Windows distribution for pka_replica.

Usage (from the repo root, with the dev Python that has all deps):
    python packaging/build_win.py            # skeleton + _ai + zips

One PyInstaller exe ("pka") is copied to pka-server.exe / pka-cli.exe — the
launcher picks the mode from the exe filename, so both share one _internal.

Products (in packaging/dist):
    pka_replica_v<ver>_lite.zip      skeleton (sites + xtb)
    pka_replica_v<ver>_full.zip      skeleton + _ai (AIMNet2 multistage)
    pka_replica_v<ver>_ai-addon.zip  the _ai folder alone (upgrade path)
    SHA256SUMS.txt
"""
import hashlib
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PY = sys.executable
SP = Path(PY).parent / "Lib" / "site-packages"
VERSION = "1.0"
XTB_SRC = Path(r"D:\GUASSIAN\xtb-6.7.1\bin")

AI_PACKAGES = [
    "torch", "warp", "ase", "aimnet",
    "sympy", "networkx", "jinja2", "filelock", "fsspec",
    "typing_extensions", "MarkupSafe", "mpmath", "torchgen",
    "functorch",
]
AI_MODEL = Path.home() / ".cache" / "aimnet" / "aimnet2_wb97m_d3_0.pt"
EXCLUDE_DIRS = {"__pycache__", "tests", "test"}


def sh(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), check=True)


def write_entries() -> None:
    # bootstrap must run before pka_replica imports: PATH/_ai/AIMNET_CACHE_DIR
    (HERE / "entry_server.py").write_text(
        "import launcher\nlauncher.run_server()\n", encoding="utf-8")
    (HERE / "entry_cli.py").write_text(
        "import launcher\nlauncher.run_cli()\n", encoding="utf-8")


def build_exe(out_dir: Path) -> None:
    work = HERE / "build_work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    sh([
        PY, "-m", "PyInstaller",
        "--noconfirm", "--clean",
        "--name", "pka",
        "--distpath", str(out_dir),
        "--workpath", str(work),
        "--specpath", str(work),
        "--console",
        "--collect-all", "rdkit",
        "--paths", str(HERE),           # launcher.py
        "--paths", str(REPO),           # pka_replica package
        "--add-data", f"{REPO / 'pka_replica' / 'web'};web",
        "--exclude-module", "torch",
        "--exclude-module", "warp",
        "--exclude-module", "aimnet",
        "--exclude-module", "ase",
        "--exclude-module", "tkinter",
        # pulled in by rdkit's Draw hooks; we only use rdMolDraw2D (SVG)
        "--exclude-module", "matplotlib",
        "--exclude-module", "PyQt5",
        "--exclude-module", "PySide2",
        "--exclude-module", "PySide6",
        "--exclude-module", "PIL",
        "--exclude-module", "IPython",
        str(HERE / "entry_pka.py"),
    ])


def write_licenses(root: Path) -> None:
    lic = root / "LICENSES"
    lic.mkdir(exist_ok=True)
    (lic / "README.txt").write_text(
        "本软件分发的第三方组件及其许可:\n\n"
        "- xtb (tight-binding calculator) — LGPL-3.0,\n"
        "  https://github.com/grimme-lab/xtb 。依 LGPL 条款以原样二进制分发,\n"
        "  源码可在上游仓库获取。\n"
        "- RDKit — BSD-3-Clause, https://www.rdkit.org/\n"
        "- NumPy — BSD-3-Clause\n"
        "- openpyxl — MIT\n"
        "- requests — Apache-2.0\n"
        "- PyTorch — BSD-3-Clause\n"
        "- warp — NVIDIA Source Code License(允许再分发)\n"
        "- aimnet / AIMNet2 模型 — https://github.com/isayevlab/aimnetcentral\n"
        "- 本项目代码 — 见仓库 LICENSE。\n",
        encoding="utf-8")


def assemble_skeleton(dist_exe: Path, out_root: Path) -> Path:
    root = out_root / f"pka_replica_v{VERSION}_lite"
    if root.exists():
        shutil.rmtree(root)
    shutil.copytree(dist_exe / "pka", root)
    # one bundle, two names: the launcher reads the exe basename
    shutil.copy2(root / "pka.exe", root / "pka-server.exe")
    shutil.copy2(root / "pka.exe", root / "pka-cli.exe")
    bin_dst = root / "bin"
    bin_dst.mkdir(exist_ok=True)
    for f in XTB_SRC.iterdir():
        if f.is_file():
            shutil.copy2(f, bin_dst / f.name)
    write_licenses(root)
    manual = HERE / "USER_MANUAL.md"
    if manual.exists():
        shutil.copy2(manual, root / "操作手册.md")
    return root


def collect_ai(out_root: Path) -> Path:
    ai = out_root / "_ai"
    if ai.exists():
        shutil.rmtree(ai)
    ai.mkdir(parents=True)
    want = {p.lower() for p in AI_PACKAGES}
    for pkg in AI_PACKAGES:
        src = SP / pkg
        if not src.is_dir():
            print(f"  [skip] {pkg} (not installed)")
            continue
        print(f"  [copy] {pkg}", flush=True)
        shutil.copytree(src, ai / pkg,
                        ignore=shutil.ignore_patterns(*EXCLUDE_DIRS, "*.pyc"))
    for d in SP.glob("*.dist-info"):
        if d.name.split("-")[0].lower().replace("_", "-") in want:
            shutil.copytree(d, ai / d.name)
    models = ai / "models"
    models.mkdir(exist_ok=True)
    if AI_MODEL.exists():
        shutil.copy2(AI_MODEL, models / AI_MODEL.name)
    else:
        print("  [warn] AIMNet2 model not found; first run will download")
    return ai


def zip_dir(src: Path, dst: Path, arc_root: str | None = None) -> None:
    print(f"  [zip] {dst.name}", flush=True)
    with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        base = src if arc_root is None else src.parent
        for p in src.rglob("*"):
            rel = p.relative_to(base)
            zf.write(p, str(rel))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    out = HERE / "dist"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    write_entries()
    print("== [1/5] PyInstaller build", flush=True)
    build_exe(out)

    print("== [2/5] assemble skeleton", flush=True)
    root = assemble_skeleton(out, out)

    print("== [3/5] collect _ai addon", flush=True)
    ai = collect_ai(out)

    print("== [4/5] full = skeleton + _ai", flush=True)
    full = out / f"pka_replica_v{VERSION}_full"
    shutil.copytree(root, full)
    shutil.copytree(ai, full / "_ai")

    print("== [5/5] zips + checksums", flush=True)
    lite_zip = out / f"pka_replica_v{VERSION}_lite.zip"
    full_zip = out / f"pka_replica_v{VERSION}_full.zip"
    ai_zip = out / f"pka_replica_v{VERSION}_ai-addon.zip"
    zip_dir(root, lite_zip)
    zip_dir(full, full_zip)
    zip_dir(ai, ai_zip, arc_root="_ai")
    sums = out / "SHA256SUMS.txt"
    with open(sums, "w", encoding="utf-8") as fh:
        for z in (lite_zip, full_zip, ai_zip):
            fh.write(f"{sha256(z)}  {z.name}\n")
    print("\nDONE. Products in", out, flush=True)
    for f in sorted(out.glob("*.zip")):
        print(f"  {f.name:40s} {f.stat().st_size/1e6:8.1f} MB", flush=True)


if __name__ == "__main__":
    main()
