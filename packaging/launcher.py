"""Portable distribution launcher for pka_replica (PyInstaller entry).

Layout expectations (relative to this exe's directory):
  _internal/          PyInstaller bundle (rdkit/numpy/openpyxl/requests + pka_replica)
  bin/xtb.exe         bundled GFN2 engine (+ its DLLs)
  _ai/                OPTIONAL: torch/warp/ase/aimnet + models  -> AIMNet2 enabled
  _ai/models/         aimnet2_wb97m_d3_0.pt (preseeded, offline capable)

Launcher duties, before importing pka_replica:
  1. put bundled xtb on PATH (energy.py discovers engines via shutil.which)
  2. if _ai exists: insert it into sys.path and register its DLL dirs,
     point AIMNET_CACHE_DIR at the preseeded models -> full protocol
  3. console: keep window open on fatal errors so the user can read them
"""
import os
import sys


def _exe_dir() -> str:
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


def bootstrap() -> None:
    root = _exe_dir()

    # 1) bundled xtb first on PATH
    xtb_bin = os.path.join(root, "bin")
    if os.path.isdir(xtb_bin):
        os.environ["PATH"] = xtb_bin + os.pathsep + os.environ.get("PATH", "")

    # 2) optional AI addon -> full multistage protocol
    ai = os.path.join(root, "_ai")
    if os.path.isdir(ai):
        sys.path.insert(0, ai)
        # torch/warp ship their DLLs inside the package and self-register via
        # os.add_dll_directory in their __init__; some extension wheels also
        # need the package root registered for sibling .pyd lookups
        for sub in ("", os.path.join("torch", "lib"), os.path.join("warp")):
            d = os.path.join(ai, sub) if sub else ai
            if os.path.isdir(d):
                try:
                    os.add_dll_directory(d)
                except OSError:
                    pass
        models = os.path.join(ai, "models")
        if os.path.isdir(models):
            os.environ["AIMNET_CACHE_DIR"] = models


bootstrap()


def run_server() -> None:
    try:
        from pka_replica.server import main

        port = int(os.environ.get("PKA_PORT", "8765"))
        main(port=port)
    except Exception as e:  # keep the console window readable on fatal errors
        print(f"\n[启动失败] {type(e).__name__}: {e}", file=sys.stderr)
        if os.environ.get("PKA_DEBUG"):
            raise
        try:
            input("\n按回车键关闭窗口…")
        except EOFError:
            pass
        sys.exit(1)


def run_cli() -> None:
    from pka_replica.cli import main

    main()


def dispatch() -> None:
    """One exe, two names: pka-cli.exe runs the CLI, everything else the UI."""
    if "cli" in os.path.basename(sys.argv[0]).lower():
        run_cli()
    else:
        run_server()
