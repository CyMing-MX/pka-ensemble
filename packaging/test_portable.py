"""Clean-environment acceptance test for the portable build.

1. extract lite zip to a scratch dir
2. launch pka-server.exe with a MINIMAL PATH (simulates a machine without
   Python/xtb on PATH) on a test port
3. hit the API: acetic acid sites + xtb (proves bundled bin/xtb.exe is found)
4. pka-cli.exe smoke test
5. extract ai-addon over it, restart, run xtb again (proves AIMNet2 chain)
"""
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
DIST = HERE / "dist"
SCRATCH = Path(r"D:\zcode_rea\pka_replica\packaging\test_env")
PORT = 8899

MINIMAL_PATH = os.pathsep.join([
    r"C:\Windows\System32", r"C:\Windows", r"C:\Windows\System32\Wbem",
    r"C:\Windows\System32\WindowsPowerShell\v1.0",
])
CLEAN_ENV = {
    "PATH": MINIMAL_PATH,
    "SYSTEMROOT": os.environ.get("SYSTEMROOT", r"C:\Windows"),
    "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
    "TMP": os.environ.get("TMP", r"C:\Windows\Temp"),
    "USERPROFILE": os.environ.get("USERPROFILE", ""),
    "HOMEDRIVE": os.environ.get("HOMEDRIVE", "C:"),
    "HOMEPATH": os.environ.get("HOMEPATH", ""),
    "COMPUTERNAME": os.environ.get("COMPUTERNAME", "PC"),
    "APPDATA": os.environ.get("APPDATA", ""),
    "LOCALAPPDATA": os.environ.get("LOCALAPPDATA", ""),
    "PROGRAMDATA": os.environ.get("PROGRAMDATA", ""),
}


def api(path, body=None, timeout=600):
    if body is None:
        r = urllib.request.urlopen(f"http://127.0.0.1:{PORT}{path}", timeout=timeout)
    else:
        req = urllib.request.Request(
            f"http://127.0.0.1:{PORT}{path}", data=json.dumps(body).encode(),
            headers={"Content-Type": "application/json"})
        r = urllib.request.urlopen(req, timeout=timeout)
    return json.load(r)


def start_server(env):
    exe = SCRATCH / "pka-server.exe"
    env = dict(env)
    env["PKA_PORT"] = str(PORT)
    return subprocess.Popen(
        [str(exe), "--no-browser"], env=env, cwd=str(SCRATCH),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        creationflags=subprocess.CREATE_NEW_CONSOLE)


def wait_up(proc, seconds=180):
    deadline = time.time() + seconds
    while time.time() < deadline:
        if proc.poll() is not None:
            out = proc.stdout.read().decode("utf-8", "replace") if proc.stdout else ""
            raise RuntimeError(f"server exited early:\n{out[-2000:]}")
        try:
            api("/api/history", timeout=5)
            return
        except Exception:
            time.sleep(2)
    raise RuntimeError("server did not come up")


def stop(proc):
    proc.terminate()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        proc.kill()


def main():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    SCRATCH.mkdir(parents=True)

    print("== [1] extract lite", flush=True)
    with zipfile.ZipFile(DIST / "pka_replica_v1.0_lite.zip") as zf:
        zf.extractall(SCRATCH)
    print("   extracted:", len(list(SCRATCH.iterdir())), "entries", flush=True)

    env = dict(CLEAN_ENV)
    print("== [2] start server (minimal PATH)", flush=True)
    proc = start_server(env)
    try:
        wait_up(proc)
        print("   server up", flush=True)

        print("== [3] sites route", flush=True)
        m = api("/api/predict", {"smiles": "CC(=O)O", "backend": "sites"})
        print("   backend =", m["macro"]["object_data"]["energy_backend"],
              " pKa =", [x["pKa"] for x in m["macro"]["pKa_values"]], flush=True)
        assert m["macro"]["object_data"]["energy_backend"] == "fragment"

        print("== [4] xtb route (bundled bin/xtb.exe)", flush=True)
        t0 = time.time()
        m = api("/api/predict", {"smiles": "CC(=O)O", "backend": "xtb"})
        bk = m["macro"]["object_data"]["energy_backend"]
        print(f"   backend = {bk}  pKa =",
              [x["pKa"] for x in m["macro"]["pKa_values"]],
              f" ({time.time()-t0:.0f}s)", flush=True)
        if "aimnet" in bk or "xtb" not in bk:
            raise RuntimeError(f"unexpected backend: {bk}")
        print("   AIMNet2 in bundle: expected NO (lite) ->", "aimnet" in bk, flush=True)
    finally:
        stop(proc)

    print("== [5] pka-cli smoke (minimal PATH)", flush=True)
    r = subprocess.run(
        [str(SCRATCH / "pka-cli.exe"), "CC(=O)O", "--backend", "sites"],
        env=env, cwd=str(SCRATCH), capture_output=True, text=True,
        timeout=300, encoding="utf-8", errors="replace")
    print("   exit", r.returncode, "|", (r.stdout or r.stderr).strip().splitlines()[-1][:80], flush=True)
    assert r.returncode == 0

    print("== [6] extract ai-addon over lite -> full", flush=True)
    with zipfile.ZipFile(DIST / "pka_replica_v1.0_ai-addon.zip") as zf:
        zf.extractall(SCRATCH)
    assert (SCRATCH / "_ai" / "models" / "aimnet2_wb97m_d3_0.pt").exists()

    print("== [7] start server again (AI enabled)", flush=True)
    proc = start_server(env)
    try:
        wait_up(proc)
        t0 = time.time()
        m = api("/api/predict", {"smiles": "CC(=O)O", "backend": "xtb"})
        print("   xtb pKa =", [x["pKa"] for x in m["macro"]["pKa_values"]],
              f" ({time.time()-t0:.0f}s; first call compiles warp kernels)",
              flush=True)
        # sanity: full-protocol acetic acid should be close to exp 4.76
        v = m["micro"]["strongest_acid"]
        print("   strongest acid =", v, "(exp 4.76)", flush=True)
        assert abs(v - 4.76) < 1.0
    finally:
        stop(proc)

    print("\nALL ACCEPTANCE TESTS PASSED", flush=True)


if __name__ == "__main__":
    main()
