import sys
import time

sys.path.insert(0, r"D:\zcode_rea\pka_replica")
import warnings

warnings.filterwarnings("ignore")

from pka_replica.workflow import run_micropka, run_macropka  # noqa: E402


def main():
    t0 = time.time()
    m = run_micropka(
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc2cccnc12", backend="xtb")
    print(f"micropka imatinib: {time.time()-t0:.1f}s, sites:",
          len(m["conjugate_bases"]) + len(m["conjugate_acids"]), flush=True)
    r = run_macropka("c1ccncc1", backend="xtb")
    print("macropka pyridine:", r["pKa_values"][0]["pKa"], flush=True)


if __name__ == "__main__":
    main()
