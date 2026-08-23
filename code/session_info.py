"""Record the computational environment (Python analogue of sessionInfo()).

Usage: python session_info.py   -> writes ../results/session_info.txt
"""
import platform, sys, importlib, datetime

MODULES = ["numpy", "scipy", "pandas", "sklearn", "matplotlib"]

def main():
    lines = [
        f"date: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"python: {sys.version.replace(chr(10), ' ')}",
        f"platform: {platform.platform()}",
        f"machine: {platform.machine()}",
    ]
    for m in MODULES:
        try:
            mod = importlib.import_module(m)
            lines.append(f"{m}: {mod.__version__}")
        except Exception as e:  # pragma: no cover
            lines.append(f"{m}: NOT AVAILABLE ({e})")
    out = "\n".join(lines) + "\n"
    with open("../results/session_info.txt", "w") as f:
        f.write(out)
    print(out)

if __name__ == "__main__":
    main()
