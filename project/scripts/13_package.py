"""
Package the data folder for submission.

Zips the whole project tree (structures, ligands, Vina configs and logs,
poses, complexes, scripts, analysis, figures and the report) into
C:/CADD/GSK3B_CADD_miniproject_data.zip, skipping caches and the raw
P2Rank visualisation payload.

Run: python 13_package.py
"""
import os, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(os.path.dirname(ROOT), "GSK3B_CADD_miniproject_data.zip")

SKIP_DIRS = {"__pycache__", ".ipynb_checkpoints", "_embed"}
SKIP_SUFFIX = (".pyc", ".tmp")
# a stub left locked open by a Word automation instance; not part of the work
SKIP_FILES = {"GSK3B_virtual_screen_report.docx"}


def main():
    n, total = 0, 0
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for dirpath, dirnames, filenames in os.walk(ROOT):
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
            for f in filenames:
                if f.endswith(SKIP_SUFFIX) or f in SKIP_FILES or f.startswith("~$"):
                    continue
                p = os.path.join(dirpath, f)
                rel = os.path.join("GSK3B_mini_project",
                                   os.path.relpath(p, ROOT)).replace("\\", "/")
                z.write(p, rel)
                n += 1
                total += os.path.getsize(p)
    print(f"{n} files, {total/1e6:.1f} MB raw -> {OUT} "
          f"({os.path.getsize(OUT)/1e6:.1f} MB zipped)")


if __name__ == "__main__":
    main()
