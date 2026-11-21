"""Make the installed MODELLER 10.7 importable from the system Python 3.13.

MODELLER ships its own extension module per Python version under
lib/x86_64-w64/python3.13 and its shared libraries under lib/x86_64-w64.
Import this module before `import modeller`.
"""
import os, sys, glob

ROOT = os.environ.get("MODELLER_ROOT", r"C:\Program Files\Modeller10.7")


def setup(root=ROOT):
    if not os.path.isdir(root):
        cands = sorted(glob.glob(r"C:\Program Files\Modeller*"))
        if not cands:
            raise RuntimeError("MODELLER installation not found")
        root = cands[-1]
    pyver = f"python{sys.version_info.major}.{sys.version_info.minor}"
    extdir = os.path.join(root, "lib", "x86_64-w64", pyver)
    libdir = os.path.join(root, "lib", "x86_64-w64")
    if not os.path.isdir(extdir):
        raise RuntimeError(f"no MODELLER extension for {pyver} in {libdir}")
    for p in (os.path.join(root, "modlib"), extdir):
        if p not in sys.path:
            sys.path.insert(0, p)
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(libdir)
    # The licence key is NOT stored here. MODELLER writes it into
    # modlib/modeller/config.py when you install it, and reads it from there;
    # set the KEY_MODELLER environment variable only if you need to override
    # that. Get your own free academic key at
    # https://salilab.org/modeller/registration.html
    return root


setup()
