from glob import glob

pyx_files = glob("*.pyx")
pyx_files.remove("__init__.pyx")
pyx_files.remove("Analysis.pyx")

for pyx in sorted(pyx_files):
    pyx = pyx.replace(".pyx", "")
    txt = "from pytraj.analyses.%s import %s" % (pyx, pyx)
    print(txt)
