#!/usr/bin/env ccp4-python

import os
import sys

sys.path.append(os.path.join(os.environ["CCP4"], "share", "smartie"))
import qtrapi

# monkey-patch put() before dimple is imported
import c4.utils
orig_put = c4.utils.put
def our_put(text):
    orig_put(text)
    qtrapi_text[0] += text
    for n in (1,2):
        if ('-> blob%dv3.png' % n) in text:
            add_blob(n)

    report.flush()
c4.utils.put = our_put

def add_blob(n):
    dfiles = qtrapi.Decorations(title=("Blob %d" % n))
    for i in range(1,4):
        path = os.path.join(jobdir, "blob%dv%d.png" % (n, i))
        dfiles.append(qtrapi.File(key="DECORATION", type="image", path=path))
        path = os.path.join(jobdir, "blob%d-coot.py" % n)
    coot_view = qtrapi.File(key="COOT_VIEW", type="coot_view", path=path,
                            title="View in")
    vfiles = qtrapi.Files()
    vfiles.append(coot_view)
    section_output.extend([dfiles, vfiles])

import dimple
opt = dimple.parse_dimple_commands()
jobdir = os.path.abspath(opt.output_dir)

print "_JOB_DIRECTORY:", jobdir
sys.stdout.flush()

input_files = qtrapi.Files()
input_files.extend([
qtrapi.File(key="HKLIN", title="Reflection Data", type="hkl:hkl",
            path=os.path.join(jobdir, opt.mtz)),
qtrapi.File(key="XYZIN", title="Input Model", type="xyz",
            path=os.path.join(jobdir, opt.pdb))
])

qtrapi_text = qtrapi.Text("")
qtrapi_text.append("")
section_info = qtrapi.Section("Output")
section_info.append(qtrapi_text)

section_input = qtrapi.Section("Input Files")
section_input.append(input_files)

section_output = qtrapi.Section("Output Files")

report = qtrapi.Report(title="Dimple", path=jobdir)
report.extend([section_info, section_input, section_output])
report.flush()

dimple.main()
