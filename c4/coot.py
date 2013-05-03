
import os
import math
from subprocess import Popen, PIPE
import textwrap
import c4.utils

M_SQRT1_2 = 0.5**0.5


def basic_script(pdb, mtz, center, toward=None):
    text = """\
           # coot script generated by dimple
           pdb = "%s"
           mtz = "%s"

           #set_nomenclature_errors_on_read("ignore")""" % (pdb, mtz)
    if pdb: text += """
           molecule = read_pdb(pdb)"""
    if center: text += """
           set_rotation_centre(%g, %g, %g)
           set_zoom(30.)""" % center
    if toward: text += """
           set_view_quaternion(%g, %g, %g, %g)""" % as_quat(center, toward)
    if mtz: text += """
           map21 = make_and_draw_map(mtz, "2FOFCWT", "PH2FOFCWT", "", 0, 0)
           map11 = make_and_draw_map(mtz, "FOFCWT", "PHFOFCWT", "", 0, 1)"""
    return textwrap.dedent(text)


def as_quat(p1, p2):
    d = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
    length = math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
    d = (d[0]/length, d[1]/length, d[2]/length)
    # ref vector: 0 0 -1 (?)
    # cross product (a2 b3 - a3 b2, ..., ...)
    prod = (d[1], -d[0], 0)
    quat = (prod[0], prod[1], prod[2], 1-d[2])
    qlen = math.sqrt(sum(a*a for a in quat))
    return (quat[0]/qlen, quat[1]/qlen, quat[2]/qlen, quat[3]/qlen)

def mult_quat(q1, q2):
    x, y, z, w = q1
    ax, ay, az, aw = q2
    return (w*ax + x*aw + y*az - z*ay,
            w*ay + y*aw + z*ax - x*az,
            w*az + z*aw + x*ay - y*ax,
            w*aw - x*ax - y*ay - z*az)


def generate_r3d(pdb, mtz, center, blobname, cwd, toward=None):
    if toward is None:
        quat0 = (0., 0., 0., 1.)
    else:
        quat0 = as_quat(center, toward)

    quaternions = [quat0, mult_quat(quat0, (0., M_SQRT1_2, 0., M_SQRT1_2)),
                          mult_quat(quat0, (M_SQRT1_2, 0., 0., M_SQRT1_2))]

    coot_path = c4.utils.find_in_path("coot")
    if not coot_path:
        c4.utils.put_error("No coot, no pictures")
        return []
    coot_process = Popen([coot_path, "--python", "--no-graphics", "--no-guano"],
                         stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=cwd)
    # In coot, raster3d() creates file.r3d, make_image_raster3d() additionally
    # calls render program and opens image (convenient for testing)
    script = basic_script(pdb=pdb, mtz=mtz, center=center)
    basenames = []
    for n, quat in enumerate(quaternions):
        script += """
set_view_quaternion(%g, %g, %g, %g)""" % quat
        basename = "%sv%d" % (blobname, n+1)
        script += """
graphics_draw() # this is needed only for coot in --no-graphics mode
raster3d("%s.r3d")""" % basename
        basenames.append(basename)
    coot_process.communicate(input=script+"""
coot_real_exit(0)
""")
    return basenames

