#!/usr/bin/python

import sys
import os
import subprocess as sp
import re
import shutil

home = os.getenv("HOME")

argv = sys.argv[1:]
if len(argv) < 4:
    print "usage: matrix_list mat_dir score_dir pic_dir"
    sys.exit()
mat_list_file = argv[0]
mat_dir = argv[1]
if not mat_dir.endswith("/"):
    mat_dir += "/"
score_dir = argv[2]
if not score_dir.endswith("/"):
    score_dir += "/"
pic_dir = argv[3]
if not pic_dir.endswith("/"):
    pic_dir += "/"


with open("riassuntivo.html", "w+") as html_fp:
    html_fp.write("<html><head><title>riassuntivo</title></head><body>\n")
    html_fp.write("<table>")
    with open(mat_list_file, "r") as mat_list_fp:
        for l in mat_list_fp:
            if "vertebrates" not in l:
                continue
            args = l.strip().split("\t")
            mat_file = mat_dir + args[0] + ".pfm"
            p = sp.Popen("Rscript logo_mix.R " +  mat_file, shell=True)
            p.wait()
            pic_file = pic_dir + args[0] + ".pfm.png"
            shutil.move(mat_file + ".png", pic_file)
            match = re.search("\.(\d)", args[0])
            score_file = score_dir + args[2] + "." + match.group(1) + ".score"
            first_col = "<tr><td colspan='3'>" + args[2] + "." + match.group(1) + "</td></tr>"

            pic_col = "<tr><td><img src='" + pic_file + "' width='360' height='180' /></td>"

            p = sp.Popen("Rscript mixture_sep.R 4 2 " + score_file, shell=True)
            p.wait()
            original_col = "<td><img src='" + score_file + "_original.png' width='700' height='350' /></td>"
            if os.path.isfile(score_file + "_mixture.png"):
                mixture_col = "<td><img src='" + score_file + "_mixture.png' width='700' height='350' /></td></tr>"
            else:
                mixture_col = "<td>NO CONVERGENCY</td></tr>"

            html_fp.write(first_col + "\n" + pic_col + original_col + mixture_col + "\n")
        html_fp.write("</table></body></html>")
