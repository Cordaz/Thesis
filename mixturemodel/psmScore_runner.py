import sys
import subprocess as sp
import re

argv = sys.argv[1:]

if not argv:
    print "usage: python psmScore_runner.py matrix_list_file matrces_dir out_dir fasta_file"
    sys.exit(1)

mat_list_file = argv[0]
del argv[0]
if not argv:
    print "usage: python psmScore_runner.py matrix_list_file matrces_dir out_dir fasta_file"
    sys.exit(1)

mat_dir = argv[0]
del argv[0]
if not argv:
    print "usage: python psmScore_runner.py matrix_list_file matrces_dir out_dir fasta_file"
    sys.exit(1)
if not mat_dir.endswith("/"):
    mat_dir += "/"

out_dir = argv[0]
del argv[0]
if not argv:
    print "usage: python psmScore_runner.py matrix_list_file matrces_dir out_dir fasta_file"
    sys.exit(1)
if not out_dir.endswith("/"):
    out_dir += "/"

fa_file = argv[0]
del argv[0]
if argv:
    print "usage: python psmScore_runner.py matrix_list_file matrces_dir out_dir fasta_file"
    sys.exit(1)

with open(mat_list_file, "r") as mat_list_fp:
    for l in mat_list_fp:
        if "vertebrates" not in l:
            continue
        args = l.strip().split("\t")
        mat_file = mat_dir + args[0] + ".pfm"
        match = re.search("\.(\d)", args[0])
        out_file = out_dir + args[2] + "." + match.group(1) + ".score"
        if '(' in out_file:
            out_file = out_file.replace('(', '_')
        if ')' in out_file:
            out_file = out_file.replace(')', '')

        p = sp.Popen("python psmScore.py " + mat_file + " " + fa_file + " " + out_file, shell=True)
        p.wait()
