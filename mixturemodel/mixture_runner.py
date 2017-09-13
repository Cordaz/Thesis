import sys
import subprocess as sp

argv = sys.argv[1:]

if not argv:
    print "usage: python mixture_runner.py score_list_file score_dir"
    sys.exit(1)

score_list_file = argv[0]
del argv[0]
if not argv:
    print "usage: python mixture_runner.py score_list_file score_dir"
    sys.exit(1)

score_dir = argv[0]
del argv[0]
if not score_dir.endswith("/"):
    score_dir += "/"

if argv:
    print "usage: python mixture_runner.py score_list_file score_dir out_dir"
    sys.exit(1)

with open(score_list_file, "r") as score_list_fp:
    for l in score_list_fp:
        l.strip()
        score_file = score_dir + l

        p = sp.Popen("Rscript mixture.R 4 2 " + score_file, shell=True)
        p.wait()
