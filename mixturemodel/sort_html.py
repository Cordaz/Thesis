#!/usr/bin/python

import sys
import re

argv = sys.argv[1:]
if (not argv) or (len(argv) != 1):
   print "usage: html_file_path"
   sys.exit(1)
html_path = argv[0]
html_sorted = "".join(html_path.rsplit('.',1)[:-1]) + "_sorted.html"

# print html_sorted

content_dict = {}

with open(html_sorted, "w+") as html_out:
   html_out.write("<html><head><title>Riassuntivo ordinato</title></head><body><table>\n")
   with open(html_path, "r") as html_in:
      content = html_in.read()
      
   matches = re.findall("<tr><td colspan='3'>(\w+.\w+)</td></tr>\n<tr>(.*)</tr>", content)
   for match in matches:
      content_dict[match[0]] = match[1]
      # print match
   
   for key in sorted(content_dict):
      html_out.write("<tr><td colspan='3'>" + key + "</td></tr>\n<tr>" + content_dict[key] + "</tr>\n")
   
   html_out.write("</table></body></html>")
