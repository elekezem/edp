#!/usr/bin/env python2

import os

size = 3.61000000000000 * 1.007431
images = 1000

for i in range(0,images):
    depth = size / float(images) * float(i)
    pos = "0,%f,0" % (depth)
    filename = "./tmp/img_%i.png" % (i)
    os.system("./bin/edp -o %s -p %s -v 1,0,0 -w 0,0,1 -s 100" % (filename, pos))
