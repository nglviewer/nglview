#!/usr/bin/env python

import subprocess
# require: ffmpeg
# e.g: conda install ffmpeg -c menpo
# 
import sys

mov_in = sys.argv[1]
mov_out = sys.argv[2]

command = f'ffmpeg -i {mov_in} -pix_fmt rgb24 -r 10 -f gif {mov_out}'
subprocess.call(command, shell=True)
