#!/usr/bin/env python

# require: ffmpeg
# e.g: conda install ffmpeg -c menpo
# 
import sys
import subprocess

mov_in = sys.argv[1]
mov_out = sys.argv[2]

command = 'ffmpeg -i {} -pix_fmt rgb24 -r 10 -f gif {}'.format(mov_in, mov_out)
subprocess.call(command, shell=True)
