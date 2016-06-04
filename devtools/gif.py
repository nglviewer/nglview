from moviepy import editor

x = editor.VideoFileClip('nglview.mov')
x.write_gif("nglview.gif", opt='nq', fps=8)
