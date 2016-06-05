from moviepy import editor

x = editor.VideoFileClip('nglview.mov')
x.write_gif("nglview.gif", fps=8, opt='nq')
