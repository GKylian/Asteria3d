from __future__ import print_function

import os
import cv2
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import glob

# Fixing random state for reproducibility
np.random.seed(19680801)


files = []

fig, ax = plt.subplots(figsize=(5, 5))
for i in range(50):  # 50 frames
    plt.cla()
    plt.imshow(np.random.rand(100, 100), interpolation='nearest')
    fname = '_tmp%03d.png' % i
    print('Saving frame', fname)
    plt.savefig(fname)
    files.append(fname)



print('Making movie animation.mpg - this may take a while')
#subprocess.call("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc "
#                "-lavcopts vcodec=wmv2 -oac copy -o animation.mpg", shell=True)
#os.system("ffmpeg -r 1 -i _tmp%01d.png -vcodec mpeg4 -y movie.mp4")

img_array = []
path = '*.png'
for filename in glob.glob(path):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)

out = cv2.VideoWriter('project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()


# cleanup
for fname in files:
    os.remove(fname)
