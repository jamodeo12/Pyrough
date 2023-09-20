# usage : python Surface_Analysis.py 'image.png' size xmin xmax
# image.png is the image to analyse, size is the lateral length of the picture and xmin and xmax are the height distribution limit values

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Func_pyrough as func
from scipy import stats
import cv2
import sys

def hann(image):
    size = np.shape(image)
    Lx = size[0]
    Ly = size[1]
    for x in range(Lx):
        for y in range(Ly):
            if (((x - 0.5 * Lx) ** 2 + (y - 0.5 * Ly) ** 2) < (0.5 * min(Lx, Ly)) ** 2):
                w = (((3 * np.pi / 8) - (2 / np.pi)) ** (-0.5)) * (1 + np.cos(
                    (2 * np.pi * np.sqrt((x - 0.5 * Lx) ** 2 + (y - 0.5 * Ly) ** 2)) / (min(Lx, Ly))))
                image[x, y] = w * image[x, y]
            else:
                image[x, y] = 0
    return (image)

print("====== > Running surface analysis algorithm")

fig = plt.figure()

# Parameters
D = cv2.imread(sys.argv[1], cv2.IMREAD_GRAYSCALE)
size_init = float(sys.argv[2])
zmax = float(sys.argv[2])
zmin = float(sys.argv[3])
npix = np.min(D.shape)
size = size_init * (npix / np.max(D.shape))

# Square the pic
D = D[0:npix, 0:npix] - np.mean(D)
size_per_pixel = size / np.shape(D)[0]

# Wavevectors range
qmax = 2 * np.pi / size_per_pixel
qmin = 2 * np.pi / (0.5 * size)

# Show the image
ax1 = fig.add_subplot(1, 3, 1)
ax1.axis('off')
ax1.imshow(D, cmap='gray')

# Apply a Hanning window to the grayscale image, and rescale as a function of zmax and zmin
image = hann(D)
image = (image - np.min(image)) / (np.max(image) - np.min(image)) * (zmax - zmin) + zmin

# print('Mean line : ', round(np.mean(D), 2))
# print('Deviation : ', func.sigma(D))
# print('RMS : ', func.rms_calc(D))
# print('Skewness : ', func.sk(D))
# print('Kurtosis : ', func.Kurto(D))

print("====== > Generation of the surface power spectrum ...")

# Fourier transform and power spectral density
fourier_image = np.fft.fftn(image)
fourier_amplitudes = np.abs(fourier_image) ** 2

# Frequencies correspondng to signal PSD
kfreq = np.fft.fftfreq(npix) * npix
kfreq2D = np.meshgrid(kfreq, kfreq)
knrm = np.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)

knrm = knrm.flatten()
fourier_amplitudes = fourier_amplitudes.flatten()

# PSD generation
kbins = np.arange(0.5, npix // 2 + 1, 1.)
kvals = 0.5 * (kbins[1:] + kbins[:-1])
Abins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes,
                                     statistic="mean",
                                     bins=kbins)
Abins *= np.pi * (kbins[1:] ** 2 - kbins[:-1] ** 2)

# Keep only the wavevectors range
index = [i for i, v in enumerate(kvals) if ((v > qmax) or (v < qmin))]
Abins[index] = 0

# Linear regression
index_reg = np.where(Abins != 0)[0]
low = np.min(index_reg) + 1
sup = np.max(index_reg)

m, b = np.polyfit(np.log(kvals[low:sup]), np.log(Abins[low:sup]), 1)
H = -1 * (0.5 * m + 1)
print("====== > H = "+str(round(H,2)))

ax2 = fig.add_subplot(1, 3, 2)
ax2.loglog(kvals, Abins, color='r', linewidth=2.5, label="$PSD$")
ax2.loglog(kvals[index_reg], np.exp(b) * np.power(kvals[index_reg], m), '--', color='k', linewidth=2.5,
           label="$H = $" + str(round(H, 2)))
ax2.legend()
ax2.set_xlabel("Spatial Frequency $k [m^{-1}]$")
ax2.set_ylabel("PSD $P(k)$")

print("====== > Construction of equivalent rough surface ...")
x = np.linspace(0, 1, 200)
y = x
xv, yv = np.meshgrid(x, y)
Z = func.rough(xv, yv, H, 1, 50, 50)

ax3 = fig.add_subplot(1, 3, 3, projection='3d')
ax3.grid(False)
ax3.axis('off')
ax3.scatter3D(xv, yv, Z, c=Z, cmap='jet', s=1)
ax3.view_init(90, -90)

plt.show()

df = pd.DataFrame(Z)
df.to_csv('Rough_data.csv')
