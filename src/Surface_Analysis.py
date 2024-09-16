# usage : python Surface_Analysis.py 'image.png' size xmin xmax
# image.png is the image to analyse, size is the lateral length of the picture and xmin and xmax are the height distribution limit values

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

print("====== > Running surface analysis algorithm ...")

# Parameters
D = cv2.imread(sys.argv[1], cv2.IMREAD_GRAYSCALE)
size_init = float(sys.argv[2])
zmin = float(sys.argv[3])
zmax = float(sys.argv[4])
npix = np.min(D.shape)
size = size_init * (npix / np.max(D.shape))

# Square the pic
D = D[0:npix, 0:npix] - np.mean(D)
size_per_pixel = size / np.shape(D)[0]

#print('RMS : ', func.rms_calc(func.rescale(D, [zmin, zmax])))

# Wavevectors range
qmax = 2 * np.pi / size_per_pixel
qmin = 2 * np.pi / (0.5 * size)

# Show the image
fig1 = plt.figure()
ax1 = fig1.add_subplot()
ax1.axis('off')
ax1.imshow(D, cmap='gray')

# Apply a Hanning window to the grayscale image, and rescale as a function of zmax and zmin
image = hann(D)
image = (image - np.min(image)) / (np.max(image) - np.min(image)) * (zmax - zmin) + zmin

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
Abins = np.delete(Abins, index)
kvals = np.delete(kvals, index)

print("Scanned spatial frequencies : ", kvals)

m, b = np.polyfit(np.log(kvals), np.log(Abins), 1)
H = -1 * (0.5 * m + 1)
eta = (H - 1)/2
#print("====== > eta = "+str(round(eta,2)))
print("====== > Extraction of rough surface statistical parameters ...")
print("RMS : {}".format(func.rms_calc(func.rescale(D, [zmin, zmax]))))
print("eta : {}, H = {}".format(eta, H))
print("A : {}".format(kvals[-1]))
print("B : {}".format(kvals[-1]))

fig2 = plt.figure()
ax2 = fig2.add_subplot()
ax2.loglog(kvals, Abins, color='r', linewidth=2.5, label="$PSD$")
ax2.loglog(kvals, np.exp(b) * np.power(kvals, m), '--', color='k', linewidth=2.5,
           label="$eta = $" + str(round(eta, 2)))
ax2.legend()
ax2.set_xlabel("Wave vector $q \ [nm^{-1}]$")
ax2.set_ylabel("PSD $C^{2D}(q) \ [nm^{4}]$")

print("====== > Construction of equivalent rough surface ...")
x = np.linspace(0, 1, 200)
y = x
xv, yv = np.meshgrid(x, y)
Z = func.rough(xv, yv, eta, func.rms_calc(D), int(np.max(kvals)), int(np.max(kvals)))

fig3 = plt.figure()
ax3 = fig3.add_subplot(projection='3d')
ax3.grid(False)
ax3.axis('off')
ax3.scatter3D(xv, yv, Z, c=Z, cmap='jet', s=1)
ax3.view_init(90, -90)

plt.show()

np.savetxt("Rough_data.csv", Z, delimiter=",")

print("====== > Data saved in Rough_data.csv")

