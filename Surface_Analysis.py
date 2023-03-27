import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Sources import Func_pyrough as func
from scipy import stats
import cv2
import sys

fig = plt.figure()

D = cv2.imread(sys.argv[1], 0)

scale = [float(sys.argv[2]), float(sys.argv[3])]
D = func.rescale(D, scale)

print('Mean line : ', round(np.mean(D),2))
print('Deviation : ', func.sigma(D))
print('RMS : ', func.RMS(D))
print('Skewness : ', func.sk(D))
print('Kurtosis : ', func.Kurto(D))

ax1 = fig.add_subplot(1, 3, 1)
ax1.axis('off')
ax1.imshow(D, cmap='gray')

L = np.shape(D)[1]
data = []
for i in range(L):
	data.append(D[i])
data = np.asarray(data)/np.max(data)

x = np.arange(0,int(np.size(data,0)),1)
y = np.arange(0,int(np.size(data,1)),1)

z = data
npix = data.shape[0]

#taking the fourier transform
fourier_image = np.fft.fft2(data)

#Get power spectral density
fourier_amplitudes = np.abs(fourier_image)**2

#calculate sampling frequency fs (physical distance between pixels)
fs = 1/npix
freq_shifted = fs/2 * np.linspace(-1,1,npix)
freq = fs/2 * np.linspace(0,1,int(npix/2))

#constructing a wave vector array
## Get frequencies corresponding to signal PSD
kfreq = np.fft.fftfreq(npix) * npix

kfreq2D = np.meshgrid(kfreq, kfreq)

knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)
knrm = knrm.flatten()

fourier_amplitudes = fourier_amplitudes.flatten()

#creating the power spectrum
kbins = np.arange(0.5, npix//2+1, 1.)
kvals = 0.5 * (kbins[1:] + kbins[:-1])

Abins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes,
                                     statistic = "mean",
                                     bins = kbins)
Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)

X = kvals
Y = Abins

sup = 0
low = 125
m, b = np.polyfit(np.log(X[sup:low]), np.log(Y[sup:low]), 1)

H = -0.5*m-1

print("Plotting power spectrum of surface ...")
#fig = plt.figure()
ax2 = fig.add_subplot(1, 3, 2)
ax2.loglog(X, Y)
ax2.loglog(X, np.exp(b)*np.power(X,m), label = 'H = '+str(round(H,2)))
ax2.legend()
#ax2.set_xlabel("Spatial Frequency $k$ [m-1]")
#ax2.set_ylabel("PSD $P(k)$")
#plt.tight_layout()

print("Construction of equivalent rough surface ...")
x = np.linspace(0,1,300)
y = x
xv, yv = np.meshgrid(x,y)
Z = func.rough(xv, yv, H, 1, 50,50)

ax3 = fig.add_subplot(1, 3, 3, projection='3d')
ax3.grid(False)
ax3.axis('off')
ax3.scatter3D(xv, yv, Z, c=Z, cmap='jet')
ax3.view_init(90, -90)

plt.show()

df = pd.DataFrame(Z)
df.to_csv('Rough_data.csv')