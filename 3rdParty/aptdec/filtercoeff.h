/* 
 *  This file is part of Aptdec.
 *  Copyright (c) 2004-2009 Thierry Leconte (F4DWV), Xerbo (xerbo@protonmail.com) 2019
 *
 *  Aptdec is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

// IQ finite impulse response filter
#define IQFilterLen 32
const float iqfilter[IQFilterLen] = { 0.0205361, 0.0219524, 0.0235785, 0.0254648, 0.0276791, 0.0303152,
0.0335063, 0.0374482, 0.0424413, 0.0489708, 0.0578745, 0.0707355, 0.0909457, 0.127324, 0.212207, 0.63662,
-0.63662, -0.212207, -0.127324, -0.0909457, -0.0707355, -0.0578745, -0.0489708, -0.0424413, -0.0374482,
-0.0335063, -0.0303152, -0.0276791, -0.0254648, -0.0235785, -0.0219524, -0.0205361 };

// Pattern of a NOAA sync line
#define SyncFilterLen 32
const float Sync[SyncFilterLen] = { -14, -14, -14, 18, 18, -14, -14, 18, 18, -14, -14, 18, 18, -14, -14,
18, 18, -14, -14, 18, 18, -14, -14, 18, 18, -14, -14, 18, 18, -14, -14, -14 };

// Gaussian finite impulse response compensation filter
#define RSFilterLen 437
const float rsfilter[RSFilterLen] = { -3.37279e-04, -8.80292e-06, -3.96418e-04, -1.78544e-04, -5.27511e-04,
-3.75376e-04, -6.95337e-04, -5.93148e-04, -8.79730e-04, -8.15327e-04, -1.05669e-03, -1.01377e-03,
-1.19836e-03, -1.15443e-03, -1.26937e-03, -1.20955e-03, -1.23904e-03, -1.15302e-03, -1.08660e-03,
-9.64235e-04, -8.02450e-04, -6.46202e-04, -3.95376e-04, -2.18096e-04, 1.11906e-04, 2.89567e-04,
6.67167e-04, 8.19039e-04, 1.21725e-03, 1.30556e-03, 1.69365e-03, 1.68588e-03, 2.03277e-03, 1.90159e-03,
2.18455e-03, 1.90833e-03, 2.12100e-03, 1.69052e-03, 1.77484e-03, 1.42542e-03, 1.18292e-03, 8.66979e-04,
5.54161e-04, 2.15793e-04, -1.11623e-04, -4.35173e-04, -7.27194e-04, -9.91551e-04, -1.20407e-03,
-1.37032e-03, -1.46991e-03, -1.51120e-03, -1.48008e-03, -1.39047e-03, -1.23115e-03, -1.02128e-03,
-7.60099e-04, -4.68008e-04, -1.46339e-04, 1.80867e-04, 5.11244e-04, 8.19243e-04, 1.09739e-03, 1.32668e-03,
1.50632e-03, 1.61522e-03, 1.66246e-03, 1.62390e-03, 1.52430e-03, 1.34273e-03, 1.10736e-03, 8.10335e-04,
4.76814e-04, 1.13622e-04, -2.64150e-04, -6.26595e-04, -9.95436e-04, -1.27846e-03, -1.54080e-03,
-1.74292e-03, -1.86141e-03, -1.89318e-03, -1.83969e-03, -1.69770e-03, -1.47938e-03, -1.18696e-03,
-8.37003e-04, -4.39507e-04, -1.56907e-05, 4.19904e-04, 8.43172e-04, 1.23827e-03, 1.58411e-03,
1.86382e-03, 2.06312e-03, 2.17177e-03, 2.18121e-03, 2.08906e-03, 1.89772e-03, 1.61153e-03,
1.24507e-03, 8.13976e-04, 3.29944e-04, -1.74591e-04, -6.83619e-04, -1.17826e-03, -1.61659e-03,
-2.00403e-03, -2.29070e-03, -2.49179e-03, -2.56546e-03, -2.53448e-03, -2.37032e-03, -2.10060e-03,
-1.72140e-03, -1.24542e-03, -7.15425e-04, -1.24964e-04, 4.83736e-04, 1.08328e-03, 1.64530e-03,
2.14503e-03, 2.55400e-03, 2.85589e-03, 3.02785e-03, 3.06271e-03, 2.95067e-03, 2.69770e-03,
2.30599e-03, 1.79763e-03, 1.18587e-03, 5.04003e-04, -2.23591e-04, -9.57591e-04, -1.66939e-03,
-2.31717e-03, -2.87636e-03, -3.31209e-03, -3.60506e-03, -3.73609e-03, -3.69208e-03, -3.44913e-03,
-3.06572e-03, -2.50229e-03, -1.80630e-03, -1.00532e-03, -1.22305e-04, 7.83910e-04, 1.69402e-03,
2.53826e-03, 3.30312e-03, 3.91841e-03, 4.38017e-03, 4.63546e-03, 4.68091e-03, 4.50037e-03,
4.09614e-03, 3.47811e-03, 2.67306e-03, 1.70418e-03, 6.20542e-04, -5.36994e-04, -1.70981e-03,
-2.84712e-03, -3.88827e-03, -4.78659e-03, -5.48593e-03, -5.95049e-03, -6.14483e-03, -6.05118e-03,
-5.65829e-03, -4.97525e-03, -4.01796e-03, -2.82224e-03, -1.43003e-03, 1.00410e-04, 1.71169e-03,
3.31983e-03, 4.87796e-03, 6.23237e-03, 7.31013e-03, 8.20642e-03, 8.67374e-03, 8.77681e-03,
8.43444e-03, 7.66794e-03, 6.46827e-03, 4.87294e-03, 2.92923e-03, 6.98913e-04, -1.72126e-03,
-4.24785e-03, -6.75380e-03, -9.13309e-03, -1.12532e-02, -1.30038e-02, -1.42633e-02, -1.49338e-02,
-1.49145e-02, -1.41484e-02, -1.25761e-02, -1.01870e-02, -6.97432e-03, -2.97910e-03, 1.75386e-03,
7.11899e-03, 1.30225e-02, 1.93173e-02, 2.58685e-02, 3.24965e-02, 3.90469e-02, 4.53316e-02,
5.11931e-02, 5.64604e-02, 6.09924e-02, 6.46584e-02, 6.73547e-02, 6.90049e-02, 6.97096e-02,
6.90049e-02, 6.73547e-02, 6.46584e-02, 6.09924e-02, 5.64604e-02, 5.11931e-02, 4.53316e-02,
3.90469e-02, 3.24965e-02, 2.58685e-02, 1.93173e-02, 1.30225e-02, 7.11899e-03, 1.75386e-03,
-2.97910e-03, -6.97432e-03, -1.01870e-02, -1.25761e-02, -1.41484e-02, -1.49145e-02, -1.49338e-02,
-1.42633e-02, -1.30038e-02, -1.12532e-02, -9.13309e-03, -6.75380e-03, -4.24785e-03, -1.72126e-03,
6.98913e-04, 2.92923e-03, 4.87294e-03, 6.46827e-03, 7.66794e-03, 8.43444e-03, 8.77681e-03,
8.67374e-03, 8.20642e-03, 7.31013e-03, 6.23237e-03, 4.87796e-03, 3.31983e-03, 1.71169e-03,
1.00410e-04, -1.43003e-03, -2.82224e-03, -4.01796e-03, -4.97525e-03, -5.65829e-03, -6.05118e-03,
-6.14483e-03, -5.95049e-03, -5.48593e-03, -4.78659e-03, -3.88827e-03, -2.84712e-03, -1.70981e-03,
-5.36994e-04, 6.20542e-04, 1.70418e-03, 2.67306e-03, 3.47811e-03, 4.09614e-03, 4.50037e-03,
4.68091e-03, 4.63546e-03, 4.38017e-03, 3.91841e-03, 3.30312e-03, 2.53826e-03, 1.69402e-03,
7.83910e-04, -1.22305e-04, -1.00532e-03, -1.80630e-03, -2.50229e-03, -3.06572e-03, -3.44913e-03,
-3.69208e-03, -3.73609e-03, -3.60506e-03, -3.31209e-03, -2.87636e-03, -2.31717e-03, -1.66939e-03,
-9.57591e-04, -2.23591e-04, 5.04003e-04, 1.18587e-03, 1.79763e-03, 2.30599e-03, 2.69770e-03,
2.95067e-03, 3.06271e-03, 3.02785e-03, 2.85589e-03, 2.55400e-03, 2.14503e-03, 1.64530e-03,
1.08328e-03, 4.83736e-04, -1.24964e-04, -7.15425e-04, -1.24542e-03, -1.72140e-03, -2.10060e-03,
-2.37032e-03, -2.53448e-03, -2.56546e-03, -2.49179e-03, -2.29070e-03, -2.00403e-03, -1.61659e-03,
-1.17826e-03, -6.83619e-04, -1.74591e-04, 3.29944e-04, 8.13976e-04, 1.24507e-03, 1.61153e-03,
1.89772e-03, 2.08906e-03, 2.18121e-03, 2.17177e-03, 2.06312e-03, 1.86382e-03, 1.58411e-03,
1.23827e-03, 8.43172e-04, 4.19904e-04, -1.56907e-05, -4.39507e-04, -8.37003e-04, -1.18696e-03,
-1.47938e-03, -1.69770e-03, -1.83969e-03, -1.89318e-03, -1.86141e-03, -1.74292e-03, -1.54080e-03,
-1.27846e-03, -9.95436e-04, -6.26595e-04, -2.64150e-04, 1.13622e-04, 4.76814e-04, 8.10335e-04,
1.10736e-03, 1.34273e-03, 1.52430e-03, 1.62390e-03, 1.66246e-03, 1.61522e-03, 1.50632e-03,
1.32668e-03, 1.09739e-03, 8.19243e-04, 5.11244e-04, 1.80867e-04, -1.46339e-04, -4.68008e-04,
-7.60099e-04, -1.02128e-03, -1.23115e-03, -1.39047e-03, -1.48008e-03, -1.51120e-03, -1.46991e-03,
-1.37032e-03, -1.20407e-03, -9.91551e-04, -7.27194e-04, -4.35173e-04, -1.11623e-04, 2.15793e-04,
5.54161e-04, 8.66979e-04, 1.18292e-03, 1.42542e-03, 1.77484e-03, 1.69052e-03, 2.12100e-03,
1.90833e-03, 2.18455e-03, 1.90159e-03, 2.03277e-03, 1.68588e-03, 1.69365e-03, 1.30556e-03,
1.21725e-03, 8.19039e-04, 6.67167e-04, 2.89567e-04, 1.11906e-04, -2.18096e-04, -3.95376e-04,
-6.46202e-04, -8.02450e-04, -9.64235e-04, -1.08660e-03, -1.15302e-03, -1.23904e-03, -1.20955e-03,
-1.26937e-03, -1.15443e-03, -1.19836e-03, -1.01377e-03, -1.05669e-03, -8.15327e-04, -8.79730e-04,
-5.93148e-04, -6.95337e-04, -3.75376e-04, -5.27511e-04, -1.78544e-04, -3.96418e-04, -8.80292e-06,
-3.37279e-04 };