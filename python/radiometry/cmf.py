#!/usr/bin/env python3

"""
:mod: `cmf` -- Colour Matching Functions handling
=================================================
.. module:: cmf
    :synopsis: This modules allows conversion from spectral image to colour
               image given a set of colour matching functions.
.. moduleauthor:: Alban Fichet <alban.fichet@gmx.fr>
"""

import numpy as np
import scipy.interpolate

class CMF:
    def __init__(self, filename, sampling=2):
        """
        Initialise an object for colour conversion.

        :param filename CSV file with the colour matching functions.
        :param sampling Used to lower the accuracy but gain speed.
        """
        x_bar = []
        y_bar = []
        z_bar = []
        
        # Read the CSV file
        with open(filename) as f:    
            for l in f:
                wl, x, y, z = [float(el) for el in l.split(',')]
                x_bar.append([wl, x])
                y_bar.append([wl, y])
                z_bar.append([wl, z])

        x_bar = np.array(x_bar)
        y_bar = np.array(y_bar)
        z_bar = np.array(z_bar)
        
        # Interpolate every 1nm to ease computations
        self.wl_start = np.min(x_bar[:, 0])
        self.wl_end   = np.max(x_bar[:, 0])
        # self.wl       = np.linspace(self.wl_start, self.wl_end, num=int(self.wl_end - self.wl_start + 1))

        # let's be a little bit less hardcore here...
        self.wl       = np.arange(self.wl_start, self.wl_end, sampling)

        x_bar_y = np.interp(self.wl, x_bar[:, 0], x_bar[:, 1], left=0, right=0)
        y_bar_y = np.interp(self.wl, y_bar[:, 0], y_bar[:, 1], left=0, right=0)
        z_bar_y = np.interp(self.wl, z_bar[:, 0], z_bar[:, 1], left=0, right=0)

        self.x_bar = np.stack((self.wl, x_bar_y), axis=1)
        self.y_bar = np.stack((self.wl, y_bar_y), axis=1)
        self.z_bar = np.stack((self.wl, z_bar_y), axis=1)
        

    def get_xyz_emissive_img(self, image_wl, spectral_image):
        """
        Converts an emissive spectral image stored in a numpy array 
        (y, x, bands) to a colour image.

        :param image_wl list containing the wavelength values in nanometers
               for the provided spectral image
        :param spectral_image numpy array storing the spectral image
        """
        interp_image_f = scipy.interpolate.interp1d(image_wl, spectral_image, bounds_error=False, fill_value=(0, 0))
        interp_image = np.maximum(interp_image_f(self.wl), 0)

        delta = self.wl[1] - self.wl[0]

        s_x = np.sum(interp_image * self.x_bar[:, 1], axis=-1) * delta
        s_y = np.sum(interp_image * self.y_bar[:, 1], axis=-1) * delta
        s_z = np.sum(interp_image * self.z_bar[:, 1], axis=-1) * delta

        return np.dstack((s_x, s_y, s_z))


    def get_xyz_reflective_img(self, wavelength_illu, spectrum_illu, image_wl, spectral_image):
        """
        Converts an reflective spectral image stored in a numpy array 
        (y, x, bands) to a colour image.

        :param wavelnegth_illu list of the wavlengths of the illuminant 
               spectrum
        :param spectrum_illu list of the corresponding radiance of the 
               illuminant spectrum
        :param image_wl list containing the wavelength values in nanometers
               for the provided spectral image
        :param spectral_image numpy array storing the spectral image
        """
        illu_values = np.interp(self.wl, wavelength_illu, spectrum_illu, left=0, right=0)

        interp_image_f = scipy.interpolate.interp1d(image_wl, spectral_image, bounds_error=False, fill_value=(0, 0))
        interp_image = np.maximum(interp_image_f(self.wl), 0)

        delta = self.wl[1] - self.wl[0]

        Y_illu = np.sum(illu_values * self.y_bar[:, 1], axis=-1) * delta

        s_x = np.sum(interp_image * illu_values * self.x_bar[:, 1], axis=-1) / Y_illu * delta
        s_y = np.sum(interp_image * illu_values * self.y_bar[:, 1], axis=-1) / Y_illu * delta
        s_z = np.sum(interp_image * illu_values * self.z_bar[:, 1], axis=-1) / Y_illu * delta

        return np.dstack((s_x, s_y, s_z))