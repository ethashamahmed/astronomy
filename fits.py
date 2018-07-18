#######################################################################
# Ethasham Ahmed, 2018
# With contributions from AJ, Despina, Conor and Isacc.
# King's College London
#
# This is a program for reading from fits files, stacking multiple fits files,
# and plotting spectra using the data.
#
# The fits files used in the code are from Sirius, Orion, And Jupiter.
# The specta for Sirius have also been calibrated by comparing the postion of
# known Hydrogen peaks to determine the wavelengths.
#######################################################################

#Importing necessary modules
#Command for installing astropy: python -m pip install astropy --user
from __future__ import print_function
from astropy.io import fits
from astropy.nddata import Cutout2D
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os 

# Enter the location of your fits files here.
location = os.getcwd()

#The main function of the program that calls all the necessary classes and functions.
def fitsfunction():

	print ('Choose between 1, 2, 3, or 4: ')
	print ('1)Sirius')
	print ('2)Jupiter')
	print ('3)Orion')
	print ('4)Use any fits files')
	fitsobject = int(input())

	#The following if statements decide which functions to call based on user input.
	if fitsobject == 1:
		fitsfile = int(input('Enter Sirius bin number: '))

		if fitsfile == 2:
			image = siriusbin2().initialise()
			xvalues = siriusbin2().calib_sb2(image)
			siriusbin2().plot_sb2(xvalues,image)

		elif fitsfile == 3:
			image = siriusbin3().initialise()
			xvalues = siriusbin3().calib_sb3(image)
			siriusbin3().plot_sb3(xvalues,image)

	elif fitsobject == 2:
		fitsfile = int(input('Enter Jupiter bin number: '))

		if fitsfile == 2:
			image = jupiterbin2().initialise()
			jupiterbin2().plot_jb2(image)

		elif fitsfile == 3:
			image = jupiterbin3().initialise()
			jupiterbin3().plot_jb3(image)
			jupiterbin3().r_velocity()

	elif fitsobject == 3:
		fitsfile = int(input('Enter Orion bin number: '))

		if fitsfile == 2:
			image = orionbin2().initialise()
			orionbin2().plot_ob2(image)

		elif fitsfile == 3:
			image = orionbin3().initialise()
			orionbin3().plot_ob3(image)

	elif fitsobject == 4:
		final_image = fitsimages().initialise()
		print ('Do you want to: ')
		print ('1) Show image')
		print ('2) Plot Histogram')
		print ('3) plot Spectra')
		choice = int(input())
		if choice == 1:
			fitsimages().show_image(final_image)
			fitsimages().write_fits(final_image)
		elif choice == 2:
			fitsimages().plot_hist(final_image)
			fitsimages().write_fits(final_image)
		elif choice == 3:
			fitsimages().plot_spectra(final_image)
			fitsimages().write_fits(final_image)

#######################################################################
# This is a generalized class for opening, stacking and plotting spectra
# for fits files given by the user. There is also a function to view the
# image.
#######################################################################
class fitsimages():
	def initialise(self):
		# Opening Fits files and loading the image data
		numfits = int(input('How many fits files do you want to use?: '))
		image = []

		for n in range(1,numfits+1):
			print (n, ') Enter file name: ')
			image_file = input()
			hdu_list = fits.open(image_file)
			hdu_list.info()
			image_data = hdu_list[0].data
			print ('Image data:')
			print (type(image_data))
			print ((image_data.shape))
			image.append(image_data)
			hdu_list.close()

		# stacking the images
		final_image = np.sum(image, axis=0)

		return final_image

	def show_image(self, final_image):
		from matplotlib.colors import LogNorm
		# The images collected were reversed hence they are flipped here
		final_image = np.flip(final_image, axis=1)
		plt.imshow(final_image, cmap='gray', norm=LogNorm())
		plt.colorbar()
		plt.savefig('stackedimage.png')
		plt.show()

	def plot_hist(self, final_image):
		numbins = int(input('Enter number of bins: '))
		histogram = plt.hist(final_image.flatten(), numbins)
		plt.savefig('histogram.png')
		plt.show()

	def plot_spectra(self, final_image):
		binim = np.sum(final_image, axis=0)
		binim = np.flip(binim, axis=0)

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Pixel number')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_xlim([0, len(binim)+1])
		ax.set_ylim([min(binim), max(binim)])
		plt.plot(binim, color='k')
		plt.savefig('spectra.png')
		plt.show()

	def write_fits(self, final_image):
		# writing new fits files using stacked data
		outfile = 'stackedimage.fits'
		hdu = fits.PrimaryHDU(final_image)
		hdu.writeto(outfile, clobber=True)
#######################################################################
# The following classes correspond to a particuler set of fitsfiles.
# Each class is named after the object the fits files data was collected from.
# Comments have been added to the class for Sirius Bin 2 to explain what the
# functions do.
# The rest of the classes contain similar functions.
#
# The class jupiterbin3 contains a function to calculate the rotational velocity
# of Jupiter.
#######################################################################
class siriusbin2():
	def initialise(self):
		# This function opens the fits files and stacks them.

		# All the names of the necessary fits files are stored in a list.
		image_list = [location+'\siriusbin2\siriusbin2_'+str(n)+'.fit' for n in range(421,440)]
		# The fits files are stacked into one 2-d array.
		spectra = [fits.getdata(image) for image in image_list]
		# This intigrates the array.
		bin2 = np.sum(spectra, axis=0)
		# The data is reversed so this flips it around.
		bin2 = np.flip(bin2, axis=0)

		# Here the stacked data is recorded into a new fits file.
		outfile = 'stacked_fits/stacked_siriusbin2.fit'
		hdu = fits.PrimaryHDU(bin2)
		hdu.writeto(outfile, clobber=True)

		# Crops the image stored - Variable 'position' (defining
		# vertical position range) should be altered slightly to cater for each image as
		# they are in slightly different vertical positions. Then the croped  image is intigrated and fliped horizontally.
		position = (475, 505)
		size = (30, 8000)
		cutout = Cutout2D(bin2, position, size)
		bin2 = np.sum(cutout.data, axis=0)
		bin2 = np.flip(bin2, axis=0)

		return bin2

	def calib_sb2(self,bin2):

		# This function calibrates the spectra by converting the x-axis from pixel number to wavelength in nanometer.
		# Takes the data aray 'bin2' as an argument which must be passed on when the function is called

		############ Calibration numbers #############
		# Following arrays contains the data for Hydrogen Beta, Gamma, delta, and epsilon peaks of the spectrum
		h_beta = bin2[1200:1400]
		h_gamma = bin2[600:800]
		h_delta = bin2[350:450]
		h_epsilon = bin2[200:350]

		# Taken the minimum value (corresponding to the peak) between each range for each H value.
		# Because we took the average value between a range of values that didnt start at 0, we have to then
		# add the minimum of each range to the corresponding x value to place it within each range
		print ('Calibration data: ')

		print ('H-beta peak is: ', np.amin(h_beta))
		print ('H-beta pixel number is: ', np.argmin(h_beta) + 1200)
		print ('H-gamma peak is: ', np.amin(h_gamma))
		print ('H-gamma pixel number is: ', np.argmin(h_gamma) + 600)
		print ('H-delta peak is: ', np.amin(h_delta))
		print ('H-delta pixel number is: ', np.argmin(h_delta) + 350)
		print ('H-epsilon peak is: ', np.amin(h_epsilon))
		print ('H-epsilon pixel number is: ', np.argmin(h_epsilon) + 200)

		# The for loop runs from 1 to the length (or number of data points) + 1 in bin2
		# (+1 appears as the index of 'bin2' starts at 0). For every data point (0, 1,
		# 2 etc.) the loop performs the operation which is the equation of the line
		# y = mx + c for the pixel number plotted against the wavelength of the H lines.
		# Defines an empty list ready to store the calibrated values of
		# pixel/Wavelength

		pixels = [np.argmin(h_beta) + 1200, np.argmin(h_gamma) + 600, np.argmin(h_delta) + 350, np.argmin(h_epsilon) + 200]
		# known Balmer series wavelengths
		wavelengths = [486.1, 434.1, 410.2, 397.0]
		m, c, r_value, p_value, std_err = stats.linregress(pixels,wavelengths)
		print ('Gradient is: ', m)
		print ('Intercept is: ', c)
		s2_lamda = []
		for x in range(1,len(bin2)+1):
			l = (m*x) + c
			s2_lamda.append(l)

		return s2_lamda

	def plot_sb2(self,xvalues,image):
		# Plots the list of wavelength values 'calib' against the intensity values
		# Takes xvalues and image as argumentd which must be passed on when the function is called

		# Defines the spectrum figure and labels the axis
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Wavelength (nm)')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_title('Sirius absorption spectra')
		# Sets the x and y axis limites on the spectrum plot
		ax.set_xlim([min(xvalues), max(xvalues)])
		ax.set_ylim([min(image), max(image)])

		# Plots and shows the spectra
		plt.plot(xvalues, image, color='k')
		plt.savefig("spectra/siriusbin2.png")
		plt.show()

class siriusbin3():
	def initialise(self):
		image_list = [location+'\siriusbin3\siriusbin3_'+str(n)+'.fit' for n in range(440,462)]
		spectra = [fits.getdata(image) for image in image_list]
		bin3 = np.sum(spectra, axis=0)
		bin3 = np.flip(bin3, axis=0)

		outfile = 'stacked_fits/stacked_siriusbin3.fit'
		hdu = fits.PrimaryHDU(bin3)
		hdu.writeto(outfile, clobber=True)

		position = (348, 378)
		size = (30, 1677)
		cutout = Cutout2D(bin3, position, size)
		bin3 = np.sum(cutout.data, axis=0)
		bin3 = np.flip(bin3, axis=0)

		return bin3

	def calib_sb3(self,bin3):

		############ Calibration numbers #############
		h_beta = bin3[800:1000]
		h_gamma = bin3[400:600]
		h_delta = bin3[200:400]
		h_epsilon = bin3[150:200]

		print ('Calibration data: ')

		print ('H-beta peak is: ', np.amin(h_beta))
		print ('H-beta pixel number is: ', np.argmin(h_beta) + 800)
		print ('H-gamma peak is: ', np.amin(h_gamma))
		print ('H-gamma pixel number is: ', np.argmin(h_gamma) + 400)
		print ('H-delta peak is: ', np.amin(h_delta))
		print ('H-delta pixel number is: ', np.argmin(h_delta) + 200)
		print ('H-epsilon peak is: ', np.amin(h_epsilon))
		print ('H-epsilon pixel number is: ', np.argmin(h_epsilon) + 150)

		pixels = [np.argmin(h_beta) + 800, np.argmin(h_gamma) + 400, np.argmin(h_delta) + 200,   np.argmin(h_epsilon) + 150]
		wavelengths = [486.1, 434.1, 410.2, 397.0]
		m, c, r_value, p_value, std_err = stats.linregress(pixels,wavelengths)
		print ('Gradient is: ', m)
		print ('Intercept is: ', c)
		s3_lamda = []
		for x in range(1,len(bin3)+1):
			l = (m*x) + c
			s3_lamda.append(l)

		return s3_lamda

	def plot_sb3(self,xvalues,image):

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Wavelength (nm)')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_title('Sirius absorption spectra')
		ax.set_xlim([min(xvalues), max(xvalues)])
		ax.set_ylim([min(image), max(image)])

		plt.plot(xvalues, image, color='k')
		plt.savefig("spectra/siriusbin3.png")
		plt.show()

class jupiterbin2():
	def initialise(self):
		image_list = [location+'\jupiterbin2\jupiterhires_'+str(n)+'.fit' for n in range(463,539)]
		spectra = [fits.getdata(image) for image in image_list]
		bin2 = np.sum(spectra, axis=0)
		bin2 = np.flip(bin2, axis=0)

		outfile = 'stacked_fits/stacked_jupiterbin2.fit'
		hdu = fits.PrimaryHDU(bin2)
		hdu.writeto(outfile, clobber=True)

		bin2 = np.sum(bin2, axis=0)
		bin2 = np.flip(bin2, axis=0)

		return bin2

	def plot_jb2(self,image):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Pixel number')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_title('Jupiter absorption spectra')
		ax.set_xlim([0, len(image)+1])
		ax.set_ylim([min(image), max(image)])

		plt.plot(image, color='k')
		plt.savefig("spectra/jupiterbin2.png")
		plt.show()

class jupiterbin3():
	def initialise(self):
		image_list = [location+'\jupiterbin3\jupiterhires_'+str(n)+'.fit' for n in range(565,573)]
		spectra = [fits.getdata(image) for image in image_list]
		bin3 = np.sum(spectra, axis=0)
		bin3 = np.flip(bin3, axis=0)

		outfile = 'stacked_fits/stacked_jupiterbin3.fit'
		hdu = fits.PrimaryHDU(bin3)
		hdu.writeto(outfile, clobber=True)

		bin3 = np.sum(bin3, axis=0)
		bin3 = np.flip(bin3, axis=0)

		return bin3

	def plot_jb3(self,image):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Pixel number')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_title('Jupiter absorption spectra')
		ax.set_xlim([0, len(image)+1])
		ax.set_ylim([min(image), max(image)])

		plt.plot(image, color='k')
		plt.savefig("spectra/jupiterbin3.png")
		plt.show()

	def r_velocity(self):
		c = 299792458.0
		fitstop = fits.open(location+"jupiterhires_571stop.fit")
		im_top = np.array(fitstop[0].data)
		bin_top = np.sum(im_top, axis=0)
		bin_top = np.flip(bin_top, axis=0)
		bin_top = np.delete(bin_top, [0,1,2], 0)
		t = [bin_top[n] for n in range (0,1674)]
		t = [ t[n]/float(np.amax(t)) for n in range (0,1674)]
		peak_top = bin_top[1000:1250]
		top = np.argmin(peak_top) + 1000
		print ('Pixel number for the peak from top part of spectra: ', top)

		fitsbottom = fits.open(location+"jupiterhires_571sbottom.fit")
		im_bottom = np.array(fitsbottom[0].data)
		bin_bottom = np.sum(im_bottom, axis=0)
		bin_bottom = np.flip(bin_bottom, axis=0)
		bin_bottom = np.delete(bin_bottom, [0,1,2], 0)
		b = [bin_bottom[n] for n in range (0,1674)]
		b = [ b[n]/float(np.amax(b)) for n in range (0,1674)]
		peak_bottom = bin_bottom[1000:1250]
		bottom = np.argmin(peak_bottom) + 1000
		print ('Pixel number for the peak from bottom part of spectra: ', bottom)

		d = bottom - top
		rest = 727.0
		r_vel = (1.0/4.0)*( float(d*(0.1311) ) / float(rest) ) * c
		print ('Radial velocity of Jupiter is: ', r_vel*1e-3 , 'km/s')

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Pixels')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_title('Absorption spectra from both edges of Jupiter')
		ax.set_ylim([0, 1.2])
		plt.plot(t, color='k')
		plt.plot(b, color='b')
		ax.legend("TB", loc="upper left")
		plt.savefig('spectra/jupiteredges.png')
		plt.show()

class orionbin2():
	def initialise(self):
		image_list = [location+'\orionbin2\orionbin2_'+str(n)+'.fit' for n in range(414,417)]
		spectra = [fits.getdata(image) for image in image_list]
		bin2 = np.sum(spectra, axis=0)
		bin2 = np.flip(bin2, axis=0)

		outfile = 'stacked_fits/stacked_orionbin2.fit'
		hdu = fits.PrimaryHDU(bin2)
		hdu.writeto(outfile, clobber=True)

		bin2 = np.sum(bin2, axis=0)
		bin2 = np.flip(bin2, axis=0)

		return bin2

	def plot_ob2(self,image):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Wavelength (nm)')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_title('Orion emission spectra')

		o2_lamda = []
		for x in range(1,len(image)+1):
			l = (0.0877*x) + 374.3
			o2_lamda.append(l)
		ax.set_xlim([min(o2_lamda), max(o2_lamda)])
		ax.set_ylim([min(image), max(image)])
		plt.plot(o2_lamda, image, color='k')
		plt.savefig("spectra/orionbin2.png")
		plt.show()

class orionbin3():
	def initialise(self):
		image_list = [location+'\orionbin3\orionbin3_'+str(n)+'.fit' for n in range(417,421)]
		spectra = [fits.getdata(image) for image in image_list]
		bin3 = np.sum(spectra, axis=0)
		bin3 = np.flip(bin3, axis=0)

		outfile = 'stacked_fits/stacked_orionbin3.fit'
		hdu = fits.PrimaryHDU(bin3)
		hdu.writeto(outfile, clobber=True)

		bin3 = np.sum(bin3, axis=0)
		bin3 = np.flip(bin3, axis=0)

		return bin3

	def plot_ob3(self,image):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xlabel('Wavelength (nm)')
		ax.set_ylabel('Intensity (arbitrary unit)')
		ax.set_title('Orion emission spectra')

		o3_lamda = []
		for x in range(1,len(image)+1):
			l = (0.1309*x) + 374.1
			o3_lamda.append(l)
		ax.set_xlim([min(o3_lamda), max(o3_lamda)])
		ax.set_ylim([min(image), max(image)])
		plt.plot(o3_lamda, image, color='k')
		plt.savefig("spectra/orionbin3.png")
		plt.show()

fitsfunction()
