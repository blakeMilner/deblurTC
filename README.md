kisip
=====

Image deblurring using a series of images (from dozens to thousands). This method is highly effective for series of many low exposure images that contiain blurring or spatial shifting of an object of interest. Written in C and utilizing MPI for parallelization, this scientific computing package is suitable for a computing cluster.

This package was originally developed and used for speckle interferometry by astrophysics studying the sun: http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1336991

Originally coded by Friedrich Wöger and Oskar von der Lühe II.

The portions of this code pertaining to astrophysics (such as telescope diffraction, modeling of atmospheric aberrations, etc.) has been stripped out to allow for general use. 

The core of this method is transforming all frames into their Bispectrum (a type of 4-D power spectrum), averaging them together, and producing a single, deblurred image that is a compositive of all frames. The 4-D Bispectrum is the Fourier Transform of the Triple Correlation (related to regular correlation used to register images) and can help reconstruct the phase of the original object, thus undoing any effects due to translation or blurring of the image. The method has been well-established since the 60's in the field of astrophysics:
http://adsabs.harvard.edu/abs/1993A%26A...278..328H
