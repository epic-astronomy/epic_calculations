import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt


class TelescopeObservation():
    """ A class which defines a telescope and observations parameters
    (e.g. frequency, bandwidth, integration, etc)
    """

    def __init__(self, layout=None, Nant=None, Darray=None, Dant=None, grid_size=None, f0=None,
                 bandwidth=None, Nchan=None, integration=None, in_bit_depth=4, out_bit_depth=32):
        """ Initialize the class

        Parameters
        ----------
        layout : array of floats, optional
            (Nant, 2 or 3) Array of antenna locations in meters. If Nant and Darray are not set,
            layout will be used to determine them.
        Nant : int
            Number of antennas in the array. If layout is also set, Nant will be set to
            min(Nant, layout.shape[0]). If layout has more antennas than Nant, the first Nant will be used.
        Darray : float
            Diameter of the array in meters. If layout is also set, Darray will be set to
            to either this argument or the Darray determined from the layout, whichever is smaller.
        Dant : float
            Diameter of the antennas in meters.
        grid_size : float
            Size of grid elements in wavelengths. If None (default), use Dant * f0 / speed_of_light.
        f0 : float or int
            Center frequency in MHz.
        bandwidth : float
            Bandwidth in MHz
        Nchan : int
            Number of frequency channels.
        integration : float
            Integration time for output data, in seconds.
        in_bit_depth : int
            Number of bits per sample of input. Default is 4 (x2 for complex).
        out_bit_depth : int
            Number of bits per sample output. Default is 32.
        """
        self.layout = layout
        if self.layout is not None:
            # center the array about zero
            self.layout -= np.mean(self.layout, axis=0).reshape(1, -1)
            # find diameter of array
            rs = np.sqrt(np.sum(np.abs(self.layout)**2, axis=1))
            max_u = 2 * np.max(rs)
            if Darray is None:
                self.Darray = max_u
            else:
                self.Darray = np.min([max_u, Darray])
            # Remove antennas outside radius
            if self.Darray < max_u:
                inds = np.where(2 * rs > self.Darray)[0]
                self.layout = np.delete(self.layout, inds, axis=0)
            if Nant is None:
                self.Nant = self.layout.shape[0]
            else:
                self.Nant = np.min([Nant, self.layout.shape[0]])
                self.layout = self.layout[:self.Nant]
        else:
            self.Nant = Nant
            self.Darray = Darray

        self.Dant = Dant
        if grid_size is None:
            grid_size = self.Dant * const.speed_of_light / (f0 * 1e6)
        self.grid_size = grid_size
        self.f0 = f0
        self.bandwidth = bandwidth
        self.Nchan = Nchan
        self.integration = integration
        self.in_bit_depth = in_bit_depth
        self.out_bit_depth = out_bit_depth

        self._set_dependents()

    def _set_dependents(self):
        # Some calculated values that are useful
        complex_factor = 2.
        bits_per_byte = 8.

        self.in_bw = (self.Nant * self.bandwidth * self.in_bit_depth
                      * complex_factor / bits_per_byte)  # MBps
        self.lambda0 = const.speed_of_light / (self.f0 * 1e-6)
        self.channels = np.linspace(self.f0 - self.bandwidth / 2,
                                    self.f0 + self.bandwidth / 2, num=self.Nchan)  # MHz
        self.lambdas = const.speed_of_light / (self.channels * 1e6)
        self.cadence = 1. / (self.bandwidth / self.Nchan * 1e6)  # Time per FFT in seconds

    def vanilla_EPIC_stats(self, padding=2, verbose=True):
        """ Calculate the computation requirement for Vanilla EPIC.

        Parameters
        ----------
        padding : float or int, optional
            Factor to pad grid. Default is 2, which will make the same size
            grid as no padding FX.
        verbose : bool
            Option to print more stats than just what is returned.

        Returns
        -------
        total_flops : float
            Number of floating point operations per second required to process the data.
        """
        max_u = self.Darray * (self.f0 + self.bandwidth / 2.) * 1e6 / const.speed_of_light
        npix = padding * 2 ** (np.ceil(np.log2(max_u / self.grid_size)))  # pixels per side

        gridding_flops_per_chan = self.Nant * (self.Dant / self.lambdas / self.grid_size)**2 / self.cadence

        fft_flops_per_chan = 5 * npix * np.log2(npix) / self.cadence

        # TODO: add squaring, averaging
        # TODO: output bandwidth
        # TODO: report FoV and resolution

        total_flops = np.sum(gridding_flops_per_chan) + fft_flops_per_chan * self.Nchan

        if verbose:
            print('Input bandwidth = ' + str(self.in_bw) + ' MBps')
            print('Npix = ' + str(npix))
            print('Cadence = ' + str(self.cadence) + ' s')
            print('Gridding GFLOPS per channel (avg) = ' + str(np.mean(gridding_flops_per_chan * 1e-9)))
            print('Total gridding GFLOPS = ' + str(np.sum(gridding_flops_per_chan * 1e-9)))
            print('FFT GFLOPS per channel = ' + str(fft_flops_per_chan * 1e-9))
            print('Total FFT GFLOPS = ' + str(fft_flops_per_chan * self.Nchan * 1e-9))
            print('All the GFLOPS = ' + str(1e-9 * total_flops))

        return total_flops

    def FX_stats(self, padding=1, verbose=True):
        """ Calculate the computation requirement for FX correlator + image.

        Parameters
        ----------
        padding : float or int, optional
            Factor to pad grid. Default is 1 (no padding).
        verbose : bool
            Option to print more stats than just what is returned.

        Returns
        -------
        total_flops : float
            Number of floating point operations per second required to process the data.
        """

        nbls = self.Nant * (self.Nant - 1.) / 2.
        corr_flops_per_channel = nbls / self.cadence

        max_u = self.Darray * (self.f0 + self.bandwidth / 2.) * 1e6 / const.speed_of_light
        npix = 2 * padding * 2 ** (np.ceil(np.log2(max_u / self.grid_size)))  # pixels per side

        gridding_flops_per_chan = self.Nant * (2 * self.Dant / self.lambdas / self.grid_size)**2 / self.integration
        fft_flops_per_chan = 5 * npix * np.log2(npix) / self.integration

        total_flops = self.Nchan * (corr_flops_per_channel + gridding_flops_per_chan.mean() + fft_flops_per_chan)

        if verbose:
            print('Input bandwidth = ' + str(self.in_bw) + ' MBps')
            print('Npix = ' + str(npix))
            print('Cadence = ' + str(self.cadence) + ' s')
            print('Corr GFLOPS per channel = ' + str(corr_flops_per_channel * 1e-9))
            print('Gridding GFLOPS per channel (avg) = ' + str(np.mean(gridding_flops_per_chan * 1e-9)))
            print('Total gridding GFLOPS = ' + str(np.sum(gridding_flops_per_chan * 1e-9)))
            print('FFT GFLOPS per channel = ' + str(fft_flops_per_chan * 1e-9))
            print('Total FFT GFLOPS = ' + str(fft_flops_per_chan * self.Nchan * 1e-9))
            print('All the GFLOPS = ' + str(1e-9 * total_flops))

        return total_flops

    def plot_layout(self):
        plt.plot(self.layout[:, 0], self.layout[:, 1], 'o')
        plt.axes().set_aspect('equal', 'datalim')
