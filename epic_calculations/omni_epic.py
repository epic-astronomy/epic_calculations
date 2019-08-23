import numpy as np
import scipy.constants as const


def omniscope_cost(njs, pad=True, corr_cost=False):
    """Calculate computational cost for omniscope given the sizes of grid levels.

    Parameters
    ----------
    nfs : list, array of ints
        The sizes of the hierarchical grid levels. For example, the example given in Fig. 1
        of Tegmark 2010, njs=[5, 3, 3, 3, 3, 3]
    pad : bool, optional
        Whether to pad each grid level. Default is True.
    corr_cost : bool, optional
        Calculate correlation cost instead of FFT cost. Default is False.

    Returns
    -------
    cost : float
        The computational cost to form images (corr_cost==False) or correlations (corr_cost==True).

    """
    fft_factor = 5. / 2.
    njs = np.array(njs)
    nant = np.prod(njs)

    if corr_cost:
        nant = np.prod(njs)
        cost = nant * (nant - 1) / 2.
    else:
        if pad:
            cost = fft_factor * 2**njs.size * nant * (njs.size + np.log2(nant))
        else:
            cost = fft_factor * nant * np.log2(njs)
    return cost


class HierarchicalGrid():
    """Class to contain a definition of a hierarchical grid."""

    def __init__(self, basis_vectors, i_maxes, origin=[0., 0.]):
        """Initialize the class.

        Parameters
        ----------
        basis_vectors : list of lists/arrays
            List of basis vectors that define the hierarchical grid (a's in Tegmark 2010, eq 3).
            For example, Fig 1 of that paper would be [[1, 0], [8, 0], [40, 10], [0, 1], [2, 8], [-10, 40]]
        i_maxes : list of ints
            The number of grid points associated with each basis vector (max i in Tegmark eq 3).
            This class assumes i ranges from 0 to i_max - 1.
            For example, Fig 1 of that paper would be [5, 3, 3, 3, 3, 3]
        origin : list of floats, optional
            The origin from which to reference the basis vectors and i's above.
            Default is [0., 0.]

        """
        # TODO: Check inputs

        self.basis_vectors = basis_vectors
        self.i_maxes = i_maxes
        self.origin = np.array(origin)
        lengths = [np.linalg.norm(vec) for vec in basis_vectors]
        self.shortest_vec = np.min(lengths)

        self.calc_grid()

    def calc_grid(self):
        """Calculate the grid pixels for the basis vectors, i_maxes, and origin."""
        self.grid = []
        lvl = 0
        self._recursive_grid_level(lvl, self.origin)
        self.grid = np.array(self.grid)

    def _recursive_grid_level(self, lvl, loc):
        """Recursively walk through grid levels to add points to self.grid."""
        for i in range(self.i_maxes[lvl]):
            new_loc = loc + i * np.array(self.basis_vectors[lvl])
            if lvl == len(self.basis_vectors) - 1:
                self.grid.append(new_loc)
            else:
                self._recursive_grid_level(lvl + 1, new_loc)

    def contains(self, loc, mindist=None, radius=0):
        """Determine whether a location is contained within the grid.

        Parameters
        ----------
        loc : length-2 list, array, or tuple
            Location to test
        mindist : float, optional
            Minimum distance from a grid point to be consider contained.
            Defaults to length of shortest basis vector.
        radius : float, optional
            Radius of object to check containment. Default is 0 (only check specific location).

        Returns
        -------
        conained : bool
            True if loc is contained in grid, otherwise False

        """
        dists = np.linalg.norm(self.grid - np.array(loc).reshape(1, 2), axis=1)
        if mindist is None:
            mindist = self.shortest_vec
        if mindist <= np.min(dists):
            return False
        if radius > 0:
            # Check four corners
            for offset in np.array([[0, radius], [0, -radius], [radius, 0], [-radius, 0]]):
                dists = np.linalg.norm(self.grid - (np.array(loc) - offset).reshape(1, 2), axis=1)
                if mindist is None:
                    mindist = self.shortest_vec
                if mindist <= np.min(dists):
                    return False
        # If we made it this far we're good.
        return True


class array1D():
    """Class to hold 1D array info."""

    def __init__(self, ant_locs):
        """Initialize the class.

        Parameters
        ----------
        ant_locs : array_like of floats
            Antenna locations. Should be shape (Nant,)

        """
        self.ant_locs = np.array(ant_locs)
        self._calc_properties()

    def _calc_properties(self):
        """Calculate some properties based on the antenna locations."""
        self.bl_max = self.ant_locs.max() - self.ant_locs.min()
        self.bl_min = np.min(np.abs(np.diff(sorted(self.ant_locs))))

    def psf(self, return_angles=True):
        """Calculate the point spread function of the array.

        Parameters
        ----------
        return_angles : bool, optional
            Whether to also return the angle axis. Default is true.

        Returns
        -------
        angles : array of floats, optional
            Angle axis for the PSF. Returned if return_angles is True (default).
        psf : array of floats
            The point spread function.

        """
        bl_max = self.ant_locs.max() - self.ant_locs.min()
        bl_min = np.min(np.abs(np.diff(sorted(self.ant_locs))))
        theta_min = 1. / (2 * self.bl_max)
        theta_max = 1. / self.bl_min
        angles = np.arange(-theta_max, theta_max, theta_min)

        psf = np.zeros(len(angles))
        for i, theta in enumerate(angles):
            psf[i] = np.abs(np.sum(np.exp(2j * np.pi * theta * self.ant_locs)))
        psf = psf / np.max(psf)
        if return_angles:
            return angles, psf
        else:
            return psf

    def _recursive_Pmax(Pmaxes, costs, currPmaxes, Pmax_brute, Nv, brute_cost):
        for i in range(1, int(np.ceil(Pmax_brute / (2**(Nv - 1)))) + 1):
            if len(currPmaxes) + 1 == Nv:
                cost = omniscope_cost(np.array(currPmaxes + [i]) + 1, pad=True)
                if cost < brute_cost:
                    Pmaxes.append(currPmaxes + [i])
                    costs.append(cost)
                else:
                    break
            else:
                _recursive_Pmax(Pmaxes, costs, currPmaxes + [i], Pmax_brute, Nv, brute_cost)

    def get_candidate_pmaxes(self, amin=None, reorder=False):
        """Get potential Pmax's which are better than brute force.

        Parameters
        ----------
        amin : float, optional
            minimum grid spacing. Default is the minimum baseline length.
        reorder : bool, optional
            Reorder the outputs by cost. Default is False

        Returns
        -------
        Pmaxes : list of lists of ints
            The list of sets of Pmaxes which have a cost lower than the brute force method (1-level).
        costs : list of floats
            The relative costs associated with each Pmax set.

        """
        if amin is None:
            amin = self.bl_min
        Pmax_brute = int(np.ceil((self.bl_max - self.bl_min) / amin))
        brute_cost = omniscope_cost([Pmax_brute], pad=True)
        # This comes from argument about the padding penalty (see notebook)
        Nv_max = int(np.floor(np.log2(Pmax_brute + 1) / 2.))

        Pmaxes = [[Pmax_brute]]
        costs = [brute_cost]
        for Nv in range(2, Nv_max + 1):
            currPmaxes = []
            self._recursive_Pmax(Pmaxes, costs, currPmaxes, Pmax_brute, Nv, brute_cost)

        if reorder:
            inds = np.argsort(costs)
            Pmaxes = list(np.array(Pmaxes)[inds])
            costs = list(np.array(costs)[inds])

        return Pmaxes, costs


class HierarchicalGrid_1D():
    """Class to contain a definition of a hierarchical grid, in 1D."""

    def __init__(self, basis_vectors, Pmaxes, origin=0.):
        """Initialize the class.

        Parameters
        ----------
        basis_vectors : list of floats
            List of basis vectors (lengths) that define the hierarchical grid
            (a's in Tegmark 2010, eq 3).
        Pmaxes : list of ints
            The largest vector coefficient for basis vectors. Because the
            coefficients start at 0, these are i_max - 1 in Tegmark eq 3.
        origin : float, optional
            The origin from which to reference the basis vectors and i's above.
            Default is 0.

        """
        # TODO: Check inputs

        self.basis_vectors = basis_vectors
        self.Pmaxes = Pmaxes
        self.origin = origin
        self.shortest_vec = np.min(self.basis_vectors)

        self.calc_grid()

    def calc_grid(self):
        """Calculate the grid pixels for the basis vectors, i_maxes, and origin."""
        self.grid = []
        lvl = 0
        self._recursive_grid_level(lvl, self.origin)
        self.grid = np.array(self.grid)

    def _recursive_grid_level(self, lvl, loc):
        """Recursively walk through grid levels to add points to self.grid."""
        for i in range(self.Pmaxes[lvl] + 1):
            new_loc = loc + i * self.basis_vectors[lvl]
            if lvl == len(self.basis_vectors) - 1:
                self.grid.append(new_loc)
            else:
                self._recursive_grid_level(lvl + 1, new_loc)

    def contains(self, loc, mindist=None, radius=0):
        """Determine whether a location is contained within the grid.

        Parameters
        ----------
        loc : length-2 list, array, or tuple
            Location to test
        mindist : float, optional
            Minimum distance from a grid point to be consider contained.
            Defaults to length of shortest basis vector / 2.
        radius : float, optional
            Radius of object to check containment. Default is 0 (only check specific location).

        Returns
        -------
        conained : bool
            True if loc is contained in grid, otherwise False

        """
        dists = np.abs(self.grid - loc)
        if mindist is None:
            mindist = self.shortest_vec / 2.
        if mindist <= np.min(dists):
            return False
        if radius > 0:
            # Check plus or minus
            for offset in [-radius, radius]:
                dists = np.abs(self.grid - (loc - offset))
                if mindist <= np.min(dists):
                    return False
        # If we made it this far we're good.
        return True

    def contains_array(self, array, mindist=None, radius=0, return_nant=False):
        """Determine whether an entire array is contained within the grid.

        Parameters
        ----------
        array : array1D
            Object containing the array in question
        mindist : float, optional
            Minimum distance from a grid point to be consider contained.
            Defaults to length of shortest basis vector / 2.
        radius : float, optional
            Radius of object to check containment. Default is 0 (only check specific location).
        return_nant : bool, optional
            Return number of antennas contained in the grid. Default is False.

        Returns
        -------
        contained : bool or int
            If return_ant is False: True if array is completely contained in grid, otherwise False.
            If return_ant is True: Return the number of antennas contained.

        """
        nant = 0
        all_contained = True
        for loc in array.ant_locs:
            contained = self.contains(loc, mindist=mindist, radius=radius)
            nant += contained
            if not contained and not return_nant:
                return False
        if return_nant:
            return nant
        else:
            return True

    def grid_vs_array_score(self, array, mindist=None, radius=0):
        """Score for the grid based on distances from antennas that are not contained.

        Parameters
        ----------
        array : array1D
            Object containing the array in question
        mindist : float, optional
            Minimum distance from a grid point to be consider contained.
            Defaults to length of shortest basis vector / 2.
        radius : float, optional
            Radius of object to check containment. Default is 0 (only check specific location).

        Returns
        -------
        score : float
            The sum of distances of each non-contained antenna to its nearest grid point.

        """
        score = 0.
        for loc in array.ant_locs:
            contained = self.contains(loc, mindist=mindist, radius=radius)
            if not contained:
                score += np.min(np.abs(self.grid - loc))
        return score


def _rand2physical(params, omax, amin, amax):
    """Quick function to translate random numbers to physical parameters.

    Parameters
    ----------
    params : array of floats
        Numbers in range [0, 1) representing the origin and basis vectors
    omax : float
        Maximum origin value
    amin : float
        Minimum basis vector length. This will be the first element.
    amax : float
        Maximum allowed basis vector length.

    Returns
    -------
    origin : float
        Origin now in physical parameter space.
    basis_vecs : array_like of floats
        Basis vectors now in physical parameter space.

    """
    basis_vecs = np.zeros(params.size - 1)
    basis_vecs[0] = amin
    for i in range(len(basis_vecs[1:])):
        amax_temp = amax / 2**(len(basis_vecs) - i - 2)
        basis_vecs[i + 1] = 2 * basis_vecs[i] + params[i + 2] * (amax_temp - 2 * basis_vecs[i])
    origin = omax - params[0] * basis_vecs[-1]

    return origin, basis_vecs


def func_to_min(params, array, Pmaxes, amin=None, mindist=None, radius=0):
    """Given list of Pmaxes and an array, try a bunch of origins and a's.

    Parameters
    ----------
    params : array of floats
        These are the parameters to be searched. This function is defined such that
        these parameters can range [0, 1], which is specified in the call to
        scipy.optimize.minimize.
    array : array1D
        Object containing the array to grid.
    Pmaxes : list of ints
        The largest vector coefficient for basis vectors.
        Assumed to be in order from smallest basis vector to largest.
    amin : float, optional
        Minimum grid spacing. Defaults to minimum baseline.
    mindist : float, optional
        Minimum distance from a grid point to be consider contained.
        Defaults to length of shortest basis vector / 2.
    radius : float, optional
        Radius of antennas to check containment. Default is 0 (only check specific location).

    Returns
    -------
    score : float
        The score associated with the grid specified by params and Pmaxes against
        the array.

    """
    if amin is None:
        amin = array.bl_min
    amax = array.bl_max / Pmaxes[-1]  # TODO: make smarter
    omax = np.min(array.ant_locs) - radius

    origin, basis_vecs = _rand2physical(params, omax, amin, amax)
    this_grid = HierarchicalGrid_1D(basis_vecs, Pmaxes, origin=origin)
    score = this_grid.grid_vs_array_score(array, mindist=mindist, radius=radius)

    return score
