#!/usr/bin/env python
import os
import tqdm

import rasterio
import numpy as np
import pandas as pd
import spectral
from typing import List, Tuple, Union, Optional


sample_axis = 0
band_axis = 1


def bandwind_from_wls(centers: List[float], wl1: float, wl2: float) -> List[int]:
    """Given a list of band centers and a pair of floats
    defining a range of wavelengths, find the bands that are covered by this
    range (but not outside). Bands are 0-based, so first band is 0"""
    ##First check if there's overlap, return None if not
    if (wl1 > centers[-1]) or (wl2 < centers[0]):
        return None
    minband = len(centers) + 1
    maxband = -1
    for ci, c in enumerate(centers):
        if (c > wl1) and (c < wl2):
            if ci < minband:
                minband = ci
            if ci > maxband:
                maxband = ci
    if (maxband < 0) or (minband > len(centers)):
        ##Should not get here if there's overlap, unless wl range is tiny
        return None
    return [minband, maxband]


class AutoMCU(object):
    """Class to setup and run AutoMCU according to specify configuration
    Inputs:
    em_csvs           : List of csv files that contain endmember spectra, one
                        for each class - typically PV, NPV, and Bare. CSVs
                        should be arranged such that samples are columns, and
                        bands are rows, with the first column containing the
                        band center wavelengths in nanometers. At least two
                        classes must be specified.
    em_counts=[]      : Either a list of integers specifying how many em
                        samples to collect from each em class, a single
                        integer specifying how many samples to take from all
                        classes, or an empty list (or equivalently None). In
                        the first case, the length of this list must match
                        the length of em_csvs. In the first and second case,
                        the integer 0 means all available samples. An empty
                        list or None is equivalent to [0]*len(em_csvs).
    wl_ranges=None    : A list or tuple of 2-element lists or tuples of
                        floats, i.e. [(600.0,750),(2030,2380)], that will be
                        used to identify the band ranges used for unmixing.
                        Data from the bands specified by these wavelengths
                        will be lumped into a single array for fitting.
                        A minimum of one band window must be specified.
    band_ranges=None  : Similar to band_ranges, but specifying starting and
                        ending integer band indexes (first band is band 0)
                        for whatever image this is applied on.
    tied=True         : Reflectance spectra from each band range are
                        subtracted from the first value in the window to
                        help relieve brightness issues. (First value in
                        unmixing array becomes 0.0 for each windows start)
    divide=False      : Similar to tied, but window values are divided by
                        the first entry. (First value in each window will
                        become 1.0)
    shade=False       : Add an endmember for shade, which is constant 0.0
    sum_to_one=True   : Add an element with value 1.0 to the end of the
                        unmixing array and the input refl to help coerce sums
                        of the endmember fractions to be near 1.0.
    iterations=30     : How many MC runs to do for each pixel.
    outscale=1000     : Scale value applied to unmixing results to make them
                        fit into data type specified for output.
    nodata=-9999      : Output no data value written into metadata. If None,
                        missing data in output will be 0, but this will not
                        be recorded in metadata.
    dtype="int16"     : Rasterio-usable specification of output map data type.
    otype="GTiff"     : GDAL shortname for output data format.
    co=[]             : GDAL data format writing options.
    emfwhm=False      : Assume second column of em csvs are fwhm values
    nointerp=False    : Assume em csvs have same count and wavelengths as
                        input reflectance data (will be checked)
    match_c           : Recreate a bug in original C code that shifted input
                        reflectance data by one index/band
    """

    def __init__(
        self,
        em_csvs: List[str],
        em_counts: List[int] = [],
        em_names: List[str] = [],
        band_ranges: List[Tuple[float, float]] = None,
        wl_ranges: List[Tuple[float, float]] = None,
        tied: bool = True,
        divide: bool = False,
        shade: bool = False,
        soil_frac: bool = False,
        sum_to_one: bool = True,
        iterations: int = 30,
        outscale: float = 1000,
        nodata: float = -9999,
        dtype: str = "int16",
        otype: str = "GTiff",
        co: List[str] = [],
        emfwhm: bool = False,
        nointerp: bool = False,
        match_c: bool = False,
    ):
        ##Check inputs

        self.tied: Optional[int] = tied
        self.divide: Optional[int] = divide
        self.shade: Optional[int] = shade
        self.soilfrac: Optional[float] = soil_frac
        self.sum1: Optional[bool] = sum_to_one
        self.niter: Optional[int] = iterations  # defined by user
        self.outscale: Optional[int] = outscale
        self.nodata: Optional[float] = nodata
        self.dtype: Optional[str] = dtype
        self.otype: Optional[str] = otype
        self.dointerp: Optional[bool] = not nointerp
        self.match_c: Optional[str] = match_c
        self.dopt: Optional[dict] = dict([s.split("=") for s in co])

        ##Filled later
        self._ems_interp = None
        self._ems_proc = None
        self._em_avg = None

        ##Check the endmember args
        ##########################
        # em_csv should be a list
        assert isinstance(em_csvs, list)
        # For now, lets say em classes should be 2 or more
        assert len(em_csvs) >= 2
        # Counts should either be an empty, a single integer, or a list of
        # integers the same size as em_csvs
        if not em_counts:
            em_counts = [0] * len(em_csvs)
        elif isinstance(em_counts, list):
            assert len(em_csvs) == len(em_counts)
        else:
            em_counts = [em_counts] * len(em_csvs)
        # Check that files exist
        for csv in em_csvs:
            assert os.path.exists(csv)
        self._ems_raw = []
        self._ems_wl = []
        self._ems_fwhm = []
        for csv, count in zip(em_csvs, em_counts):
            ##Here csv should be formatted with a header row and samples listed
            ## in columns, bands as rows. Wl should be the first columns, and
            ## if --emfwhm is specified, fwhm should be column 2
            tmpdf = pd.read_csv(csv)
            self._ems_wl.append(tmpdf.iloc[:, 0].to_numpy())
            st = 1
            if emfwhm:
                self._ems_fwhm.append(tmpdf.iloc[:, 1].to_numpy())
                st = 2
            else:
                self._ems_fwhm.append(None)
            if count > 0:
                maxcount = tmpdf.shape[1] - st
                if count > maxcount:
                    count = maxcount
                end = st + count
                print(f"Reading {count} samples from {csv}")
                self._ems_raw.append(tmpdf.iloc[:, st:end].to_numpy().T)
            else:
                print(f"Reading all samples of {csv}")
                self._ems_raw.append(tmpdf.iloc[:, st:].to_numpy().T)

        # Save the number of em classes
        self.num_class = len(self._ems_raw)

        ##Check for em_names
        if em_names:
            self._em_names = em_names
        else:
            self._em_names = [f"Class{n:d}" for n in range(1, self.num_class + 1)]

        if self.shade:
            self.num_class += 1
            self._em_names.append("Shade")

        ##Check the specified ranges
        ############################
        # One and only one should be not None
        if band_ranges is not None:
            # There should be at least one window
            self._wl_ranges = None
            assert isinstance(band_ranges, list) or isinstance(band_ranges, tuple)
            assert len(band_ranges) > 0
            # If the first element is a float, then assume only one range
            if not (
                isinstance(band_ranges[0], list) or isinstance(band_ranges[0], tuple)
            ):
                self._band_ranges = [band_ranges[:2]]
            else:
                self._band_ranges = [list(p[:2]) for p in band_ranges]
            # Test that all range entries can be floats
            for pr in self._band_ranges:
                for it in pr:
                    try:
                        int(it)
                    except ValueError:
                        raise RuntimeException(
                            f"Window range entry {it} could not be converted to int"
                        )
                    except TypeError:
                        raise RuntimeException(f"Invalid window range entry: {it}")
            pass
        elif wl_ranges is not None:
            # There should be at least one window
            self._band_ranges = None
            assert isinstance(wl_ranges, list) or isinstance(wl_ranges, tuple)
            assert len(wl_ranges) > 0
            # If the first element is a float, then assume only one range
            if not (isinstance(wl_ranges[0], list) or isinstance(wl_ranges[0], tuple)):
                self._wl_ranges = [wl_ranges[:2]]
            else:
                self._wl_ranges = [list(p[:2]) for p in wl_ranges]
            # Test that all range entries can be floats
            for pr in self._wl_ranges:
                for it in pr:
                    try:
                        float(it)
                    except ValueError:
                        raise RuntimeException(
                            f"Window range entry {it} could not be converted to float"
                        )
                    except TypeError:
                        raise RuntimeException(f"Invalid window range entry: {it}")
        else:
            raise ArgumentError("One of wl_ranges or band_ranges must be specified")
        return

    def _find_hdr_name(self, filepath: str) -> str:
        """Find the name of the hdr file associated with the input img file

        Args:
            filepath (str): The filepath of the input img file

        Returns:
            str: The filepath of the hdr file
        """
        inhdr = None
        ##Go through common hdr name practices looking for one that exists
        inhdropts = [os.path.splitext(filepath)[0] + ".hdr", filepath + ".hdr"]
        for opt in inhdropts:
            if os.path.exists(opt):
                print(f"Found input hdr file at {opt}")
                inhdr = opt
                break
        ##Will be None if no options exist
        return inhdr

    def interpolate_ems(self, wl, fwhm=None):
        """
        Interpolate end member spectra to the specified wavelengths.

        Parameters
        ----------
        wl : array-like
            Wavelengths to interpolate to. Should be in the same units as
            self._ems_wl.
        fwhm : array-like or None
            FWHM to use for the interpolation. If None, the average spacing of
            `wl` will be used.

        Returns
        -------
        None

        """
        if fwhm is None:
            avg_spacing = np.diff(self.wl).mean()
            fwhm = [avg_spacing] * len(wl)
        self._ems_interp = []
        for rawmat, rawwl, rawfwhm in zip(self._ems_raw, self._ems_wl, self._ems_fwhm):
            ##Build fwhm if needed
            if rawfwhm is None:
                avg_spacing = np.diff(rawwl).mean()
                rawfwhm = [avg_spacing] * len(rawwl)
            ##Interpolate
            brs = spectral.BandResampler(rawwl, wl, rawfwhm, fwhm)
            self._ems_interp.append(
                np.concatenate([brs(v)[np.newaxis,:] for v in rawmat], \
                               axis=sample_axis))
            ##Check that interpolation did not result in nans in region bands
            for bandlist in self.regdefs:
                bands_with_nan = \
                    np.array(bandlist)[np.any(~np.isfinite(
                                    self._ems_interp[-1][:,bandlist],
                                    ),axis=sample_axis)]
                if len(bands_with_nan) > 0:
                    raise RuntimeError("Interpolator produced nans for bands: "+\
                                       f"{bands_with_nan} - check overlap")
        return

    def average_interpolated_ems(self):
        if not self._ems_interp:
            raise RuntimeError("interpolated_ems() has not been called yet")
        use_bands = []
        for rb in self.regdefs:
            use_bands.extend(rb)
        self._em_avg = np.concatenate(
            [
                emarr[:, use_bands].mean(axis=sample_axis, keepdims=True)
                for emarr in self._ems_interp
            ],
            axis=sample_axis,
        )
        return

    def make_unmixing_array(self, arr, match_c=False):
        """Select correct bands from arr, a NxB array, where n is number of
        samples and B is number of bands, and perform tying or dividing if
        requested"""
        dat_list = [arr[:, bandlist] for bandlist in self.regdefs]

        if self.divide:
            for dat in dat_list:
                dat /= dat[:, 0][:, np.newaxis]
        elif self.tied:
            for dat in dat_list:
                dat -= dat[:, 0][:, np.newaxis]

        dat = np.concatenate(dat_list, axis=band_axis)
        ##Bug in automcu.c makes Y(a.k.a., b) data shifted by one index
        ## This results in first band being omitted, and lastband being 0, since
        ##  fsubimgslice array is 0 after last svd band:
        ##   /* This is a bug since fsubimgslice is 0-based indexing, and it  */
        ##   /*     should be fsubimgslice[i][k-1] */
        ##   for (k=1; k<=num_prim_imgbands; k++) imgspec[k] = fsubimgslice[i][k];
        ##
        ##   if (sum_to_one == TRUE) imgspec[num_prim_svdbands] = 1.0;
        ## Re-enact this bug if requested
        if match_c:
            dat = np.roll(dat, -1, axis=band_axis)
            dat[:, -1] = 0.0

        if self.sum1:
            return np.concatenate(
                [dat, np.full((dat.shape[0], 1), 1.0)], axis=band_axis
            )
        else:
            return dat

    def process_iterpolated_ems(self):
        """With em data now matching application bands, process for unmixing"""
        if self._ems_interp is None:
            raise RuntimeError("interpolate_ems() has not been called yet")

        self._ems_proc = []
        for dat_i, dat in enumerate(self._ems_interp):
            procdat = self.make_unmixing_array(dat)
            ##BandResampler returns NaN for wl with no overlap
            if np.any(~np.isfinite(procdat)):
                raise RuntimeError("Processed em data contains non-finite values")
            self._ems_proc.append(procdat)
        return

    def _random_unmixing_libs(self, n):
        """
        This function is used to pick n random samples from the EMs_proc list
        and stack them together to form a tensor of size n x num_class x num_bands.
        """
        picks = np.concatenate(
            [
                s[np.random.choice(s.shape[0], size=n, replace=True), :][
                    :, np.newaxis, :
                ]
                for s in self._ems_proc
            ],
            axis=band_axis,
        )
        return picks  # n x num_class x num_bands

    def unmix_array(self, arr):
        """Actually do the unmixing for an array of valid refl data of shape
        N x B, where N is number of samples and B is number of raw spectral
        bands"""

        ##Will use pseudoinverse to linearly solve Y = X*coef
        ##Done for all samples simultaneously thanks to numpy magic

        ##Check that some members are not None
        if self._ems_proc is None:
            raise RuntimeError("process_interpolated_ems() has not been called yet")

        if self._em_avg is None:
            raise RuntimeError("average_interpolated_ems() has not been called yet")

        #########################
        ##Prepare Y array (N x B)
        #########################

        ##Process data in to unmixing array by applying tieing/dividing/sum1
        Y = self.make_unmixing_array(arr, self.match_c)
        num_samples = Y.shape[sample_axis]
        num_bands = Y.shape[band_axis]
        if num_samples < 1:
            return None
            

        # Y is shape (num_samples, num_bands)

        ##Also get original refl spectrum for the selected bands
        ##Shape (num_samples, num_bands)
        Yorig = np.concatenate(
            [arr[:, bandlist] for bandlist in self.regdefs], axis=band_axis
        )

        ######################################################
        ##Prepare X matrices from randomly selected endmembers
        ######################################################

        ##Collect needed number of random draws from SPECs
        ## One randomly-selected sample from each class for each iteration for
        ##    each sample
        ##shape (num_samples, self.niter, num_class, num_bands)
        ######################################################################
        random_selections = self._random_unmixing_libs(
            self.niter * num_samples
        ).reshape(num_samples, self.niter, self.num_class, num_bands)

        #################
        ##Do the unmixing
        #################

        ##Collect pseudo_inverses into a new array
        # shape (num_samples, self.niter, num_bands, num_class)
        ####################################################
        pinvs = np.linalg.pinv(random_selections)

        ##Get optimal coefficients for each sample and iteration
        ########################################################
        # shape num_samples, self.niter, num_class)
        ########################################################
        coefs = np.concatenate(
            [
                np.dot(pinvs[c, i, :].T, Y[c, :])[np.newaxis, :]
                for c in range(num_samples)
                for i in range(self.niter)
            ],
            axis=sample_axis,
        )
        coefs = coefs.reshape(num_samples, self.niter, self.num_class)

        #######################
        ##Compute trimmed stats
        #######################

        ##Compute quantiles across just iterations
        ## Each is shape (num_samples, num_class)
        ##########################################
        q10, q90 = np.percentile(coefs, [10, 90], axis=1)

        ##Trim 10% of coefs from each end by sample and class
        #####################################################
        coefs[
            np.logical_or(coefs < q10[:, np.newaxis, :], coefs > q90[:, np.newaxis, :])
        ] = np.nan

        ##Compute mean and stdev of middle 80%
        ## each shape (num_samples, num_class)
        ######################################
        means = np.nanmean(coefs, axis=band_axis)  # num_samples x num_class
        sds = np.nanstd(coefs, axis=band_axis)  # num_samples x num_class

        ######################################################################
        ##Remix expected spectra for each sample and compute sample-level RMSE
        ######################################################################

        ##Build mixed spectra (in original refl scale for selected bands)
        ##(num_samples, num_class) ⋅ (num_class, num_bands) =
        ##  (num_samples, num_bands)
        ##  Does not include band of 1s from sum1
        #################################################################
        preds = np.dot(means, self._em_avg)

        ##RMSE
        ## Shape (num_samples)
        ######
        rmses = np.sqrt(
            np.sum((preds - Yorig) ** 2, axis=band_axis, keepdims=True) / num_bands
        )

        return np.concatenate([means, sds, rmses], axis=band_axis)

    def apply_image(self, refl_path, input_scale, output_path, hdrpath=None, nblocks=0):

        ##Get info from input image
        ###########################

        ## Grab shape and metadata
        with rasterio.open(refl_path) as inref:
            imgshape = inref.shape

            imgoffset = (0, 0)
            outmeta = inref.meta.copy()
            input_blocks = list(inref.block_windows())

        ## Modify if we're not using all blocks
        if len(nblocks) > 0:
            nblocks = [int(x) for x in nblocks[0].split(",")]  ##Ealhe Change
            input_blocks = input_blocks[nblocks[0] : nblocks[1]]  ##Elahe change
            print(f"Unmixing image rows within range of: {nblocks}")  ##Elahe change
        minr, minc = imgshape
        maxr, maxc = (0, 0)
        for w_i, wind in input_blocks:
            minr = min(minr, wind.row_off)
            minc = min(minc, wind.col_off)
            maxr = max(maxr, wind.row_off + wind.height - 1)
            maxc = max(maxc, wind.col_off + wind.width - 1)
        imgshape = (maxr - minr + 1, maxc - minc + 1)
        imgoffset = (minr, minc)

        ##Try to find hdr if not specified
        if not hdrpath:
            hdrpath = self._find_hdr_name(refl_path)
        if not hdrpath:
            raise RuntimeError("No matching hdr found")
        if not os.path.exists(hdrpath):
            raise RuntimeError(f"HDR {hdrpath} does not exist")
        input_hdr = spectral.envi.read_envi_header(hdrpath)

        ##Get wavelengths from header
        if "wavelength" in input_hdr:
            self.wl = [float(wl) for wl in input_hdr["wavelength"]]
        else:
            raise RuntimeError(f"Could not find wavelengths in hdr")
        ##Get fwhm from header or build
        if "fwhm" in input_hdr:
            self.fwhm = [float(hm) for hm in input_hdr["fwhm"]]
        else:
            ## If fwhm missing, just assume equivalent to center spacing
            avg_spacing = np.diff(self.wl).mean()
            self.fwhm = [avg_spacing] * len(self.wl)

        ##Check that self.wl and self.fwhm lengths match up
        if len(self.wl) != len(self.fwhm):
            raise RuntimeError(
                f"Number of hdr wavelengths {len(self.wl)} != "
                # "number of fwhm {len(self.fwhm)}" ##Elahe comment this line
            )

        ##Check that number of self.wl is the same as the number of bands
        num_image_bands = int(input_hdr["bands"])
        if len(self.wl) != int(input_hdr["bands"]):
            raise RuntimeError(
                f"Number of hdr wavelengths {len(self.wl)} != "
                "number of image bands {num_image_bands}"
            )

        ##Figure out band ranges
        ########################

        ## Build _band_ranges if wl_ranges provided
        if self._wl_ranges is not None:
            self._band_ranges = []
            for wlr in self._wl_ranges:
                br = bandwind_from_wls(self.wl, *wlr[:2])
                if br is None:
                    raise RuntimeError(
                        "Could not get band indices from wl" f" range {wlr}"
                    )
                self._band_ranges.append(br)
                print(
                    f"Computed band range {br}"
                    f" ({self.wl[br[0]]},{self.wl[br[1]]})"
                    f" from wavelength window {wlr}"
                )

        ##Convert _band_ranges to lists of indices
        self.regdefs = []
        print(f"Found {len(self._band_ranges)} band regions:")
        for br_i, br in enumerate(self._band_ranges):
            self.regdefs.append(list(range(br[0], br[1] + 1)))
            print(f"{self.regdefs[-1]}")

        num_selected_bands = sum([len(rd) for rd in self.regdefs])
        num_unmixing_bands = num_selected_bands
        if self.sum1:
            num_unmixing_bands += 1
            print(f"Total unmixing bands (incl. sum_to_one): {num_unmixing_bands}")
        else:
            print(f"Total unmixing bands: {num_unmixing_bands}")

        ##Interpolate and process em libraries for use in building A matrices
        #####################################################################

        ##Build self._ems_interp
        ##Will create a list of three NxB matrices, each representing one of the
        ## endmember classes. N is the number of samples in the given class and
        ## B is the number of input image bands
        if self.dointerp:
            self.interpolate_ems(self.wl, self.fwhm)
        else:
            self._ems_interp = self._ems_raw

        ##Check that number of bands match
        for emnum, emdat in enumerate(self._ems_interp):
            if emdat.shape[band_axis] != len(self.wl):
                raise RuntimeError(
                    "Number of bands in interpolated em data"
                    f" for {self._em_names[emnum]}"
                    f" ({emdat.shape[band_axis]})"
                    " does not match number of image bands"
                    f" ({len(self.wl)})"
                )

        ##Compute averaged endmembers (must be done before processing)
        ##############################################################
        self.average_interpolated_ems()

        ##Build self._ems_proc
        ##Will create a list of three NxB matrices, each representing one of the
        ## endmember classes. N is the number of samples in the given class and
        ## B is the number of unmixing bands (sum of band regions + any extra
        ## for self.sum1 and self.shade)
        ########################################################################
        self.process_iterpolated_ems()

        ##Add a single-entry extra class for shade
        if self.shade:
            ##Shade is just a row of 0's
            tmpshade = np.full((1, num_unmixing_bands), 0.0)
            ##Last entry still 1.0 if sum1
            if self.sum1:
                tmpshade[0, -1] = 1.0
            ##Append it to the _ems_proc list
            self._ems_proc.append(tmpshade)
            ##Add an average shade to the _em_avg array
            self._em_avg = np.concatenate(
                [self._em_avg, np.full((1, num_selected_bands), 0.0)], axis=sample_axis
            )

        ##Prepare output map
        ####################

        ##Update the output meta to match args
        outmeta["dtype"] = self.dtype
        outmeta["driver"] = self.otype
        if self.nodata:
            outmeta["nodata"] = self.nodata
        else:
            outmeta.pop("nodata", None)
        outmeta["count"] = 2 * self.num_class + 1

        ##Fix output width and height based on nblocks
        outmeta["height"] = imgshape[0]
        outmeta["width"] = imgshape[1]
        outmeta["transform"] = inref.window_transform(
            rasterio.windows.Window(*reversed(imgoffset), *reversed(imgshape))
        )
        output_blocks = [
            (
                w_i,
                rasterio.windows.Window(
                    wind.col_off - imgoffset[1],
                    wind.row_off - imgoffset[0],
                    wind.width,
                    wind.height,
                ),
            )
            for w_i, wind in input_blocks
        ]

        ##Run AutoMCU by image blocks and write to output
        #################################################
        ##Zip the input and output blocks so we can iterate over pairs
        comb_blocks = list(
            zip([w for _, w in input_blocks], [w for _, w in output_blocks])
        )

        with rasterio.open(output_path, "w", **outmeta, **self.dopt) as oref:

            ##Update band descriptions
            ##########################
            banddesc = [*self._em_names]
            banddesc.extend([f"sd{cls}" for cls in self._em_names])
            banddesc.append("RMSE")
            print(f"Band names will be: {banddesc}")
            oref.descriptions = banddesc

            with rasterio.open(refl_path) as inref:
                for inwind, outwind in tqdm.tqdm(
                    comb_blocks, desc="Processing blocks", total=len(comb_blocks)
                ):
                    ###############################################
                    ##Read image data and filter out no data pixels
                    ###############################################

                    ##Make output shape tuple and count total block pixels
                    ######################################################
                    out_bands = 2 * self.num_class + 1
                    out_shape = (out_bands, inwind.height, inwind.width)
                    num_pix = inwind.height * inwind.width

                    ##Read the data for this block
                    ##############################
                    val_dat = (
                        inref.read(window=inwind).reshape(inref.count, -1) / input_scale
                    )

                    ##Filter for pixels with constant values across all
                    ## bands (i.e. no data)
                    ###################################################
                    isval = ~np.all(
                        val_dat == val_dat[0, :][np.newaxis, :], axis=sample_axis
                    )
                    ##How many valid pixels in this block window
                    num_valid = isval.sum()

                    ##Drop the nodata pixels
                    ########################
                    val_dat = val_dat[:, isval]
                    ##Shape (image_bands, num_valid)

                    ###########################################################
                    ##Unmix the valid data
                    ##Returns an array of concatenated mean coefficients, stdev
                    ## coefficients, and rmse
                    ## or None id val_dat is empty
                    ###########################################################
                    results = self.unmix_array(val_dat.T).T

                    ##########################
                    ##Prepare and write output
                    ##########################
                    out_dat = np.zeros((out_bands, num_pix), dtype=np.int16)
                    if results is not None:
                        out_dat[:, isval] = results * self.outscale
                    oref.write(out_dat.reshape(out_shape), window=outwind)
