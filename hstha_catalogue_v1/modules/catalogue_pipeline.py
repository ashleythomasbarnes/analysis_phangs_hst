from cat_imports import *

class PyHSTHACat:
    def __init__(self, galaxy, galaxy_hst=None, galaxyname_cat=None, output_dir_append=''):
        """
        Initialize the PyHSTHACat class with default parameters and file paths.
        """
        # ---------------------------------------------------------------------
        # 1. BASIC GALAXY NAMES & ROOT DIRECTORY
        # ---------------------------------------------------------------------
        self.galaxy = galaxy
        # If galaxy_hst is not provided, we'll assume it's the same as galaxy
        self.galaxy_hst = galaxy_hst if galaxy_hst else galaxy
        # Some association file naming uses 'galaxy_cat'
        self.galaxy_cat = galaxyname_cat if galaxyname_cat else galaxy
        # The main root directory where all data is stored
        self.rootdir = '/Users/abarnes/Dropbox/work/Smallprojects/galaxies'

        # ---------------------------------------------------------------------
        # 2. HST/MUSE DATA FILES
        # ---------------------------------------------------------------------
        self.hstha_file = (
            f"{self.rootdir}/data_hstha/{self.galaxy_hst}/hst_contsub/"
            f"{self.galaxy_hst}_hst_ha.fits"
        )
        self.hstha_err_file = (
            f"{self.rootdir}/data_hstha/{self.galaxy_hst}/hst_contsub/"
            f"{self.galaxy_hst}_hst_ha_err.fits"
        )
        self.muscat_file = (
            f"{self.rootdir}/data_hstha/{self.galaxy_hst}/muse/"
            f"{self.galaxy.upper()}_nebmask.fits"
        )

        # For the MUSE H-alpha flux map, we pick the first match from glob
        musha_glob = (
            f"{self.rootdir}/data_hstha/{self.galaxy_hst}/muse/"
            f"{self.galaxy.upper()}-*_MAPS.fits"
        )
        self.musha_file = glob(musha_glob)[0]

        # ---------------------------------------------------------------------
        # 3. OUTPUT DIRECTORIES FOR CUTOUTS, CATALOG, HDUs
        # ---------------------------------------------------------------------

        self.output_dir_append = output_dir_append

        self.cutout_dir = (
            f"{self.rootdir}/data_hstha_nebulae_catalogue/{self.galaxy_hst}/cutouts{self.output_dir_append}"
        )
        self.catalogue_dir = (
            f"{self.rootdir}/data_hstha_nebulae_catalogue/{self.galaxy_hst}/catalogue{self.output_dir_append}"
        )
        # self.cutouts_hdus_dir = (
        #     f"{self.rootdir}/data_hstha_nebulae_catalogue/{self.galaxy_hst}/cutouts_hdus{self.output_dir_append}"
        # )

        # ---------------------------------------------------------------------
        # 4. RERUN FLAGS
        # ---------------------------------------------------------------------
        # These flags let you skip or repeat certain steps.
        self.rerun_all = False
        self.rerun_regions = False
        self.rerun_masking = False
        self.rerun_cutouts = False
        self.rerun_cutouts_associations = False

        # ---------------------------------------------------------------------
        # 5. KEY ANCILLARY FILES FOR REGIONS, TABLES
        # ---------------------------------------------------------------------
        # Region definitions for the H-alpha cutouts
        self.regions_file = f"{self.cutout_dir}/sample.reg"
        self.regions_pickle_file = f"{self.cutout_dir}/sample.pickel"
        
        self.regions_file_hstonly = f"{self.cutout_dir}/hdus_catalog_mask_hstonly.pickel"
        # MUSE sample table & the main galaxy properties table
        self.sample_table_file = f"{self.rootdir}/data_misc/sample_table/phangs_sample_table_v1p6.fits"
        self.muscat_table_file = f"{self.rootdir}/data_misc/Nebulae_catalogue_v3/Nebulae_catalogue_v3.fits"

        # ---------------------------------------------------------------------
        # 6. PLACEHOLDERS FOR LOADED HDUs, TABLES, FINAL PROPERTIES
        # ---------------------------------------------------------------------
        self.hstha_hdu = None
        self.hstha_err_hdu = None
        self.musha_hdu = None
        self.muscat_hdu = None
        self.hdus = []
        self.regions = None

        self.hstha_err = None
        self.props = None  # final nebula catalog (props_all.fits)
        self.props_associations = None  # final association-augmented catalog

        self.muscat_table = None
        self.sample_table = None

        # ---------------------------------------------------------------------
        # 7. STELLAR ASSOCIATION FILES (MULTI-SCALE, VISUAL/NUV SELECT)
        # ---------------------------------------------------------------------
        # "v" stands for "visually selected" associations; "nuv" for NUV-selected
        # "08pc", "16pc", "32pc", "64pc" are different physical scale windows.

        def assoc_path(fmt, scale, selection):
            """
            Returns the absolute path for a given scale (e.g. '8pc', '16pc')
            and selection type ('v' -> vselect, 'nuv' -> nuvselect).
            `fmt` can be 'table' or 'mask'.
            
            Examples:
                assoc_path('table', '8pc', 'v')
                assoc_path('mask',  '32pc','nuv')
            """
            if fmt == 'table':
                return (
                    f"{self.rootdir}/data_misc/multiscale_stellar_associations/"
                    f"{self.galaxy_cat}/{selection}select/ws{scale}/"
                    f"PHANGS_IR4_hst_wfc3_{self.galaxy_cat}_v1p3_multi_assoc-"
                    f"{selection}select-ws{scale}-main.fits"
                )
            elif fmt == 'mask':
                return (
                    f"{self.rootdir}/data_misc/multiscale_stellar_associations/"
                    f"{self.galaxy_cat}/{selection}select/ws{scale}/"
                    f"PHANGS_IR4_hst_wfc3_{self.galaxy_cat}_v1p3_multi_assoc-"
                    f"{selection}select-ws{scale}-idmask.fits"
                )
            else:
                raise ValueError("Unknown fmt. Use 'table' or 'mask'.")

        # Tables for visually selected associations (v_*) & NUV selected (nuv_*)
        self.association_table_file_v_08pc   = assoc_path('table', '8pc',  'v')
        self.association_table_file_v_16pc   = assoc_path('table', '16pc', 'v')
        self.association_table_file_v_32pc   = assoc_path('table', '32pc', 'v')
        self.association_table_file_v_64pc   = assoc_path('table', '64pc', 'v')

        self.association_table_file_nuv_08pc = assoc_path('table', '8pc',  'nuv')
        self.association_table_file_nuv_16pc = assoc_path('table', '16pc', 'nuv')
        self.association_table_file_nuv_32pc = assoc_path('table', '32pc', 'nuv')
        self.association_table_file_nuv_64pc = assoc_path('table', '64pc', 'nuv')

        # Masks for visually selected associations & NUV selected
        self.association_mask_file_v_08pc   = assoc_path('mask', '8pc',  'v')
        self.association_mask_file_v_16pc   = assoc_path('mask', '16pc', 'v')
        self.association_mask_file_v_32pc   = assoc_path('mask', '32pc', 'v')
        self.association_mask_file_v_64pc   = assoc_path('mask', '64pc', 'v')

        self.association_mask_file_nuv_08pc = assoc_path('mask', '8pc',  'nuv')
        self.association_mask_file_nuv_16pc = assoc_path('mask', '16pc', 'nuv')
        self.association_mask_file_nuv_32pc = assoc_path('mask', '32pc', 'nuv')
        self.association_mask_file_nuv_64pc = assoc_path('mask', '64pc', 'nuv')

        # The final H-alpha nebula mask that will be used to index association data
        self.catalog_mask_file = (
            f"{self.catalogue_dir}/{self.galaxy}_mask.fits"
        )



    # -------------------------------------------------------------------------
    # A. MAKE PATHS
    # -------------------------------------------------------------------------
    def make_paths(self):
        """
        Create necessary directories for outputs. 
        Deletes and recreates root directory if `rerun_all` is True.
        """
        root_dir = f"{self.rootdir}/data_hstha_nebulae_catalogue/{self.galaxy_hst}/"
        print('[Info] Outputting to the following:')
        print(root_dir)

        # Remove and recreate root directory if rerun_all is set
        if self.rerun_all:
            os.system(f'rm -rf {root_dir}')

        # Create root directory if it doesn't exist
        if not os.path.isdir(root_dir):
            os.mkdir(root_dir)

        # Create subdirectories for cutouts, catalogues, and HDUs
        for path in [self.cutout_dir, self.catalogue_dir]:
            if not os.path.isdir(path):
                os.mkdir(path)

    # -------------------------------------------------------------------------
    # B. LOAD FILES
    # -------------------------------------------------------------------------
    def load_files(self):
        """
        Load HST & MUSE FITS files, handle NaNs, convert to float32.
        Also loads the sample and MUSE property tables into memory.
        """
        # -- Load HST and MUSE data --
        self.hstha_hdu     = fits.open(self.hstha_file)[0]
        self.hstha_err_hdu = fits.open(self.hstha_err_file)[0]
        self.musha_hdu     = fits.open(self.musha_file)['HA6562_FLUX']
        self.muscat_hdu    = fits.open(self.muscat_file)[0]

        # Handle NaNs in the MUSE H-alpha flux
        self.musha_hdu.data[np.isnan(self.musha_hdu.data)] = -100

        # Replace -1 with NaN in the MUSE nebmask
        self.muscat_hdu.data = np.array(self.muscat_hdu.data, dtype=float)
        self.muscat_hdu.data[self.muscat_hdu.data == -1] = np.nan

        # Convert all HDUs to float32 format
        hdus = [self.hstha_hdu, self.musha_hdu, self.muscat_hdu, self.hstha_err_hdu]
        hdus_converted = []
        for hdu in hdus:
            hdus_converted.append(cat_misc.convert_to_float32(hdu.copy()))
        self.hstha_hdu, self.musha_hdu, self.muscat_hdu, self.hstha_err_hdu = hdus_converted

        # Load the sample (galaxy) table & MUSE property table
        self.muscat_table = cat_misc.get_museprops(self.galaxy, self.muscat_table_file)
        self.sample_table = cat_misc.get_galaxyprops(self.galaxy, self.sample_table_file)
        
        # Remove HDUS to free up space
        del hdus, hdus_converted
        gc.collect()

    # -------------------------------------------------------------------------
    # C. GET REGIONS
    # -------------------------------------------------------------------------
    def get_regions(self):
        """
        Retrieve DS9 regions from disk (if they exist) or create them 
        from the MUSE table and store them as a pickle for future use.
        """
        if os.path.exists(self.regions_file) and not self.rerun_regions:
            print(f'[INFO] Using existing region file: {self.regions_file}')
            self.regions = cat_misc.load_pickle(self.regions_pickle_file)
        else:
            print('[INFO] Generating new DS9 region file...')
            cat_cutouts.get_ds9regions_all(self.muscat_table, outputfile=self.regions_file)
            self.regions = cat_cutouts.get_regions(self.regions_file)
            cat_misc.save_pickle(self.regions, self.regions_pickle_file)

    # -------------------------------------------------------------------------
    # D. GET CUTOUTS
    # -------------------------------------------------------------------------
    def get_cutouts(self):
        """
        Generate cutouts for each HDU and save them to disk. 
        Loads them into a single dictionary for quick reuse later.
        """
        names = ['hstha_hdu', 'musha_hdu', 'muscat_hdu']
        self.hdus_cutouts = {}

        combined_path = f"{self.cutout_dir}/hdus_all.pickel"
        if os.path.exists(combined_path) and not self.rerun_cutouts:
            print(f'[INFO] All cutouts already exist. Loading from file {combined_path}')
            self.hdus_cutouts = cat_misc.load_pickle(combined_path)
        else:
            # Generate cutouts and pickle them individually
            for hdu, name in zip([self.hstha_hdu, self.musha_hdu, self.muscat_hdu], names):
                
                print(f'[INFO] Generating cutouts for {name}...')
                hdu_cutouts = cat_cutouts.get_croppeddata_all(hdu, self.regions)
                pickle_path = f"{self.cutout_dir}/{name}.pickel"
                cat_misc.save_pickle(hdu_cutouts, pickle_path)
                
                # Save space
                del hdu_cutouts
                gc.collect()

            # Reload them back into a single dictionary
            for name in names:
                pickle_path = f"{self.cutout_dir}/{name}.pickel"
                self.hdus_cutouts[name] = cat_misc.load_pickle(pickle_path)

            # Store the combined dictionary for future runs
            print(f'[INFO] All HDU cutouts saved to {combined_path}')
            cat_misc.save_pickle(self.hdus_cutouts, combined_path)

    # -------------------------------------------------------------------------
    # E. MAKE CATALOGUE
    # -------------------------------------------------------------------------
    def make_catalogue(self):
        """
        Make a final catalogue (props_all.fits) by:
          - Interpolating masks,
          - Computing noise properties from the H-alpha error map,
          - Loading region-based cutouts / applying masks,
          - Computing dendrogram-based properties,
          - Merging with MUSE tables,
          - Correcting fluxes/luminosities,
          - Saving outputs (catalog and mask FITS files).
        """
        print("[CATALOGUE] Running catalogue generation for:", self.galaxy)

        # Ensure the catalogue directory exists
        cat_misc.checkmakedir(self.catalogue_dir)

        # Get the median noise from the HST H-alpha error map if not done before
        if self.hstha_err is None:
            self.hstha_err = np.nanmedian(self.hstha_err_hdu.data)

        # Decide whether to load or create masked HDUs
        hdus_file = f'{self.cutout_dir}/hdus_all_withmasked.pickel'
        if os.path.exists(hdus_file) and (not self.rerun_masking):
            muscat_regionIDs = self.muscat_table['region_ID']
            hdus = cat_misc.load_pickle(hdus_file)
        else:
            muscat_regionIDs = self.muscat_table['region_ID']
            hdus = cat_mask.get_maskedhdus(self.hdus_cutouts, self.regions, muscat_regionIDs)
            cat_misc.save_pickle(hdus, hdus_file)

        props_all       = []
        hdus_mask       = []
        hdus_mask_id    = []
        hdus_data_masked= []

        # Loop over each region in the MUSE table
        for i in tqdm(range(len(muscat_regionIDs)), desc='Get sources:', position=0):
            regionID = np.int16(muscat_regionIDs[i])
            
            data   = hdus['hstha_hdu_masked'][i].data.copy()
            header = hdus['hstha_hdu_masked'][i].header.copy()

            # 1) Create threshold-based masks
            mask_low       = cat_mask.get_threshmask(data, self.hstha_err, thresh=1)
            mask_low_prune = cat_mask.get_prunemask(mask_low, thresh=50)
            mask_high      = cat_mask.get_threshmask(data, self.hstha_err, thresh=3)
            mask_grow      = ndimage.binary_dilation(mask_high, iterations=-1, mask=mask_low_prune)
            mask_prune     = cat_mask.get_prunemask(mask_grow, thresh=9)
            mask_filled    = cat_mask.get_filled_outer_contour(mask_prune)
            mask_final     = mask_filled.copy()

            # 2) Save the mask
            hdu_mask = fits.PrimaryHDU(mask_final.astype(np.float32), header)
            hdus_mask.append(hdu_mask)

            # 3) Create an ID-based mask
            hdu_mask_id = hdu_mask.copy()
            hdu_mask_id.data = (hdu_mask_id.data * i).astype(np.float32)
            hdu_mask_id.data[~mask_final] = -1
            hdus_mask_id.append(hdu_mask_id)

            # 4) Mask the data itself
            data_masked = data.copy()
            data_masked[~mask_final] = np.nan
            hdu_data_masked = fits.PrimaryHDU(data_masked, header)
            hdus_data_masked.append(hdu_data_masked)

            # If everything is masked, skip
            if np.nansum(~np.isnan(data_masked)) == 0:
                continue

            # 5) Determine pixel scale
            if 'CDELT1' in header and 'CDELT2' in header:
                pixsize = np.array([abs(header['CDELT1']), abs(header['CDELT2'])]).mean() * au.degree
            else:
                pixsize = 1

            # fallback if these are placeholders
            if pixsize.value == 1:
                if 'CD1_1' in header and 'CD2_2' in header:
                    pixsize = np.array([abs(header['CD1_1']), abs(header['CD2_2'])]).mean() * au.degree
                    # print('[INFO] Pixel scale taken as CD1_1, CD2_2')
                elif 'PC1_1' in header and 'PC2_2' in header:
                    pixsize = np.array([abs(header['PC1_1']), abs(header['PC2_2'])]).mean() * au.degree
                    # print('[INFO] Pixel scale taken as PC1_1, PC2_2')
            # else: 
            #     print('[INFO] Pixel scale taken as CDELT')

            # 6) Basic region properties (flux, area, etc.)
            npix      = np.nansum(mask_final) * au.pix
            flux      = np.nansum(data_masked) * au.erg/au.s/au.cm**2
            flux_err  = np.sqrt(npix.value)*self.hstha_err * flux.unit

            area_exact  = npix.value * np.abs(pixsize.to(au.arcsec)**2)
            radius_circ = np.sqrt(area_exact / np.pi).to('arcsec')

            flux_max  = np.nanmax(data_masked) * flux.unit
            flux_min  = np.nanmin(data_masked) * flux.unit
            flux_mean = np.nanmean(data_masked) * flux.unit

            x_max, y_max = np.where(data_masked == flux_max.value)
            if len(x_max) > 1 or len(y_max) > 1:
                x_max, y_max = np.nanmean(x_max), np.nanmean(y_max)
            else:
                x_max, y_max = x_max[0], y_max[0]

            data_zeros = data_masked.copy()
            data_zeros[np.isnan(data_zeros)] = 0
            x_com, y_com = ndimage.center_of_mass(data_zeros)

            # Convert array indices to RA,Dec
            wcs = WCS(header)
            ra_max, dec_max = wcs.array_index_to_world_values([[y_max, x_max]])[0] * au.deg
            ra_com, dec_com = wcs.array_index_to_world_values([[y_com, x_com]])[0] * au.deg

            # 7) Assemble into a table row
            table_data = [
                regionID, x_max, y_max, x_com, y_com, 
                ra_max, dec_max, ra_com, dec_com,
                npix, flux, flux_err, area_exact, radius_circ, 
                flux_max, flux_min, flux_mean
            ]
            table_data = [np.array(td) for td in table_data]
            table_names = [
                'region_ID', 'x_max', 'y_max', 'x_com', 'y_com',
                'ra_max', 'dec_max', 'ra_com', 'dec_com',
                'npix', 'flux', 'flux_err', 'area_exact', 'radius_circ',
                'flux_max', 'flux_min', 'flux_mean'
            ]
            props = QTable(np.array(table_data), names=table_names)
            props['x_max'].unit       = au.pix
            props['y_max'].unit       = au.pix
            props['x_com'].unit       = au.pix
            props['y_com'].unit       = au.pix
            props['ra_max'].unit      = au.deg
            props['dec_max'].unit     = au.deg
            props['ra_com'].unit      = au.deg
            props['dec_com'].unit     = au.deg
            props['npix'].unit        = au.pix
            props['flux'].unit        = flux.unit
            props['flux_err'].unit    = flux.unit
            props['area_exact'].unit  = au.arcsec**2
            props['radius_circ'].unit = au.arcsec
            props['flux_max'].unit    = flux_max.unit
            props['flux_min'].unit    = flux_min.unit
            props['flux_mean'].unit   = flux_mean.unit

            pcperarcsec = cat_props.get_pcperarcsec(self.sample_table)
            props['radius_circ_pc'] = props['radius_circ'] * pcperarcsec

            # 8) Dendrogram analysis (lenient params)
            dendro = Dendrogram.compute(
                data_masked, 
                min_delta=-1000, 
                min_value=0, 
                min_npix=0, 
                wcs=wcs
            )
            metadata = {
                'data_unit': au.Jy/au.beam,
                'spatial_scale': pixsize.to('arcsec'),
                'beam_major': 0.1*au.arcsec,
                'beam_minor': 0.1*au.arcsec
            }
            props_dendro = pp_catalog(dendro.trunk, metadata, verbose=False)
            props_dendro = QTable(props_dendro)
            props_dendro = props_dendro[np.nanargmax(props_dendro['flux'])]  # largest flux structure

            # Save dendrogram peak coords
            ra_dendro, dec_dendro = wcs.array_index_to_world_values(
                [[props_dendro['x_cen'].value, props_dendro['y_cen'].value]]
            )[0] * au.deg

            props['x_mom']          = props_dendro['x_cen']
            props['y_mom']          = props_dendro['y_cen']
            props['ra_mom']         = ra_dendro
            props['dec_mom']        = dec_dendro
            props['area_ellipse']   = props_dendro['area_ellipse']
            props['major_sigma']    = props_dendro['major_sigma']
            props['minor_sigma']    = props_dendro['minor_sigma']
            props['mean_sigma']     = props_dendro['radius']
            props['position_angle'] = props_dendro['position_angle']
            props['mean_sigma_pc']  = props['mean_sigma'] * pcperarcsec

            # 9) Dendrogram complexity with higher thresholds
            dendro = Dendrogram.compute(
                data_masked,
                min_delta=(self.hstha_err*3),
                min_value=(self.hstha_err),
                min_npix=9,
                wcs=wcs
            )
            dendro_IDs = np.unique(dendro.index_map.data)
            dendro_IDs = [dID for dID in dendro_IDs if dID > -1]

            dendro_complex = len(dendro_IDs)
            dendro_complex_leaves = len(dendro.leaves)
            complexity_rms = np.sqrt(np.nanmean((data_masked.copy())**2))
            complexity_std = np.nanstd(data_masked.copy())

            props['complexity_score'] = dendro_complex
            props['complexity_score_leaves'] = dendro_complex_leaves
            props['complexity_rms'] = complexity_rms
            props['complexity_std'] = complexity_std

            # 10) Edge & Touch flags
            flag_edge = np.isnan(hdus['hstha_hdu_masked_ones'][i].data).any()*1
            props.add_column(Column(flag_edge, name='flag_edge_hst'))

            mask_touch = (ndimage.binary_dilation(mask_final) & ~mask_final)
            flag_touch = (np.nansum(np.isnan(data[mask_touch])) > 0)*1
            props.add_column(Column(flag_touch, name='flag_touch_hst'))

            # Add this region's measurements
            props_all.append(props)

        # Clean up memory
        del data, header, dendro, metadata, props_dendro, dendro_IDs, hdu_data_masked, hdu_mask, mask_low, mask_low_prune, mask_high, mask_grow, mask_prune, mask_filled
        gc.collect()

        # Combine all region properties into one table
        props_all = vstack(props_all)

        # Merge with MUSE data on region_ID
        self.muscat_table = cat_misc.get_museprops(self.galaxy, self.muscat_table_file)
        self.muscat_table_rename = self.muscat_table.copy()
        for column in self.muscat_table_rename.colnames:
            self.muscat_table_rename.rename_column(column, column+'_MUSE')
        self.muscat_table_rename.rename_column('gal_name_MUSE', 'gal_name')
        self.muscat_table_rename.rename_column('region_ID_MUSE', 'region_ID')
        self.muscat_table_rename.rename_column('Lum_HA6562_CORR_MUSE', 'HA6562_LUMINOSITY_MUSE')

        props_all_final = join(props_all, self.muscat_table_rename, keys='region_ID')

        # Correct fluxes & luminosities
        props_all_final['flux_corr'] = cat_props.correct_ha_flux(props_all_final, props_all_final['flux'])
        props_all_final['flux_err_corr'] = cat_props.correct_ha_flux(props_all_final, props_all_final['flux_err'])
        dist_val = self.sample_table['dist'][0]
        props_all_final['lum_hst'] = cat_props.calculate_luminosity(
            props_all_final['flux_corr'] * 1e-20, dist_val
        )
        props_all_final['lum_err_hst'] = cat_props.calculate_luminosity(
            props_all_final['flux_err_corr'] * 1e-20, dist_val
        )
        props_all_final['region_circ_rad_pc_MUSE'] = cat_props.calculate_radius(
            props_all_final['region_circ_rad_MUSE'], dist_val
        )

        # Rename columns for clarity
        props_all_final.rename_column('flux',          'HA6562_FLUX_HST')
        props_all_final.rename_column('flux_err',      'HA6562_FLUX_ERR_HST')
        props_all_final.rename_column('flux_corr',     'HA6562_FLUX_CORR_HST')
        props_all_final.rename_column('flux_err_corr', 'HA6562_FLUX_ERR_CORR_HST')
        props_all_final.rename_column('lum_hst',       'HA6562_LUMINOSITY_HST')
        props_all_final.rename_column('lum_err_hst',   'HA6562_LUMINOSITY_ERR_HST')

        # Store final catalog in memory & on disk
        self.props = props_all_final

        props_file = f'{self.catalogue_dir}/props_all.fits'
        self.props.write(props_file, overwrite=True)
        print(f'[CATALOGUE] Final properties table saved to: {props_file}')

        # Save masks
        mask_file = f'{self.catalogue_dir}/{self.galaxy}_mask.fits'
        cat_mask.get_hdumask(self.hstha_hdu, hdus_mask_id, outputfile=mask_file)

        complexity_file = f'{self.catalogue_dir}/{self.galaxy}_complexity.fits'
        cat_mask.get_hducomplex(self.props, inputfile=mask_file, outputfile=complexity_file)

        regions_out = f'{self.catalogue_dir}/{self.galaxy}_mask'
        cat_mask.get_ds9regions(props_all, outputfile=regions_out)

        print(f'[CATALOGUE] Masks saved to {mask_file} and {complexity_file}...')
        print(f'[CATALOGUE] DS9 regions output to {regions_out}...')
        print('[CATALOGUE] Done!')

    # -------------------------------------------------------------------------
    # F. MAKE ASSOCIATIONS
    # -------------------------------------------------------------------------

    def make_associations(self):
        """
        Build an extended table that associates each nebula region_ID with 
        stellar associations at multiple physical scales (v_08pc, nuv_16pc, etc.).
        """
        # 1. Define all association files
        assoc_files = {
            "v_08pc": (self.association_table_file_v_08pc, self.association_mask_file_v_08pc),
            "v_16pc": (self.association_table_file_v_16pc, self.association_mask_file_v_16pc),
            "v_32pc": (self.association_table_file_v_32pc, self.association_mask_file_v_32pc),
            "v_64pc": (self.association_table_file_v_64pc, self.association_mask_file_v_64pc),
            "nuv_08pc": (self.association_table_file_nuv_08pc, self.association_mask_file_nuv_08pc),
            "nuv_16pc": (self.association_table_file_nuv_16pc, self.association_mask_file_nuv_16pc),
            "nuv_32pc": (self.association_table_file_nuv_32pc, self.association_mask_file_nuv_32pc),
            "nuv_64pc": (self.association_table_file_nuv_64pc, self.association_mask_file_nuv_64pc),
        }

        print(f"[ASSOCIATIONS] Checking association files for galaxy {self.galaxy} ...")

        # Check for missing files
        valid_files = {}
        for scale, (table_file, mask_file) in assoc_files.items():

            cat_misc.get_unpack(mask_file)  # Unpack mask if needed

            if os.path.isfile(table_file) and os.path.isfile(mask_file):
                valid_files[scale] = (table_file, mask_file)
                print(f"[INFO] Files found for {scale}.")

            else:
                print(f"[WARNING] Files missing for {scale}. Skipping.")

        if not valid_files:
            print("[ERROR] No valid association files found. Exiting.")
            return

        # 2. Read each association FITS table & mask for valid files
        props_associations = {}
        hdus_association_masks = {}

        for scale, (table_file, mask_file) in valid_files.items():
            props_associations[scale] = QTable.read(table_file)
            hdus_association_masks[scale] = fits.open(mask_file)[0]

        # TODO - remove this.. 
        # 3. Load the main nebula catalog (props_all.fits) & region definitions
        props_all_path = f"{self.catalogue_dir}/props_all.fits"
        props_all = QTable.read(props_all_path)
        regions = cat_misc.load_pickle(self.regions_pickle_file)
        hdu_catalog_mask = fits.open(self.catalog_mask_file)[0]

        # Potentially filter 'regions' to only those in props_all
        region_IDs = props_all['region_ID'].data
        for key in regions.keys():
            subset = [int(rid) for rid in region_IDs]
            regions[key] = [regions[key][s] for s in subset]

        # 4. Create or load association-mask cutouts
        hdus_association_masks_new = {}
        if self.rerun_cutouts_associations:
            print("[ASSOCIATIONS] Generating new cutouts for association masks...")
            
            # Save the new catalog mask cutout with only regions in HST catalogue - not all MUSE 
            hdus_catalog_mask_new = cat_cutouts.get_croppeddata_all(hdu_catalog_mask, regions)
            cat_misc.save_pickle(hdus_catalog_mask_new, f"{self.cutout_dir}/hdus_catalog_mask_new.pickel")

            # Save all association masks at valid scales
            for scale, hdu_mask in hdus_association_masks.items():
                hdus_association_masks_new[scale] = cat_cutouts.get_croppeddata_all(hdu_mask, regions)
                cat_misc.save_pickle(hdus_association_masks_new[scale], f"{self.cutout_dir}/hdus_association_mask_new_{scale}.pickel")

        else:
            print("[ASSOCIATIONS] Loading existing association mask cutouts...")

            # Load the new catalog mask cutout with only regions in HST catalogue - not all MUSE
            hdus_catalog_mask_new = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_catalog_mask_new.pickel")

            # Load all association masks
            for scale in valid_files.keys():
                hdus_association_masks_new[scale] = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_{scale}.pickel")

        # 5. For each region, determine which association IDs overlap
        props_associations_all = []
        for i in tqdm(range(len(props_all['region_ID'])), desc='Match associations:', position=0):
            region_ID = props_all['region_ID'][i]

            # Gather matching IDs for all valid scales
            props_associations_ = {}
            for j, scale in enumerate(valid_files.keys()):

                association_IDs = cat_associations.get_association_IDs(region_ID, i, hdus_catalog_mask_new, hdus_association_masks_new[scale])
                props_associations_[scale] = cat_associations.get_all_associations(region_ID, association_IDs, props_associations[scale], scale)

                if j == 0:
                    props_associations_all_ = props_associations_[scale]

                else:
                    props_associations_all_ = join(props_associations_all_, props_associations_[scale], keys='region_ID')

            props_associations_all.append(props_associations_all_)

        # 6. Stack them into one big table, then join with the nebula table
        props_associations_stacked = vstack(props_associations_all)
        self.props_associations = join(props_all, props_associations_stacked, keys='region_ID')

        # 7. Save final association-augmented catalog
        out_fits = f"{self.catalogue_dir}/props_all_association.fits"
        self.props_associations.write(out_fits, overwrite=True)
        print(f"[ASSOCIATIONS] Final association table saved to: {out_fits}")

        # Optionally repack the files to save space
        for scale, (table_file, mask_file) in valid_files.items():
            cat_misc.clean_unpack(mask_file)
    
    # # F. MAKE ASSOCIATIONS
    # # -------------------------------------------------------------------------
    # def make_associations(self):
    #     """
    #     Build an extended table that associates each nebula region_ID with 
    #     stellar associations at multiple physical scales (v_08pc, nuv_16pc, etc.).
        
    #     Steps:
    #       1. Check if association files exist, unpack them if needed.
    #       2. Read in the association tables & masks (via fits).
    #       3. Load the main nebula props table ('props_all.fits') and region definitions.
    #       4. Create or load cutouts for each association mask if 'self.rerun_cutouts_associations' is True.
    #       5. For each nebula region, compute which association IDs fall inside.
    #       6. Join them into one final table, 'props_all_association.fits'.
    #       7. Clean up (e.g., remove unpacked temp files).
    #     """
    #     # 1. Define all association files
    #     assoc_files = [
    #         self.association_table_file_v_08pc,   self.association_table_file_v_16pc,
    #         self.association_table_file_v_32pc,   self.association_table_file_v_64pc,
    #         self.association_table_file_nuv_08pc, self.association_table_file_nuv_16pc,
    #         self.association_table_file_nuv_32pc, self.association_table_file_nuv_64pc,
    #         self.association_mask_file_v_08pc,    self.association_mask_file_v_16pc,
    #         self.association_mask_file_v_32pc,    self.association_mask_file_v_64pc,
    #         self.association_mask_file_nuv_08pc,  self.association_mask_file_nuv_16pc,
    #         self.association_mask_file_nuv_32pc,  self.association_mask_file_nuv_64pc
    #     ]

    #     print(f"[ASSOCIATIONS] Checking association files for galaxy {self.galaxy} ...")
    #     for file in assoc_files:
    #         cat_misc.get_unpack(file)  # Unpack if needed
            
    #         # if os.path.isfile(file):
    #         #     print(f"[INFO] File found: {file}")
    #         # else:
    #         #     print(f"[WARNING] File not found: {file}")

    #     # 2. Read each association FITS table & mask
    #     props_associations_v_08pc   = QTable.read(self.association_table_file_v_08pc)
    #     props_associations_v_16pc   = QTable.read(self.association_table_file_v_16pc)
    #     props_associations_v_32pc   = QTable.read(self.association_table_file_v_32pc)
    #     props_associations_v_64pc   = QTable.read(self.association_table_file_v_64pc)

    #     props_associations_nuv_08pc = QTable.read(self.association_table_file_nuv_08pc)
    #     props_associations_nuv_16pc = QTable.read(self.association_table_file_nuv_16pc)
    #     props_associations_nuv_32pc = QTable.read(self.association_table_file_nuv_32pc)
    #     props_associations_nuv_64pc = QTable.read(self.association_table_file_nuv_64pc)

    #     hdu_association_mask_v_08pc   = fits.open(self.association_mask_file_v_08pc)[0]
    #     hdu_association_mask_v_16pc   = fits.open(self.association_mask_file_v_16pc)[0]
    #     hdu_association_mask_v_32pc   = fits.open(self.association_mask_file_v_32pc)[0]
    #     hdu_association_mask_v_64pc   = fits.open(self.association_mask_file_v_64pc)[0]

    #     hdu_association_mask_nuv_08pc = fits.open(self.association_mask_file_nuv_08pc)[0]
    #     hdu_association_mask_nuv_16pc = fits.open(self.association_mask_file_nuv_16pc)[0]
    #     hdu_association_mask_nuv_32pc = fits.open(self.association_mask_file_nuv_32pc)[0]
    #     hdu_association_mask_nuv_64pc = fits.open(self.association_mask_file_nuv_64pc)[0]

    #     # 3. Load the main nebula catalog (props_all.fits) & region definitions
    #     props_all_path = f"{self.catalogue_dir}/props_all.fits"
    #     props_all = QTable.read(props_all_path)
    #     regions = cat_misc.load_pickle(self.regions_pickle_file)
    #     hdu_catalog_mask = fits.open(self.catalog_mask_file)[0]

    #     # Potentially filter 'regions' to only those in props_all
    #     region_IDs = props_all['region_ID'].data
    #     for key in regions.keys():
    #         subset = [int(rid) for rid in region_IDs]
    #         regions[key] = [regions[key][s] for s in subset]

    #     # 4. Create or load association-mask cutouts
    #     if self.rerun_cutouts_associations:
    #         print("[ASSOCIATIONS] Generating new cutouts for association masks...")

    #         hdus_association_mask_new_v_08pc   = cat_cutouts.get_croppeddata_all(hdu_association_mask_v_08pc,   regions)
    #         hdus_association_mask_new_v_16pc   = cat_cutouts.get_croppeddata_all(hdu_association_mask_v_16pc,   regions)
    #         hdus_association_mask_new_v_32pc   = cat_cutouts.get_croppeddata_all(hdu_association_mask_v_32pc,   regions)
    #         hdus_association_mask_new_v_64pc   = cat_cutouts.get_croppeddata_all(hdu_association_mask_v_64pc,   regions)

    #         hdus_association_mask_new_nuv_08pc = cat_cutouts.get_croppeddata_all(hdu_association_mask_nuv_08pc, regions)
    #         hdus_association_mask_new_nuv_16pc = cat_cutouts.get_croppeddata_all(hdu_association_mask_nuv_16pc, regions)
    #         hdus_association_mask_new_nuv_32pc = cat_cutouts.get_croppeddata_all(hdu_association_mask_nuv_32pc, regions)
    #         hdus_association_mask_new_nuv_64pc = cat_cutouts.get_croppeddata_all(hdu_association_mask_nuv_64pc, regions)

    #         hdus_catalog_mask_new = cat_cutouts.get_croppeddata_all(hdu_catalog_mask, regions)

    #         # Save all these new cutouts
    #         cat_misc.save_pickle(hdus_association_mask_new_v_08pc,   f"{self.cutout_dir}/hdus_association_mask_new_v_08pc.pickel")
    #         cat_misc.save_pickle(hdus_association_mask_new_v_16pc,   f"{self.cutout_dir}/hdus_association_mask_new_v_16pc.pickel")
    #         cat_misc.save_pickle(hdus_association_mask_new_v_32pc,   f"{self.cutout_dir}/hdus_association_mask_new_v_32pc.pickel")
    #         cat_misc.save_pickle(hdus_association_mask_new_v_64pc,   f"{self.cutout_dir}/hdus_association_mask_new_v_64pc.pickel")

    #         cat_misc.save_pickle(hdus_association_mask_new_nuv_08pc, f"{self.cutout_dir}/hdus_association_mask_new_nuv_08pc.pickel")
    #         cat_misc.save_pickle(hdus_association_mask_new_nuv_16pc, f"{self.cutout_dir}/hdus_association_mask_new_nuv_16pc.pickel")
    #         cat_misc.save_pickle(hdus_association_mask_new_nuv_32pc, f"{self.cutout_dir}/hdus_association_mask_new_nuv_32pc.pickel")
    #         cat_misc.save_pickle(hdus_association_mask_new_nuv_64pc, f"{self.cutout_dir}/hdus_association_mask_new_nuv_64pc.pickel")

    #         cat_misc.save_pickle(hdus_catalog_mask_new,              f"{self.cutout_dir}/hdus_catalog_mask_new.pickel")
    #     else:
    #         print("[ASSOCIATIONS] Loading existing association mask cutouts...")

    #         hdus_association_mask_new_v_08pc   = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_v_08pc.pickel")
    #         hdus_association_mask_new_v_16pc   = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_v_16pc.pickel")
    #         hdus_association_mask_new_v_32pc   = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_v_32pc.pickel")
    #         hdus_association_mask_new_v_64pc   = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_v_64pc.pickel")

    #         hdus_association_mask_new_nuv_08pc = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_nuv_08pc.pickel")
    #         hdus_association_mask_new_nuv_16pc = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_nuv_16pc.pickel")
    #         hdus_association_mask_new_nuv_32pc = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_nuv_32pc.pickel")
    #         hdus_association_mask_new_nuv_64pc = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_association_mask_new_nuv_64pc.pickel")

    #         hdus_catalog_mask_new = cat_misc.load_pickle(f"{self.cutout_dir}/hdus_catalog_mask_new.pickel")

    #     # 5. For each region, figure out which association IDs fall inside
    #     props_associations_all = []
    #     for i in tqdm(range(len(props_all['region_ID'])), desc='Match associations:', position=0):
    #         region_ID = props_all['region_ID'][i]

    #         # v_08pc scale
    #         association_IDs_v_08pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_v_08pc
    #         )
    #         # ... similarly for all scales (v_16pc, v_32pc, ...)
    #         association_IDs_v_16pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_v_16pc
    #         )
    #         association_IDs_v_32pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_v_32pc
    #         )
    #         association_IDs_v_64pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_v_64pc
    #         )

    #         association_IDs_nuv_08pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_nuv_08pc
    #         )
    #         association_IDs_nuv_16pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_nuv_16pc
    #         )
    #         association_IDs_nuv_32pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_nuv_32pc
    #         )
    #         association_IDs_nuv_64pc = cat_associations.get_association_IDs(
    #             region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new_nuv_64pc
    #         )

    #         # Now gather up the matching association properties in each scale/selection
    #         props_associations_new_v_08pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_v_08pc, props_associations_v_08pc, 'v_08pc'
    #         )
    #         props_associations_new_v_16pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_v_16pc, props_associations_v_16pc, 'v_16pc'
    #         )
    #         props_associations_new_v_32pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_v_32pc, props_associations_v_32pc, 'v_32pc'
    #         )
    #         props_associations_new_v_64pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_v_64pc, props_associations_v_64pc, 'v_64pc'
    #         )

    #         props_associations_new_nuv_08pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_nuv_08pc, props_associations_nuv_08pc, 'nuv_08pc'
    #         )
    #         props_associations_new_nuv_16pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_nuv_16pc, props_associations_nuv_16pc, 'nuv_16pc'
    #         )
    #         props_associations_new_nuv_32pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_nuv_32pc, props_associations_nuv_32pc, 'nuv_32pc'
    #         )
    #         props_associations_new_nuv_64pc = cat_associations.get_all_associations(
    #             region_ID, association_IDs_nuv_64pc, props_associations_nuv_64pc, 'nuv_64pc'
    #         )

    #         # Combine them all into one row for this region
    #         props_associations_all_ = join(
    #             props_associations_new_v_08pc,
    #             props_associations_new_v_16pc,
    #             keys='region_ID'
    #         )
    #         props_associations_all_ = join(
    #             props_associations_all_,
    #             props_associations_new_v_32pc,
    #             keys='region_ID'
    #         )
    #         props_associations_all_ = join(
    #             props_associations_all_,
    #             props_associations_new_v_64pc,
    #             keys='region_ID'
    #         )
    #         props_associations_all_ = join(
    #             props_associations_all_,
    #             props_associations_new_nuv_08pc,
    #             keys='region_ID'
    #         )
    #         props_associations_all_ = join(
    #             props_associations_all_,
    #             props_associations_new_nuv_16pc,
    #             keys='region_ID'
    #         )
    #         props_associations_all_ = join(
    #             props_associations_all_,
    #             props_associations_new_nuv_32pc,
    #             keys='region_ID'
    #         )
    #         props_associations_all_ = join(
    #             props_associations_all_,
    #             props_associations_new_nuv_64pc,
    #             keys='region_ID'
    #         )

    #         props_associations_all.append(props_associations_all_)

    #     # 6. Stack them into one big table, then join with the nebula table
    #     props_associations_stacked = vstack(props_associations_all)
    #     self.props_associations = join(props_all, props_associations_stacked, keys='region_ID')

    #     # 7. Save final association-augmented catalog & optional cleanup
    #     out_fits = f"{self.catalogue_dir}/props_all_association.fits"
    #     self.props_associations.write(out_fits, overwrite=True)
    #     print(f"[ASSOCIATIONS] Final association table saved to: {out_fits}")

    #     # Optionally repack the files to save space
    #     for file in assoc_files:
    #         cat_misc.clean_unpack(file)

    #     print("[ASSOCIATIONS] Done building the multi-scale association catalog!") 