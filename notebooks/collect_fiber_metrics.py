#!/usr/bin/env python

# find $SAS_BASE_DIR/sdsswork/lvm/spectro/redux/1.1.0 -type f -name 'lvmSFrame-*.fits' > drp_1.1.0_SFrames_list.txt

import os
import click
import numpy as np
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import warnings

# Ignore warning in np.nanmedian if all values are nan in some slice
warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="RuntimeWarning: Mean of empty slice.")
warnings.filterwarnings("ignore", message="RuntimeWarning: Degrees of freedom <= 0 for slice.")

# This regions defined to compute the continuum and line metrics.
# See visual inspection notebook for details `Which_metrics_use.ipynb`.
WLINES = [
    dict(name='OI_5577', wave=5577.34, line_region=[-3, 3], cont_regions=[[-7, -3], [3, 7]]),
    dict(name='Na_5893', wave=5893, line_region=[-7, 6], cont_regions=[[-11, -7], [6, 10]]),
    dict(name='OI_6300', wave=6300.3, line_region=[-3, 3], cont_regions=[[-8, -5], [10, 14]]),
    dict(name='OH_7341', wave=7340.75, line_region=[-3, 3], cont_regions=[[-7, -3], [3, 7]]),
    dict(name='O2_8650', wave=8650, line_region=[-42, 62], cont_regions=[[-48, -42], [62, 72]]),
]

WCONTS = [
    dict(name='Bcont', wave=4620, cont_regions=[[-20, 20]]),
    dict(name='Rcont', wave=6765, cont_regions=[[-15, 15]]),
    dict(name='Zcont', wave=9193, cont_regions=[[-15, 15]]),
]


@click.command()
@click.argument('input_list', type=click.Path(exists=True, readable=True), required=True)
@click.option('--output-dir', '-o', type=click.Path(file_okay=False, writable=True), 
              default='./processed_files', show_default=True, help="Directory to save the output files.")
@click.option('--num-workers', '-w', type=int, default=10, show_default=True, help="Number of parallel workers.")
def process_multiple_files(input_list, output_dir, num_workers):
    """
    CLI tool to process multiple lvm-SFrame FITS files and compute night sky metrics in parallel.
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the list of files
    with open(input_list, 'r') as f:
        files = [line.strip() for line in f if line.strip()]

    try:
        # Process files in parallel
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_file, file, output_dir): file for file in files}

            for future in tqdm(futures, desc="Processing files", total=len(futures)):
                try:
                    future.result()
                    # click.echo(f"Processed: {futures[future]}")
                except Exception as e:
                    click.echo(f"Error processing {futures[future]}: {e}", err=True)
    except KeyboardInterrupt:
        click.echo("Interrupted! Cleaning up...")
        executor.shutdown(wait=False)  # Force shutdown of pool
        os._exit(1)  # Immediate termination


def process_file(input_file, output_dir):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # try:
    #     # Open the FITS file
    #     with fits.open(input_file) as hdul:
    #         # Generate output filename and table
    #         output_filename, output_table = generate_output(hdul)
    #         output_path = os.path.join(output_dir, output_filename)
    #         output_table.write(output_path, overwrite=True)

    # except Exception as e:
    #     click.echo(f"Error processing FITS file: {e}", err=True)

    with fits.open(input_file) as hdul:
        # Generate output filename and table
        output_filename, output_table = generate_output(hdul)
        output_path = os.path.join(output_dir, output_filename)
        output_table.write(output_path, overwrite=True)


def generate_output(hdul):
    """
    Generate an output filename and table
    """
    hdr = hdul[0].header
    filename = f"tbl_{hdr['EXPOSURE']:08d}.fits"

    wave = hdul['WAVE'].data
    flux = hdul['FLUX'].data
    sky = hdul['SKY'].data
    skyflux = sky + flux
    slitmap = Table.read(hdul, hdu='SLITMAP')

    wave_step = wave[1] - wave[0]

    # select needed columns
    tout = slitmap['fiberid', 'spectrographid', 'blockid', 'finblock', 'targettype',
                   'ifulabel', 'finifu','telescope', 'fibstatus', 'ra', 'dec']

    # adding id into the table
    tout_size = len(tout)
    tout.add_column([ hdr['EXPOSURE'] ] * tout_size, name='exposure', index=0)
    tout.add_column([ hdr['MJD'] ] * tout_size, name='mjd', index=1)
    tout.add_column([ hdr['TILE_ID'] ] * tout_size, name='tile_id', index=2)
    tout.add_column([ hdr['DPOS'] ] * tout_size, name='dpos', index=3)

    # adding night sky lines metrics
    for wline in WLINES:
        msk_line = ( wave >= (wline['wave'] + wline["line_region"][0]) ) & \
                   ( wave <= (wline['wave'] + wline["line_region"][1]) )

        msk_cont = ( ( wave >= wline['wave'] + wline['cont_regions'][0][0] ) & ( wave < wline['wave'] + wline['cont_regions'][0][1] ) ) | \
                   ( ( wave > wline['wave'] + wline['cont_regions'][1][0] ) & ( wave <= wline['wave'] + wline['cont_regions'][1][1] ) )

        for spec, suffix in zip([flux, skyflux], ['s', 'c']):
            cont_1d_median = np.nanmedian(spec[:, msk_cont], axis=1)
            cont_1d_mean = np.nanmean(spec[:, msk_cont], axis=1)
            cont = np.broadcast_to(cont_1d_median[:, np.newaxis], spec.shape)
            spec_nocont = spec - cont
            line_sum = np.nansum(spec_nocont[:, msk_line], axis=1) * wave_step
            line_sum_trapz = np.trapz(spec_nocont[:, msk_line], x=wave[msk_line], axis=1)
            line_mean = np.nanmean(spec_nocont[:, msk_line], axis=1)
            line_std = np.nanstd(spec_nocont[:, msk_line], axis=1)
            tout.add_column(line_sum, name=f"{suffix}_{wline['name']}_sum")
            tout.add_column(line_sum_trapz, name=f"{suffix}_{wline['name']}_sumtrapz")
            tout.add_column(line_mean, name=f"{suffix}_{wline['name']}_mean")
            tout.add_column(line_std, name=f"{suffix}_{wline['name']}_std")
            tout.add_column(cont_1d_median, name=f"{suffix}_{wline['name']}_cont_median")
            tout.add_column(cont_1d_mean, name=f"{suffix}_{wline['name']}_cont_mean")

    # adding continuum metrics
    for wcont in WCONTS:
        msk_cont = ( ( wave >= wcont['wave'] + wcont['cont_regions'][0][0] ) & ( wave < wcont['wave'] + wcont['cont_regions'][0][1] ) )

        for spec, suffix in zip([flux, skyflux], ['s', 'c']):
            cont_median = np.nanmedian(spec[:, msk_cont], axis=1)
            cont_mean = np.nanmean(spec[:, msk_cont], axis=1)
            cont_std = np.nanstd(spec[:, msk_cont], axis=1)
            tout.add_column(cont_median, name=f"{suffix}_{wcont['name']}_median")
            tout.add_column(cont_mean, name=f"{suffix}_{wcont['name']}_mean")
            tout.add_column(cont_std, name=f"{suffix}_{wcont['name']}_std")

    return filename, tout


if __name__ == '__main__':
    process_multiple_files()