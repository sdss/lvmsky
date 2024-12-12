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

    try:
        # Open the FITS file
        with fits.open(input_file) as hdul:
            # Generate output filename and table
            output_filename, output_table = generate_output(hdul)
            output_path = os.path.join(output_dir, output_filename)
            output_table.write(output_path, overwrite=True)

    except Exception as e:
        click.echo(f"Error processing FITS file: {e}", err=True)


def generate_output(hdul):
    """
    Generate an output filename and table
    """
    hdr = hdul[0].header
    filename = f"tbl_{hdr['EXPOSURE']:08d}.fits"

    wave = hdul['WAVE'].data
    flux = hdul['FLUX'].data
    sky = hdul['SKY'].data
    slitmap = Table.read(hdul, hdu='SLITMAP')

    spec = sky + flux

    w = 5
    msk_cont = ((wave >= 5574 - w) & (wave < 5574)) | ((wave > 5580) & (wave <= 5580 + w)) 
    msk_line = (wave >= 5574) & (wave <= 5580)

    spec_2d = spec[:, msk_line]

    cont_1d = np.nanmedian(spec[:, msk_cont], axis=1)
    cont_2d = np.broadcast_to(cont_1d[:, np.newaxis], spec_2d.shape)
    line_2d = spec_2d - cont_2d
    wave_step = wave[1] - wave[0]
    flux_2d_sum = np.nansum(line_2d, axis=1) * wave_step
    flux_2d_trapz = np.trapz(line_2d, x=wave[msk_line], axis=1)

    # import ipdb; ipdb.set_trace()
    tout = slitmap['fiberid', 'spectrographid', 'blockid', 'finblock', 'targettype',
                   'ifulabel', 'finifu','telescope', 'fibstatus', 'ra', 'dec']

    tout['median'] = np.nanmedian(spec, axis=1)
    tout['cont_OI5577'] = cont_1d
    tout['flux_OI5577_sum'] = flux_2d_sum
    tout['flux_OI5577_trapz'] = flux_2d_trapz

    tout_size = len(tout)
    tout.add_column([ hdr['EXPOSURE'] ] * tout_size, name='exposure', index=0)
    tout.add_column([ hdr['MJD'] ] * tout_size, name='mjd', index=1)
    tout.add_column([ hdr['TILE_ID'] ] * tout_size, name='tile_id', index=2)
    tout.add_column([ hdr['DPOS'] ] * tout_size, name='dpos', index=3)

    return filename, tout


if __name__ == '__main__':
    process_multiple_files()