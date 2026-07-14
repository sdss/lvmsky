def fit_chunk_worker(args, progress_queue):
    """Worker function: builds its own SkyDecomp and fits a chunk of spectra."""
    chunk_idxs, flx_sci_chunk, flx_sky1_chunk, flx_sky2_chunk, flx_ivar_chunk, wave, lsf_sigma, base_dir = args

    from sky_decomp.fit import SkyDecomp
    local_decomposer = SkyDecomp(wave, lsf_sigma=lsf_sigma, base_dir=base_dir)

    local_sci, local_sky1, local_sky2 = [], [], []
    for i in range(len(chunk_idxs)):
        local_sci.append(local_decomposer.fit(flx_sci_chunk[i],  flx_ivar_chunk[i], verbose=False, n_lsf_refits=3))
        local_sky1.append(local_decomposer.fit(flx_sky1_chunk[i], flx_ivar_chunk[i], verbose=False, n_lsf_refits=3))
        local_sky2.append(local_decomposer.fit(flx_sky2_chunk[i], flx_ivar_chunk[i], verbose=False, n_lsf_refits=3))
        progress_queue.put(1)  # Report one spectrum completed

    return chunk_idxs, local_sci, local_sky1, local_sky2