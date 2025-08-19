#!/usr/bin/env python3
"""
Examples for eso_skymodel_eval.py.
"""

from eso_skymodel_eval import evaluate_skymodel, evaluate_batch

# Single evaluation examples

def example_single_table():
    tab = evaluate_skymodel(ra=121.75, dec=-29.7, obstime="2012-07-17T21:12:14")
    print(tab)


def example_single_arrays():
    arrays, meta = evaluate_skymodel(
        ra=121.75,
        dec=-29.7,
        obstime="2012-07-17T21:12:14",
        return_format="arrays",
    )
    print(list(arrays.keys()), meta)


def example_single_save():
    evaluate_skymodel(
        ra=121.75,
        dec=-29.7,
        obstime="2012-07-17T21:12:14",
        save_fits="sky_single.fits",
        hdu_extname="SKY",
    )
    print("Saved sky_single.fits")


# Batch evaluation examples

def example_batch_in_memory():
    records = [
        (121.75, -29.7, "2012-07-17T21:12:14"),
        (200.1, 10.5, "2019-05-20T08:01:02.123"),
    ]
    results = evaluate_batch(records, ncpu=1)
    print([len(r.table) for r in results])


def example_batch_save():
    records = [
        (121.75, -29.7, "2012-07-17T21:12:14"),
        (200.1, 10.5, "2019-05-20T08:01:02.123"),
    ]
    evaluate_batch(records, ncpu=1, save_fits="sky_batch.fits")
    print("Saved sky_batch.fits")


if __name__ == "__main__":
    example_single_table()
    example_single_arrays()
    example_single_save()
    example_batch_in_memory()
    example_batch_save()
