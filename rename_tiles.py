#!/usr/bin/env python3

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy_healpix import HEALPix
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import sys
import argparse
import logging
import os


logging.basicConfig(level=logging.INFO)


"""
This scripts renames a given tile. For a given tile fits data, it will read the
position, Stokes, beam size, frequency information to derive/assign the needed pixel ID, 
Stokes parameter and central frequency to the output file. The output file will take 
the form:

    prefix_cenfreq_resolution_RADEC_TileID_Stokes_version.fits

where:
    1. prefix is user given parameter (the prefix to use for the outname e.g PSM_pilot1
    2. cenfreq is the central frequency derived within the code.
    3. resolution is beam resolution extracted within code.
    4. RADEC are RA and DEC of the tile in hours-minutes and degrees-minutes, derived
       within code.
    5. TileID is the HPX ID for a specific sky region for a given NSIDE. Also derived
       within the code.
    6. Stokes paramater extracted within the code.
    7. Version number to be specified by a user, e.g. v1 for version 1 of the data products.
      
"""



def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", dest="fitsimage", help="Image file", required=True)
    parser.add_argument("-pf", dest="prefix", help="Output prefix", required=True)
    parser.add_argument("-v", dest="version", help="Output version", required=True)
    parser.add_argument("-mfs", dest="do_mfs", help="Input image is an mfs, so add 't0'",
    action='store_true')

    args = parser.parse_args(argv)

    return args
    
    
def reference_header(naxis, cdelt):
    
    """
       Reference header centred at (0, 0) in a healpix grid. 
       This is important as it allows us to properly determine 
       correct pixel central pixel anywhere within the grid.
       
       NB: We use this header to convert the crpix1/2 in the header to tile ID,
       then degrees. 
       
       cdelt : the pixel size of the image in the grid. Must be the 
               same as the one used for tiling.
       naxis : number of pixels within each axis.
       
    """
   
    
    hdr  = "SIMPLE  =                    T / file does conform to FITS standard \n"
    hdr += "BITPIX  =                  -32 / number of bits per data pixel \n"
    hdr += "NAXIS   =                    2 / number of data axes  \n"
    hdr += "NAXIS1  =                %d / length of data axis 1 \n" %naxis
    hdr += "NAXIS2  =                %d / length of data axis 2 \n" %naxis
    hdr += "EXTEND  =                    F / No FITS extensions are present \n"
    hdr += "CRPIX1  =             %r / Coordinate reference pixel \n"%((naxis/2.0) + 0.0)
    hdr += "CRPIX2  =             %r / Coordinate reference pixel \n"%((naxis/2.0) + 0.5)
    hdr +=  "PC1_1   =           0.70710677 / Transformation matrix element \n"
    hdr +=  "PC1_2   =           0.70710677 / Transformation matrix element \n"
    hdr +=  "PC2_1   =           -0.70710677 / Transformation matrix element \n"
    hdr +=  "PC2_2   =           0.70710677 / Transformation matrix element \n"
    hdr += "CDELT1  =            -%r  / [deg] Coordinate increment \n" %cdelt
    hdr += "CDELT2  =             %r  / [deg] Coordinate increment \n" %cdelt 
    hdr += "CTYPE1  = 'RA---HPX'           / Right ascension in an HPX projection \n"
    hdr += "CTYPE2  = 'DEC--HPX'           / Declination in an HPX projection \n"
    hdr += "CRVAL1  =                   0. / [deg] Right ascension at the reference point \n"
    hdr += "CRVAL2  =                   0. / [deg] Declination at the reference point \n"
    hdr += "PV2_1   =                    4 / HPX H parameter (longitude) \n"
    hdr += "PV2_2   =                    3 / HPX K parameter  (latitude) \n"

    return hdr


def name(fitsimage, prefix, version="v1", mfs=False):
    """Setting up the name to be used for tiles. The script reads
    the bmaj and stokes from the fits header. The rest of the parameters are
    flexible to change.

    fitsimage: tile image
    prefix   : prefix to use. E.g. PSM for full survey,
               PSM_pilot1 for POSSUM pilot 1
               PSM_pilot2 for POSSUM pilot 2
    tileID   : tile pixel (Healpix pixel)

    version  : version of the output product. Version 1 is v1, version is v2,
               and so forth.

    """

    logging.info(f"Reading {fitsimage} header")

    hdr = fits.getheader(fitsimage)

    # get bmaj.
    bmaj = round(hdr["BMAJ"] * 3600.0)
    bmaj =  '%dasec'%bmaj

    # extract stokes parameter. It can be in either the 3rd or fourth axis.

    if hdr["CTYPE3"] == "STOKES":
        stokes = hdr["CRVAL3"]
        # if Stokes is axis 3, then frequency is axis 4. 
        freq0 = hdr['CRVAL4']
        dfreq = hdr['CDELT4']
        N = hdr['naxis4']
        if N > 1:
            cenfreq = round((freq0 + (freq0 + N * dfreq))/(2.0 * 1e6))
        else:
            cenfreq = round(freq0/1e6)
        

    elif hdr["CTYPE4"] == "STOKES":
        stokes = hdr["CRVAL4"]
        # if Stokes is axis 4, then frequency is axis 3. If we have >4 axis, the script will fail.
        freq0 = hdr['CRVAL3']
        dfreq = hdr['CDELT3']
        N = hdr['naxis3']
        if N > 1:
            cenfreq = round((freq0 + (freq0 + N * dfreq))/(2.0 * 1e6))
        else:
            cenfreq = round(freq0/1e6)

    else:
        sys.exit(">>> Cannot find Stokes axis on the 3rd/4th axis")
        
    cenfreq = '%dMHz'%cenfreq

    # stokes I=1, Q=2, U=3 and 4=V
    if int(stokes) == 1:
        stokesid = "i"

    elif int(stokes) == 2:
        stokesid = "q"

    elif int(stokes) == 3:
        stokesid = "u"

    elif int(stokes) == 4:
        stokesid = "v"

    logging.info("Define healpix grid for nside 32")
    # define the healpix grid
    hp = HEALPix(nside=32, order="ring", frame="icrs")

    # read the image crpix1 and crpix2 to determine the tile ID, and coordinates in degrees.
    naxis = hdr['naxis1']
    cdelt = abs(hdr['cdelt1'])
    hpx_ref_hdr = reference_header(naxis=naxis, cdelt=cdelt)
    hpx_ref_hdr = fits.Header.fromstring("""%s"""%hpx_ref_hdr, sep='\n')
    hpx_ref_wcs = WCS(hpx_ref_hdr)
   
    crpix1 = hdr['crpix1']
    crpix2 = hdr['crpix2']
    crval1, crval2 = hpx_ref_wcs.wcs_pix2world(-crpix1, -crpix2 , 0)    
    tileID = hp.lonlat_to_healpix(crval1 * u.deg, crval2 * u.deg, return_offsets=False)
    tileID = tileID - 1 #shifts by 1. 
    
    # extract the RA and DEC for a specific pixel
    center = hp.healpix_to_lonlat(tileID) * u.deg
    RA, DEC = center.value

    logging.info(f"Derived RA is {RA} degrees and DEC is {DEC} degrees")
    c = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree, frame="icrs")

    h, hm, hs = c.ra.hms
    hmhs = "%s" % round(hm + (hs / 60.0))
    hmhs = hmhs.zfill(2)
    if h < 0.0:
        hm = f"{int(h):03d}{hmhs}"
    else:
        hm = f"{int(h):02d}{hmhs}"
    

    d, dm, ds = c.dec.dms
    # if dec is in the southern sky leave as is. If northen add a +.
    dmds = "%s" % round(abs(dm) + (abs(ds) / 60.0))
    dmds = dmds.zfill(2)
    dm = f"{int(d):02d}{dmds}"
    if (c.dec < 0) and (dm[0] != '-'):
        dm = '-'+dm
    if c.dec > 0:
        dm = "+" + dm

    RADEC = "%s%s" % (hm, dm)


    if mfs:
        outname = (
            prefix
            + "_%s" % cenfreq
            + "_%s" % bmaj
            + "_%s" % (RADEC)
            + "_%s" % (tileID)
            + "_t0"
            + "_%s" % stokesid
            + "_%s" % version
            + ".fits"
         )

    else:
        outname = (
            prefix
            + "_%s" % cenfreq
            + "_%s" % bmaj
            + "_%s" % (RADEC)
            + "_%s" % (tileID)
            + "_%s" % stokesid
            + "_%s" % version
            + ".fits"
         )
    print(outname, end="")
    
    os.system('mv %s %s'%(fitsimage, outname))



def main(argv):
    """ Renames tile images.

    Usage:
        python renaming_tiles -i <image_name> -pf <output_prefix> 
        -id <tile_ID> -v <version number>

    Returns a new name for a tile

    """

    args = parse_args(argv)

    name(
        fitsimage=args.fitsimage,
        prefix=args.prefix,
        version=args.version,
        mfs = args.do_mfs
    )


if __name__ == "__main__":
    argv = sys.argv[1:]
    main(argv)
