#!/usr/bin/env python
import numpy
import pandas as pd
from astropy_healpix import HEALPix
from astropy.coordinates import Angle
from astropy import units as u
import pyregion
from astropy.io import fits
from astropy.wcs import WCS
import glob
import json
import sys
import argparse


"""

    This script generates tile footprints containing HealPix ID, CRPIX1/2 and CRVAL1/2. 
    Inputs:
	1. Footprint region files. Supports only .txt or reg. The region file name
           should be something like name..._ID.reg/txt, where name can be anything the
           user chose to name the footprints, and ID is the SB ID e.g. '2156-54' (corresponds
           sb10040). 
        2. Tile setting. HealPix Nside, image CDELT. See input json file 'AllSky-TileConfig.json.

    Output:
        1. A xlsx and csv file containing a list tile IDs and the corresponding SBs. This file is used 
           to identify HealPix pixels that finds contribution from multiple SBs. 
        2. A xlsv and csv file for each SB containg pixes parameters (CRPIX) used for generating tiles.    
     
Author: Lerato

"""


def hms2deg(hms): 
    """ Converts unit from hms to degrees"""
    
    h, m, s = hms.split(':')
    a = Angle('%sh%sm%s'%(h, m, s)) 
    return a.deg


def dms2deg(dms): 
    
    """Converts unit from dms to degrees"""
    
    d, m, s = dms.split(':')
    a = Angle('%sd%sm%s'%(d, m, s)) 
    return a.deg


def get_deg(hms_array, dms_array):
    
    """Converts unit of an array from hms to degrees"""
    
    conversion_hms = []
    conversion_dms = []
    for (hms, dms) in zip(hms_array, dms_array):
        conversion_hms.append(hms2deg(hms))
        conversion_dms.append(dms2deg(dms))
    return conversion_hms, conversion_dms
        

def points_within_circle(x0, y0, radius, num_points=4):
    
    """
       Samples the primary beams to determine the tiles.
       
       x0, y0    :  beam center in degrees.
       radius    :  the beam radius.
       num_points: is a number of sample points within a circular beam.
       
    """
    
    angle = numpy.linspace(0, 360, num_points)
    angle = numpy.deg2rad(angle)
    
    if len(x0) and len(y0) > 1:
        x_corners = []
        y_corners = []
        for (x_0, y_0) in zip(x0, y0):
            x_corners.append( radius * numpy.cos(angle) + x_0 )
            y_corners.append( radius * numpy.sin(angle) + y_0 )
    else:
        x_corners =  radius * numpy.cos(angle) + x0
        y_corners =  radius * numpy.sin(angle) + y0
            
    return x_corners, y_corners


def get_HealPix_Tiles(ra_deg, dec_deg):
    
    """Determines healpix tiles for the beam sampled points"""
    
    SB_index = hp.lonlat_to_healpix(ra_deg * u.deg, dec_deg * u.deg, return_offsets=False)
    SB_index_unique = numpy.unique(SB_index) 
    
    return SB_index_unique


def generate_DS9_polygons(healpix_pixel, nside, outname_prefix ):

    """
       Option function: for the extracted tiles, it generates the 
       ds9 region files.
    """  

    regions = []
    centers = []
    for pixel in healpix_pixel:
        
        corner = hp.boundaries_lonlat(pixel, step=1) * u.deg
        RA, DEC = corner.value
        RA = RA[0]
        DEC = DEC[0]
        polygon_string = 'polygon(%f, %f, %f, %f, %f, %f, %f, %f)'%(RA[0], DEC[0], 
              RA[1], DEC[1], RA[2], DEC[2], RA[3], DEC[3])        
        regions.append(polygon_string)
        
        center = hp.healpix_to_lonlat(pixel) * u.deg
        center_RA, center_DEC = center.value
        
        circle_string = 'circle(%f, %f, %f)'%(center_RA, center_DEC, 0.1)
        centers.append(circle_string)
    first_line = "#Region file format: DS9 version 4.1 \n"
    second_line = 'global color=black dashlist=8 3 width=2 font="helvetica 10 normal roman" \
    select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n'
    third_line = 'fk5 \n'
    
    SB = outname_prefix
    region_file = open('%s-boundary-%d.reg'%(SB, nside), 'w')

    region_file.write(first_line)
    region_file.write(second_line)
    region_file.write(third_line)

    for region in regions:
        region_file.write(region +' \n')
    region_file.close()

    center_file = open('%s-center-%d.reg'%(SB, nside), 'w')
    center_file.write(first_line)
    center_file.write(second_line)
    center_file.write(third_line)

    for center in centers:
        center_file.write(center +' \n')
    center_file.close()
    
    

def reference_header(naxis, cdelt):
    
    """
       Reference header centred at (0, 0) in a healpix grid. 
       This is important as it allows us to properly determine 
       correct pixel central pixel anywhere within the grid.
       
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
    hdr += "CRPIX2  =             %r / Coordinate reference pixel \n"%((naxis/2.0) + 0.0)
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
    


def HPX_in_degrees(HPX, HPX_wcs):
    
    """
        Determines HealPix pixel 
        reference center coordinates.     
    """
    
    hpx_ra, hpx_dec = [], []
    crpix_ra, crpix_dec = [], []

    for SBs_pixels in HPX:
    
        hpx2lon = hp.healpix_to_lonlat(SBs_pixels) * u.deg
        HPX_RA, HPX_DEC = hpx2lon.value
        crpix_RA, crpix_DEC = HPX_wcs.wcs_world2pix(HPX_RA, HPX_DEC, 1)
        
        crpix_ra.append(-crpix_RA)
        crpix_dec.append(-crpix_DEC)
        hpx_ra.append(HPX_RA)
        hpx_dec.append(HPX_DEC)
    
    return crpix_ra, crpix_dec, hpx_ra, hpx_dec



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser("Generate files defining Healpix tiles.")
    parser.add_argument('-j', dest='json', help='The Json file.')
    args = parser.parse_args()
    config_json_file = args.json

    # read settings from the config json file.
    with open(config_json_file, "r") as read_file:
        tile_config = json.load(read_file)

    Nside           = tile_config['HPX_nside']
    beam_radius     = tile_config['beam_radius']
    beam_sample_points = tile_config['beam_sample_points']
    Naxis = tile_config['tile_naxis']
    Cdelt = tile_config['tile_cdelt']
    number_of_beams = tile_config['number_of_beams']
    
    footprint_files  = tile_config['footprint_files']
    path_footprints  = tile_config['path_footprints']
    run_all_footprints  = tile_config['run_all_footprints']
        
    outfile_json_prefix = tile_config['outfile_prefix']
    generate_ds9regions = tile_config['generate_ds9regions']
    
    # globalizing Healpix definition.
    global hp
    
    hp = HEALPix(nside=Nside, order='ring', frame='icrs')
    
    
    print('Number of pixels for nside %d is %d. '%(Nside, hp.npix))
    print('HealPix pixel resolution for nside %d is %s.' %(Nside, hp.pixel_resolution))
    

    if run_all_footprints:
        
        footprints = glob.glob(path_footprints + '*')
        footprints_ID = [regions.split('_')[-1].split('.')[0] for regions in footprints]
        extension = footprints[0].split('.')[-1]

    else:
        
        if len(footprint_files) > 0:
            footprints = [path_footprints + foot for foot in footprint_files]
            footprints_ID =  [regions.split('_')[-1].split('.')[0] for regions in footprints]
            extension = footprints[0].split('.')[-1]
            
        if len(footprint_files) == 0:
            sys.exit('No footprint file provided. Provide atleast one footprint name inside path_footprints.')
    
                
    HPX_ID_Pixels = {}
    HPX_Pixels = []
    HPX_ID = []

    for foot, footprint in enumerate(footprints):

        if extension == 'reg':
            footprint_region = pyregion.open(footprint)
            x_deg = [r.coord_list[0] for r in footprint_region]
            y_deg = [r.coord_list[1] for r in footprint_region]
        
        if extension == 'txt':
            footprint_region = pd.read_csv(footprint, header=None,  sep=")")
            x_deg = [ footprint_region[1][i].split(',')[0] for i in range(number_of_beams)]
            y_deg = [ footprint_region[1][i].split(',')[1] for i in range(number_of_beams)]
            x_deg, y_deg = get_deg(x_deg, y_deg)
            

        beam_x_corner_sample, beam_y_corner_sample = points_within_circle(
                  x_deg, y_deg, beam_radius, num_points=beam_sample_points)
        
        HealPixels = get_HealPix_Tiles(ra_deg=beam_x_corner_sample, dec_deg=beam_y_corner_sample)
     
        
        HPX_Pixels.append(HealPixels)
        HPX_ID.append(footprints_ID[foot])
        HPX_ID_Pixels.update({'%s'%footprints_ID[foot]: HealPixels})
        
    # read the reference header to estimate pixel centers in degrees, J2000. 
    HPX_hdr = reference_header(naxis=Naxis, cdelt=Cdelt)
    HPX_hdr = fits.Header.fromstring("""%s"""%HPX_hdr, sep='\n')
    HPX_wcs = WCS(HPX_hdr)
    
    crpix_ra, crpix_dec, hpx_ra, hpx_dec = HPX_in_degrees(HPX_Pixels, HPX_wcs)

    for i in range(len(HPX_ID)):        
    
        print(footprints_ID[i])
        tile_output_filename = outfile_json_prefix + '-TileConfig-%s.xlsx'%footprints_ID[i]
        csv_tile_output =  outfile_json_prefix + '-TileConfig-%s.csv'%footprints_ID[i]
        
        with pd.ExcelWriter(tile_output_filename) as writer:
        
            data = {'PIXELS': HPX_Pixels[i].tolist() ,
                'CRPIX_RA' : crpix_ra[i].tolist(),
                'CRPIX_DEC' : crpix_dec[i].tolist(),
                'CRVAL_RA [deg]'  : hpx_ra[i].tolist(),
                'CRVAL_DEC [deg]' : hpx_dec[i].tolist()}
            
            df = pd.DataFrame(data, columns = ['PIXELS', 'CRPIX_RA', 'CRPIX_DEC', 
                                              'CRVAL_RA [deg]', 'CRVAL_DEC [deg]'])
            df.to_excel(writer, sheet_name='%s'%footprints_ID[i])
            df.to_csv(csv_tile_output, index=False)
    
            column_width = 15
            workbook  = writer.book
            worksheet = writer.sheets['%s'%footprints_ID[i]]
            cell_format = workbook.add_format({'font_name':'Calibri', 'font_size':10})
            worksheet.set_column(0, 5, column_width, cell_format)
            header_format = workbook.add_format({
                    'bold': True,
                    'text_wrap': True,
                    'valign': 'top',
                    'fg_color': '#D7E4BC',
                    'border': 1})
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num + 1, value, header_format)
     

    HPX_PIXELS = numpy.unique(numpy.hstack(HPX_Pixels)) #all HPX tiles.
 
    SBsID = []
    SBs_HPX = []
    
   
    for hpxs in HPX_PIXELS:
        for SBid in HPX_ID:
            if hpxs in set(HPX_ID_Pixels[SBid]):
                SBs_HPX.append((hpxs))
                SBsID.append(SBid)    
            
    repeat_tiles_outfile = outfile_json_prefix + '-TileRepeat.xlsx'
    csv_repeat_tiles = outfile_json_prefix + '-TileRepeat.csv'
    
    with pd.ExcelWriter(repeat_tiles_outfile) as writer:
        data = []
        j = 0 
        for hpx in numpy.unique(SBs_HPX):
        
            count = SBs_HPX.count(hpx)
            if count is 1:
                j+=1
                
            if count is 2:
                ind = numpy.where(hpx == SBs_HPX)[0].tolist()
                SBs_data = numpy.asarray(SBsID)[ind]
                SBs_temp = SBs_data.tolist() 
                data.append([hpx, SBs_temp[0], SBs_temp[1]])
        
            if count is 3:
                ind = numpy.where(hpx == SBs_HPX)[0].tolist()
                SBs_data = numpy.asarray(SBsID)[ind]
                SBs_temp = SBs_data.tolist()                         
                data.append([hpx, SBs_temp[0], SBs_temp[1], SBs_temp[2]])
        
            if count is 4:
                ind = numpy.where(hpx == SBs_HPX)[0].tolist()
                SBs_data = numpy.asarray(SBsID)[ind]
                SBs_temp = SBs_data.tolist()                                       
                data.append([hpx, SBs_temp[0], SBs_temp[1], SBs_temp[2], SBs_temp[3]])
            
        df2 = pd.DataFrame(data, columns = ['PIXEL', 'SB1', 'SB2', 'SB3', 'SB4'])
        df2.to_excel(writer, sheet_name='Tiles')
        df2.to_csv(csv_repeat_tiles, index=False)
    
        column_width = 12
        workbook  = writer.book
        worksheet = writer.sheets['Tiles']
        cell_format = workbook.add_format({'font_name':'Calibri', 'font_size':10})
        worksheet.set_column(0, 5, column_width, cell_format)
        header_format = workbook.add_format({
               'bold': True,
               'text_wrap': True,
               'valign': 'top',
               'fg_color': '#D7E4BC',
               'border': 1})
        for col_num, value in enumerate(df2.columns.values):
            worksheet.write(0, col_num + 1, value, header_format)
    
    if generate_ds9regions:
        generate_DS9_polygons(healpix_pixel=HPX_PIXELS, nside=Nside, outname_prefix=outfile_json_prefix)

