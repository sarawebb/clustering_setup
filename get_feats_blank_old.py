"""get_feats_blank.py -- Input a field, year, month, day and min and max numbers of the lightcurves you want to run the feature extraction on, will return sub feature catalog. 

Usage: get_feats_blank [-h] [-v] [--debug] <field> <year> <month> <day> <min_num> <max_num>
Arguments:
    field (string)
        The DWF field name. 
    year (string)
	The year.
    month (string)
	The month. 
    day (string)
	The day of the lightcurves you want to make.
    min_num (integer)
        The min number of lightcurves to check. 
    max_num (integer)
	The max number of the lightcurves to check.

Options:
    -h  Show this screen
    -v	Show extra information [default: False]     
    --debug	Output more for debugging [default: False]
Example
    python get_feats_blank.py -v Antlia 2017 02 170204 0 10000
 """

import docopt 
import FATS
import os 
import math 
import numpy as np 

from astropy.table import Table, Column, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier

def get_feats_blank(field, year, month, day, min_num, max_num, verbose=False, debugmode=False):
	path ="/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/"+year +"/"+month+"/"+ field +"/g_band/single/lightcurves/files/"
	print(path)

	lc_with_only_zeros = []
	used_lcs = []
	filenames = []

	Autocor_length = []
	Beyond1Std = []
	CAR_sigma = []
	CAR_mean = []
	CAR_tau =[]
	Con =[]
	Eta_e = []
	LinearTrend = []
	MaxSlope =[]
	Mean = []
	Meanvariance = []
	MedianAbsDev =[]
	MedianBRP =[]
	PairSlopeTrend =[]
	PercentAmplitude =[]
	Q31 =[]
	Rcs =[]
	Skew =[]
	SlottedA_length =[]
	SmallKurtosis =[]
	Std = []
	StetsonK_AC = []
	Amplitudes = []
	VariabilityIndex = []
	pmra = []
	pmde = []
	gaia_G_RP = []
	gaia_BP_G = []
	detection_fraction = []

	test_filename = os.listdir(path)
	#print(test_filename[0:1])

	print(len(test_filename))

	for filename in test_filename[min_num:max_num]:
	    print filename
	    if filename.endswith(day):
		print filename
		#print '-------'
		#print path + i
		#print '-------'
		try:
		    mjd, mag, emag, uplim = np.loadtxt(path + filename, unpack = True, skiprows=1)
		    #print mjd
		except:
		    print filename
		#mjd = np.loadtxt(path+i, usecols =1, skiprows =1)
		sum_mag = np.sum(mag)
		#print(sum_mag)
		if sum_mag == 0:
		    lc_with_only_zeros.append(filename)
		elif sum_mag != 0:
		    used_lcs.append(filename)


	#print(flux)


		    # Remove the non detections
		    clean_mjd = []
		    clean_mag = []
		    clean_emag = []

		    for l,m,n in zip(mjd,mag, emag):
			if m != 0:
				if m < 25 and n < 0.8:
					clean_mjd.append(l)
					clean_mag.append(m)
					clean_emag.append(n)
			elif m == 0:
			    pass

		    print(clean_mjd, clean_mag, clean_emag)

		    if len(clean_mjd) > 3 :
			lc = np.array([clean_mag, clean_mjd, clean_emag])
			a=FATS.FeatureSpace(featureList=['Autocor_length', 'Beyond1Std', 'CAR_sigma', 'CAR_mean',
						     'CAR_tau', 'Con', 'Eta_e', 'LinearTrend', 'MaxSlope',
						     'Mean', 'Meanvariance', 'MedianAbsDev', 'MedianBRP',
						     'PairSlopeTrend', 'PercentAmplitude', 'Q31', 'Rcs', 'Skew',
						     'SlottedA_length', 'SmallKurtosis', 'Std',
						     'StetsonK_AC' ], Data=['magnitude','time', 'error'])
			a=a.calculateFeature(lc)
			results = a.result(method='dict')
			#print results
			filenames.append(filename)
			Autocor_length.append(results['Autocor_length'])
			Beyond1Std.append(results['Beyond1Std'])
			CAR_sigma.append(results['CAR_sigma'])
			CAR_mean.append(results['CAR_mean'])
			CAR_tau.append(results['CAR_tau'])
			Con.append(results['Con'])
			Eta_e.append(results['Eta_e'])
			LinearTrend.append(results['LinearTrend'])
			MaxSlope.append(results['MaxSlope'])
			Mean.append(results['Mean'])
			Meanvariance.append(results['Meanvariance'])
			MedianAbsDev.append(results['MedianAbsDev'])
			MedianBRP.append(results['MedianBRP'])
			PairSlopeTrend.append(results['PairSlopeTrend'])
			PercentAmplitude.append(results['PercentAmplitude'])
			Q31.append(results['Q31'])
			Rcs.append(results['Rcs'])
			Skew.append(results['Skew'])
			SlottedA_length.append(results['SlottedA_length'])
			SmallKurtosis.append(results['SmallKurtosis'])
			Std.append(results['Std'])
			StetsonK_AC.append(results['StetsonK_AC'])



			# Find Amplitude because FATS is unable to
			N = len(clean_mag)
			sorted_mag = np.sort(clean_mag)
			amp = (np.median(sorted_mag[-int(math.ceil(0.05*N)):]) - np.median(sorted_mag[0:int(math.ceil(0.05*N))])) / 2
			#print amp
			Amplitudes.append(amp)

	 # Find VariabilityIndex beacuase FATS has it removed
			N = len(clean_mag)
			clean_mag_array = np.asarray(clean_mag)
			sigma2 = np.var(clean_mag)
			VarIndex = 1/((N-1)*sigma2)*np.sum(np.power(clean_mag_array[1:]-clean_mag_array[:-1], 2))
			VariabilityIndex.append(VarIndex)

			# Find GAIA colors
			    #First get the RA and DEC seperated
			print('-------------------------------------')
			RA_h = filename[3:5]
			RA_m = filename[5:7]
			RA_s = filename[7:13]
			DEC_d = filename[13:16]
			DEC_m = filename[16:18]
			DEC_s = filename[18:24]
			print('-------------------------------------')

			coords = SkyCoord(RA_h + ' '+RA_m + ' '+RA_s +' '+str(DEC_d)+' '+ DEC_m +' '+ DEC_s + ' ', unit=(u.hourangle, u.deg))
			print(coords)
			RA = coords.ra.degree
			DEC = coords.dec.degree

			print('-------------------------------------')

			GAIA_DR2 = 'I/345'
			result = Vizier.query_region(str(RA)+' '+str(DEC), radius="0d0m2s", catalog=GAIA_DR2, cache=False)
			#print(result)
			try:
			    pmra.append(float(result[0]['pmRA']))
			except:
			    pmra.append(0)

			try:
			    pmde.append(float(result[0]['pmDE']))
			except:
			    pmde.append(0)

			try:
			    gaia_G_RP.append(float(result[0]['G-RP']))
			except:
			    gaia_G_RP.append(0)

			try:
			    gaia_BP_G.append(float(result[0]['BP-G']))
			except:
			    gaia_BP_G.append(0)

			# find fraction of detections to non detections
			detects = 0
			non_detects = 0

			for mags in mag:
			    if mags != 0:
				detects += 1
			    else:
				pass
			for uplims in uplim:
			    if uplims != 0:
				non_detects += 1
			    else:
				pass
			print('detect len: ' + str(detects))
			print('number of points: '+ str(len(mag)))


			detection_fraction.append(detects/len(mag))

		else:
			print('Not enough data points')

	print 'Not used'
	print len(lc_with_only_zeros)
	print 'Used'
	print len(used_lcs)

	feature_table = Table()
	feature_table['LC_name'] = filenames
	feature_table['Autocor_length'] = Autocor_length
	feature_table['Beyond1Std'] = Beyond1Std
	feature_table['CAR_sigma'] = CAR_sigma
	feature_table['CAR_mean'] = CAR_mean
	feature_table['CAR_tau'] = CAR_tau
	feature_table['Con'] = Con
	feature_table['Eta_e'] = Eta_e
	feature_table['LinearTrend'] = LinearTrend
	feature_table['MaxSlope'] = MaxSlope
	feature_table['Mean'] =Mean
	feature_table['Meanvariance'] =Meanvariance
	feature_table['MedianAbsDev'] =MedianAbsDev
	feature_table['MedianBRP'] =MedianBRP
	feature_table['PairSlopeTrend'] =PairSlopeTrend
	feature_table['PercentAmplitude'] =PercentAmplitude
	feature_table['Q31'] =Q31
	feature_table['Rcs'] =Rcs
	feature_table['Skew'] =Skew
	feature_table['SlottedA_length'] =SlottedA_length
	feature_table['SmallKurtosis'] =SmallKurtosis
	feature_table['Std'] =Std
	feature_table['StetsonK_AC'] =StetsonK_AC
	feature_table['Amplitudes'] = Amplitudes
	feature_table['VariabilityIndex'] = VariabilityIndex
	feature_table['pmra'] = pmra
	feature_table['pmde'] = pmde
	feature_table['gaia_G_RP'] = gaia_G_RP
	feature_table['gaia_BP_G'] = gaia_BP_G
	feature_table['detection_fraction'] = detection_fraction

	print feature_table
	output = '/fred/oz100/LC_SUBSETS/feats_lists/'+year+'_'+month+'_'+field+'_'+day+'_feat_'+str(max_num)+'.ascii'
	feature_table.write(output, format = 'ascii', overwrite = True)
	return featue_table

if __name__ == "__main__":
	arguments = docopt.docopt(__doc__)
	field = arguments['<field>']
	year = arguments['<year>']
	month = arguments['<month>']
	day = arguments['<day>']
	min_num = arguments['<min_num>']
	max_num = arguments['<max_num>']

	verbose = arguments['--verbose']
	debugmode = arguments['--debug']
	print('--------')
	print('---TESTING----')
	print(field, year, month, day)
	get_feats_blank(field, year, month, day, min_num, max_num, verbose=verbose, debugmode=debugmode)
