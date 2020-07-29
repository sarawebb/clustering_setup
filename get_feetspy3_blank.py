import docopt 
import feets
import os 
import math 
import numpy as np 
import argparse 
from astropy.table import Table, Column, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier
import sys
import pandas
fft = False
np.seterr(divide='ignore', invalid='ignore')

def get_feats_blank(field, year, month, day, min_num, max_num):
	path ="/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/"+year +"/"+month+"/"+ field +"/g_band/single/lightcurves/files/"
	#print(ath)

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
	hl_ratio = []
	test_filename = os.listdir(path)
	#print(test_filename[0:1])
	amp1_val = []
	amp_2_1_ratio = []
	amp_3_1_ratio = []
	phase_2_1_ratio = []
	phase_3_1_ratio = []
	
	#print(len(test_filename))

	for filename in test_filename[min_num:max_num]:
		if filename.endswith(day):
			try:
				mjd, mag, emag, uplim = np.loadtxt(path + filename, unpack = True, skiprows=1)
		    #print mjd
			except:
                    		print('FILE EMPTY)')
		#mjd = np.loadtxt(path+i, usecols =1, skiprows =1)
			sum_mag = np.sum(mag)
			print(sum_mag)
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

		    #print(clean_mjd, clean_mag, clean_emag)
			print(clean_mjd)
			if len(clean_mjd) > 3 :
				print('in cleaned loop')
				lc = np.array([clean_mag, clean_mjd, clean_emag])
			
				fs=feets.FeatureSpace(only=['Autocor_length', 'Beyond1Std', 'CAR_sigma', 'CAR_mean',
						     'CAR_tau', 'Con', 'LinearTrend', 'MaxSlope',
						     'Mean', 'Meanvariance', 'MedianAbsDev', 'MedianBRP',
						     'PairSlopeTrend', 'PercentAmplitude', 'Q31', 'Rcs', 'Skew',
						     'SlottedA_length', 'SmallKurtosis', 'Std',
						     'StetsonK_AC' ])
				features, values = fs.extract(*lc)
				results = dict(zip(features,values)) 
				#except:
				#	print filename 
				print('after feets')
				filenames.append(filename)
				Autocor_length.append(results['Autocor_length'])
				Beyond1Std.append(results['Beyond1Std'])
				CAR_sigma.append(results['CAR_sigma'])
				CAR_mean.append(results['CAR_mean'])
				CAR_tau.append(results['CAR_tau'])
				Con.append(results['Con'])
				#Eta_e.append(results['Eta_e'])
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
				N = len(clean_mag)
				sorted_mag = np.sort(clean_mag)
				amp = (np.median(sorted_mag[-int(math.ceil(0.05*N)):]) - np.median(sorted_mag[0:int(math.ceil(0.05*N))])) / 2
				Amplitudes.append(amp)
				N = len(clean_mag)
				clean_mag_array = np.asarray(clean_mag)
				sigma2 = np.var(clean_mag)
				VarIndex = 1/((N-1)*sigma2)*np.sum(np.power(clean_mag_array[1:]-clean_mag_array[:-1], 2))
				VariabilityIndex.append(VarIndex)
				non_detects = []
				for lim in uplim:
					if lim != 0:
						non_detects.append(lim)
						len_mags = len(clean_mag)
						len_uplim = len(non_detects)
						if len_uplim == 0: 
							detection_fraction.append(1)
						elif len_uplim != 0: 
							detection_fraction.append(len_mags/len_uplim)
				fft = np.fft.rfft(clean_mag)
				amps = np.sqrt(fft.real**2+fft.imag**2)
				amp1 = amps[0]
				amp1_val.append(amp1)
				amp2 = amps[1]
				amp3 = amps[2]
				amp_2_1 = amp2/amp1
				amp_3_1 = amp3/amp1 
				amp_2_1_ratio.append(amp_2_1)
				amp_3_1_ratio.append(amp_3_1)
				phases = np.arctan2(fft.imag, fft.real)
				phase1 = phases[0]
				phase2 = phases[1]
				phase3 = phases[2]
				phase_2_1 = phase2/phase1
				phase_3_1 = phase3/phase1  
				phase_2_1_ratio.append(phase_2_1)
				phase_3_1_ratio.append(phase_3_1)

			else:
				print('Not enough data points')

	#print 'Not used'
	#print len(lc_with_only_zeros)
	#print 'Used'
	#print len(used_lcs)
	feature_table = Table()
	feature_table['LC_name'] = filenames
	feature_table['Autocor_length'] = Autocor_length
	feature_table['Beyond1Std'] = Beyond1Std
	feature_table['CAR_sigma'] = CAR_sigma
	feature_table['CAR_mean'] = CAR_mean
	feature_table['CAR_tau'] = CAR_tau
	feature_table['Con'] = Con
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
	#feature_table['DetectionFraction'] = detection_fraction
	feature_table['amp1'] = amp1_val
	feature_table['amp_2_1_ratio'] = amp_2_1_ratio 
	feature_table['amp_3_1_ratio'] = amp_3_1_ratio
	feature_table['phase_2_1_ratio'] = phase_2_1_ratio
	feature_table['phase_3_1_ratio'] = phase_3_1_ratio
	output = '/home/swebb/oz100/LC_CLUSTERS/feature_lists/'+year+'_'+month+'_'+field+'_'+day+'/'+year+'_'+month+'_'+field+'_'+day+'_feat_'+str(max_num)+'.csv'
	print(output)
	df = feature_table.to_pandas()
	df.to_csv(output)
	#feature_table.write(output, format = 'ascii', overwrite = True)
	return feature_table

if len(sys.argv)>1:
	field = str(sys.argv[1])
	year = str(sys.argv[2])
	month = str(sys.argv[3])
	day = str(sys.argv[4])
	min_num = int(sys.argv[5])
	max_num = int(sys.argv[6])

get_feats_blank(field, year, month, day, min_num, max_num)
