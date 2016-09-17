import numpy as np
import sys
import re
import itertools
import subprocess as sub

import matplotlib.pyplot as plt


import pandas as pd

allele_pattern = re.compile("L(\d+)\.A(\d+)")



def parse_haplos(filename):
	possible_haplos_array = []
	with open(filename,'r') as fp:
		lines = fp.readlines()
		for idx,line in enumerate(lines):
			possible_haplos_array.append([])
			haps_descr = line.strip().split(' ')
			for item in haps_descr:
				# assuming the input file is formatted ideally ...
				res = allele_pattern.match(item)
				locus,allele = map(int,res.groups())
				# populate haplotypes array ...
				possible_haplos_array[idx].append(allele)
	return possible_haplos_array





def parse_haplolimit_output(outlines):
	#
	index_pattern = re.compile("H(\d+)")
	haplo_pattern = re.compile("L(\d+)A(\d+)")
	# skip 7 lines and then go on
	skip_line_number = 7
	# start parsing the thing ...
	haplo_out = []
	for line in outlines[skip_line_number:]:
		index,fmin,fmax,haplo = line.strip().split()
		fmin,fmax = float(fmin),float(fmax)
		index, = index_pattern.match(index).groups()
		index = int(index)
		haplo = haplo.strip(':').split(':')
		haplo_parsed = [haplo_pattern.match(_hap).groups() for _hap in haplo]
		recovered_haplo = [int(allele) if ((idx+1)==int(lid)) else None for idx,(lid,allele) in enumerate(haplo_parsed) ]
		if None not in recovered_haplo:
			entry = (index,fmin,fmax,recovered_haplo)
		else:
			print "could not parse the haplo output file ... exit!"
			sys.exit(1)
		# appending parsed info ...
		haplo_out.append(entry)
	return haplo_out




# phaps = parse_haplos("./input/possibleHaps.txt")

# # freqs = parse_pop("good_test.dat")
# freqs = parse_pop("./input/lociInfo.txt")
cols = "haplo_idx haplo fmin fmax scan_idx fscan ratio_prod".split(' ')
dat = pd.read_csv('violin.dat',header=None,names=cols,sep=' ')
# dat = pd.read_csv('violin.dat',sep=' ')


dat_grp = dat.groupby('haplo')










plt.clf()
# let's draw stuff ...
fig,ax = plt.subplots(figsize=(9,10))





for scan_hap_idx,haplo in enumerate(dat['haplo'].unique()):
# # scan haplotype number XXX :
# for scan_hap_idx in range(21):
	# scan_hap_freq = [0.272,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.381]
	#
	#
	haplo_dat = dat_grp.get_group(haplo)
	fmin = haplo_dat['fmin'].mean()
	fmax = haplo_dat['fmax'].mean()
	#
	distro = list(haplo_dat[['fscan','ratio_prod']].itertuples(index=False))

	# now plot violin-like plot ...
	xarr = haplo_dat['fscan']
	max_ratio = haplo_dat['ratio_prod'].max()
	if max_ratio == 0.0:
		yup = [ (21-scan_hap_idx)+0.3 for yy in haplo_dat['ratio_prod']]
		ydown = [ (21-scan_hap_idx)-0.3 for yy in haplo_dat['ratio_prod']]
	else:
		yup = [ (21-scan_hap_idx)+0.4*yy/max_ratio for yy in haplo_dat['ratio_prod']]
		ydown = [ (21-scan_hap_idx)-0.4*yy/max_ratio for yy in haplo_dat['ratio_prod']]
	# ax.plot(xarr,y)
	ax.fill_between(xarr, yup, ydown, facecolor="blue", alpha=0.75)

ax.set_ylim(0,22)


plt.show()










def plot_haplo_range(ax,fmin_arr,fmax_arr,fillcolor='red'):
	# def nonlinear(x):
	# 	return np.sqrt(x)*np.sqrt(1.0-x)
	# transform = np.vectorize(lambda x: np.power(x,0.75)*np.power((1.0-x),0.75))
	transform = np.vectorize(lambda x: x)
	# 
	xxx = range(1,22)
	#
	ax.plot(xxx,transform(fmin_arr),color='blue',linestyle='-',marker='_',mec='black',mew=2.0,ms=10)
	ax.plot(xxx,transform(fmax_arr),color='red',linestyle='-',marker='_',mec='black',mew=2.0,ms=10)
	#
	ax.fill_between(xxx, transform(fmax_arr), transform(fmin_arr), facecolor=fillcolor, alpha=0.5)
	#
	ax.set_xlim(0,22)
	# ax.set_ylim(0,1.0)
	# plt.show()

# # uncomment for massive plotting ...
# # plot_haplo_range(ax,fmin_arr_init,fmax_arr_init,fillcolor='gray')
# # scan haplotype number XXX :
# scan_hap_idx = 13
# #
# #
# print "haplotype N %d looking like "%(scan_hap_idx+1)
# print phaps_list[scan_hap_idx]
# # get its global index ...
# global_haplo_index = phaps_indices[scan_hap_idx]
# # get interval initial of this guy ...
# num = 12
# _,fmin,fmax,_ = parsed_out_initial[global_haplo_index]
# #
# print "freq_range: ",fmin,fmax
# #
# scan_hap_freq = [ round(_,3) for _  in np.linspace(fmin,fmax,num) ]
# #
# #
# # scan_hap_freq = [0.272,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.381]
# # # get its global index ...
# # global_haplo_index = phaps_indices[scan_hap_idx]
# # 
# # 
# for hap_fixed_freq in scan_hap_freq:
# 	# add 1 more link to fixate that haplotype ...
# 	linkages = link
# 	links_freqs_str_ADD = "LINK%d %d %.3f\n" % ( linkages, loci, hap_fixed_freq )
# 	linked_la_str_ADD = ""
# 	for lid in range(loci):
# 		linked_la_str_ADD += "LINK%d L%d A%d\n"%( linkages, lid+1, phaps_list[scan_hap_idx][lid])
# 	link_info_out = "LINKAGE %d\n"%linkages + (links_freqs_str+links_freqs_str_ADD) + (linked_la_str+linked_la_str_ADD)
# 	# pop_info_out
# 	# output file ...
# 	filename = "input_test_%d_%.3f.dat"%(scan_hap_idx+1,hap_fixed_freq) 
# 	with open(filename,'w') as fp:
# 		fp.write( pop_info_out + link_info_out )
# 	#
# 	# after the file is wsritten, we can go ahead and launch haplotest program!!!
# 	output_lines = sub.check_output("./bin/testhaplo %s"%filename,shell=True).strip().split('\n')
# 	parsed_out = parse_haplolimit_output(output_lines)
# 	intervals = get_intervals(parsed_out)
# 	ratios = get_interval_ratios(intervals,intervals_initial)
# 	# calculate ratios excluding that haplotype (it if fixed and always result in a zero result!) ...
# 	ratios_product = prod(ratios[:global_haplo_index]+ratios[(global_haplo_index+1):])
# 	# print hap_fixed_freq,ratios_product
# 	#
# 	#
# 	fmin_arr_constr = []
# 	fmax_arr_constr = []
# 	for phap_idx in range(21):
# 		# get its global index ...
# 		global_haplo_index = phaps_indices[phap_idx]
# 		# get interval initial of this guy ...
# 		_,fmin,fmax,_ = parsed_out[global_haplo_index]
# 		# ...
# 		fmin_arr_constr.append(fmin)
# 		fmax_arr_constr.append(fmax)
# 	#
# 	#
# 	# uncomment if plotting is needed ...
# 	ax.cla()
# 	plot_haplo_range(ax,fmin_arr_init,fmax_arr_init,fillcolor='gray')
# 	plot_haplo_range(ax,fmin_arr_constr,fmax_arr_constr,fillcolor='red')
# 	plt.savefig("scan_h%d_f%.3f.png"%((scan_hap_idx+1),hap_fixed_freq),dpi=150)


# plt.show()
# plt.savefig("scan_h%d_f%.3f.pdf"%((scan_hap_idx+1),scan_hap_freq[0]))

# # print "LINKAGE %d"%(link-1)
# # print links_freqs_str + linked_la_str
# # print linked_la_str

# # link_info_out = "LINKAGE %d\n"%(link-1) + links_freqs_str + linked_la_str







