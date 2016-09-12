import numpy as np
import sys
import re
import itertools
import subprocess as sub

import matplotlib.pyplot as plt


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



def parse_pop(filename):
	with open(filename,'r') as fp:
		# first line ...
		line = fp.readline()
		loci, = re.match("LOCI (\d+)",line.strip()).groups()
		loci = int(loci)
		# loci lines ...
		alleles = [0,]*loci
		for locus in range(loci):
			idx,num_al = re.match("L(\d+) (\d+)",fp.readline().strip()).groups()
			alleles[int(idx)-1] = int(num_al)
		# total number of alleles ...
		tot_al = sum(alleles)
		# allele frequencies ...
		freq = []
		for locus in range(loci):
			freq.append([])
			for allele in range(alleles[locus]):
				freq[locus].append(0.0)
		########################################
		for allele in range(tot_al):
			lidx,aidx,fin = re.match("L(\d+) A(\d+) (\d+\.\d+|AUTO)",fp.readline().strip()).groups()
			lidx = int(lidx)-1
			aidx = int(aidx)-1
			freq[lidx][aidx] = float(fin) if fin != "AUTO" else round(1.0 - sum(freq[lidx]),3)
	return freq


def check_pop(freq_arr):
	for idx,locus_freqs in enumerate(freq_arr):
		freq_norm = round(sum(locus_freqs),3)
		if freq_norm != 1.0:
			print "Troubles parsing population info: normalization failed @ locus %d"%(idx+1)
			return False
	return True


def generate_compliant_popinfo(freq_arr):
	out_loci_string = ""
	out_freq_string = ""
	nloci = len(freq_arr)
	out_loci_string += "LOCI %d\n"%nloci
	for lid in range(nloci):
		alleles_at_locus = len(freq_arr[lid])
		out_loci_string += "L%d %d\n"%(lid+1,alleles_at_locus)
		for aid in range(alleles_at_locus-1):
			out_freq_string += "L%d A%d %.3f\n"%( lid+1, aid+1, round(freq_arr[lid][aid],3) )
		out_freq_string += "L%d A%d AUTO\n"%( lid+1, alleles_at_locus )
	return out_loci_string + out_freq_string



def parse_haplolimit_outfile(fname):
	#
	index_pattern = re.compile("H(\d+)")
	haplo_pattern = re.compile("L(\d+)A(\d+)")
	# skip 7 lines and then go on
	with open(fname,'r') as fp:
		skip_line_number = 7
		for skip_line in range(skip_line_number):
			fp.readline()
		# start parsing the thing ...
		haplo_out = []
		for line in fp.readlines():
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



# get intervals for parsed haplo out ...
# (5, 0.0, 0.014, [1, 1, 1, 5]),
def get_intervals(haplo_out):
	to_ret = []
	for entry in haplo_out:
		_,fmin,fmax,_ = entry
		to_ret.append(fmax-fmin)
	return to_ret


def get_interval_ratios(int_numer,int_denom):
	to_ret = []
	# index = 0
	for intn,intd in zip(int_numer,int_denom):
		ratio = intn/intd if intd else 1.0
		to_ret.append(ratio)
		# index += 1
	return to_ret


def prod(array):
	prod_accum = 1.0
	for item in array:
		prod_accum *= item
	return prod_accum



phaps = parse_haplos("./input/possibleHaps.txt")

# freqs = parse_pop("good_test.dat")
freqs = parse_pop("./input/lociInfo.txt")


check_pop(freqs)

pop_info_out = generate_compliant_popinfo(freqs)

phaps = np.asarray(phaps)
# print phaps


haps,loci = phaps.shape
#reconstruct alleles:
alleles = []
for locus in range(loci):
	alleles.append( phaps[:,locus].max() )

children_at_locus = []
for locus in range(loci-1):
	children_at_locus.append([])
	for allele in range(alleles[locus]):
		children_at_locus[locus].append([])
#################################
# groups_at_locus = [[]]*loci


# fill out those children of at different loci ...
for locus in range(loci-1):
	for hap in range(haps):
		children_at_locus[locus][ phaps[hap,locus]-1 ].append( phaps[hap,locus+1] )
	# print children_at_locus[locus][ phaps[hap,locus]-1 ]
	for hap in range(haps):
		children_at_locus[locus][ phaps[hap,locus]-1 ] = tuple(sorted(set(children_at_locus[locus][ phaps[hap,locus]-1 ])))
	children_at_locus[locus] = list(enumerate(children_at_locus[locus]))

# first key in a tuple ...
key = lambda tup: tup[1]

link = 1
links_freqs_str = ""
linked_la_str = ""
for locus in range(loci-1):
	for shared_children,groups in itertools.groupby(sorted(children_at_locus[locus],key=key),key=key):
		# print [parent[0]+1 for parent in groups],shared_children
		# each of this things is a link!!!
		parents = [parent[0]+1 for parent in groups] # @ locus
		children = shared_children # @ locus+1
		# let's calculate some frequencies ...
		freq1 = round(sum([freqs[locus][ppp-1] for ppp in parents]),3)
		freq2 = round(sum([freqs[locus+1][ccc-1] for ccc in children]),3)
		# if (freq1 == freq2):
		actual_freq = min(freq1,freq2)
		links_freqs_str += "LINK%d %d %.3f\n" % (link,len(parents)+len(children),actual_freq)
		# else:
		# 	print "ACHTUNG ACHTUNG ACHTUNG!!!!!!!!! %.2f != %.2f"%(freq1,freq2) 
		for ppp in parents:
			linked_la_str += "LINK%d L%d A%d\n"%(link,locus+1,ppp)
		for ccc in children:
			linked_la_str += "LINK%d L%d A%d\n"%(link,locus+2,ccc)
		link += 1



# let's get the intervals first - without ever restricting any haplotypes ...
linkages = link-1
link_info_out = "LINKAGE %d\n"%linkages + links_freqs_str + linked_la_str
filename = "input_test_ASIS.dat" 
with open(filename,'w') as fp:
	fp.write( pop_info_out + link_info_out )
#
# after the file is written, we can go ahead and launch haplotest program!!!
output_lines = sub.check_output("./bin/testhaplo %s"%filename,shell=True).strip().split('\n')
parsed_out_initial = parse_haplolimit_output(output_lines)
intervals_initial = get_intervals(parsed_out_initial)

# let's get bare haplotypes ...
haplotypes_only_initial = [_[3] for _ in parsed_out_initial]

# transform phaps back to python lists (silly but ...)
phaps_list = [list(_) for _ in phaps]


def index_possible_haps(pos_haps,haps_initial):
	indices = []
	for the_hap in pos_haps:
		the_hap_index = haplotypes_only_initial.index(the_hap)
		indices.append(the_hap_index)
	return indices

# now get the indices of those haplotypes ...
phaps_indices = enumerate(index_possible_haps(phaps_list,haplotypes_only_initial))
# once ready, sort em by the fmax!
indices_sorted_by_fmax = list(sorted(
							[ (_pos,_glob, parsed_out_initial[_glob][2]) for _pos,_glob in phaps_indices],
							key=lambda x: x[2],
							reverse=True ))


# extracting re-sorted data here ...
phaps_list = [ phaps_list[idx] for idx,_,_ in indices_sorted_by_fmax ]
phaps_indices = [ idx for _,idx,_ in indices_sorted_by_fmax ]


# let's draw stuff ...
fig,ax = plt.subplots(figsize=(10,4))


fmin_arr_init = []
fmax_arr_init = []
for phap_idx in range(21):
	# get its global index ...
	global_haplo_index = phaps_indices[phap_idx]
	# get interval initial of this guy ...
	_,fmin,fmax,_ = parsed_out_initial[global_haplo_index]
	# ...
	fmin_arr_init.append(fmin)
	fmax_arr_init.append(fmax)
	#


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

# scan haplotype number XXX :
for scan_hap_idx in range(21):
	print "haplotype N %d looking like "%(scan_hap_idx+1)
	print phaps_list[scan_hap_idx]
	print
	# scan_hap_freq = [0.272,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.381]
	# get its global index ...
	global_haplo_index = phaps_indices[scan_hap_idx]
	# get interval initial of this guy ...
	num = 12
	_,fmin,fmax,_ = parsed_out_initial[global_haplo_index]
	scan_hap_freq = [ round(_,3) for _  in np.linspace(fmin,fmax,num) ]
	# 
	#
	distro = []
	#
	#
	for hap_fixed_freq in scan_hap_freq:
		# add 1 more link to fixate that haplotype ...
		linkages = link
		links_freqs_str_ADD = "LINK%d %d %.3f\n" % ( linkages, loci, hap_fixed_freq )
		linked_la_str_ADD = ""
		for lid in range(loci):
			linked_la_str_ADD += "LINK%d L%d A%d\n"%( linkages, lid+1, phaps_list[scan_hap_idx][lid])
		link_info_out = "LINKAGE %d\n"%linkages + (links_freqs_str+links_freqs_str_ADD) + (linked_la_str+linked_la_str_ADD)
		# pop_info_out
		# output file ...
		filename = "input_test_%d_%.3f.dat"%(scan_hap_idx+1,hap_fixed_freq) 
		with open(filename,'w') as fp:
			fp.write( pop_info_out + link_info_out )
		#
		# after the file is written, we can go ahead and launch haplotest program!!!
		output_lines = sub.check_output("./bin/testhaplo %s"%filename,shell=True).strip().split('\n')
		parsed_out = parse_haplolimit_output(output_lines)
		intervals = get_intervals(parsed_out)
		ratios = get_interval_ratios(intervals,intervals_initial)
		# calculate ratios excluding that haplotype (it if fixed and always result in a zero result!) ...
		# #
		# # using all ~200-300 "nonzero" haplotypes (including "prohibitet ones") 
		# ratios_product = prod(ratios[:global_haplo_index]+ratios[(global_haplo_index+1):])
		# #
		# # alternative - using just 21 possible haps excluding the one that is fixed ...
		prod_hhh = 1.0
		for hidx in phaps_indices:
			if hidx != global_haplo_index:
				prod_hhh *= ratios[hidx]
		#
		ratios_product = prod_hhh
		#
		#
		# print hap_fixed_freq,ratios_product
		#
		#
		distro.append((hap_fixed_freq,ratios_product))
	#
	# now plot violin-like plot ...
	xarr = [x for x,_ in distro]
	max_ratio = max([y for _,y in distro])
	if max_ratio == 0.0:
		yup = [ (21-scan_hap_idx)+0.3 for yy in [y for _,y in distro]]
		ydown = [ (21-scan_hap_idx)-0.3 for yy in [y for _,y in distro]]
	else:
		yup = [ (21-scan_hap_idx)+0.4*yy/max_ratio for yy in [y for _,y in distro]]
		ydown = [ (21-scan_hap_idx)-0.4*yy/max_ratio for yy in [y for _,y in distro]]
	# ax.plot(xarr,y)
	ax.fill_between(xarr, yup, ydown, facecolor="blue", alpha=0.75)

ax.set_ylim(0,22)


plt.show()

# # print "LINKAGE %d"%(link-1)
# # print links_freqs_str + linked_la_str
# # print linked_la_str

# # link_info_out = "LINKAGE %d\n"%(link-1) + links_freqs_str + linked_la_str







