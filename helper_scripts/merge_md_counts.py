import pandas as pd
import os
import sys, getopt

argv=sys.argv[1:]
opts,args = getopt.getopt(argv,"hi:m:p:n:r:",["in=", "md-ref=", "project=", "nfilt=", "rfilt="])	

nfilt = 100
rfilt = 100

for opt,arg in opts:
    if opt == '-h':
        print("merge_md.py -i <input path> -p <project name> -n <filter contig counts by this number> -r <filter read counts by this number>")
        sys.exit()
    elif opt in ("-i", "--in"):
    	summary_count_path = arg    	
    elif opt in ("-m", "--md-ref"):
        main_ref_file = arg
    elif opt in ("-p", "--project"):
        project = arg
    elif opt in ("-n", "--nfilt"):
        nfilt = arg
    elif opt in ("-r", "--rfilt"):
        rfilt = arg


meta_contig_files_tmp = os.popen("ls "+summary_count_path+"/*_metaspades_contigs_blastx_daa_summary_count.tsv").read().split('\n')
read_files_tmp = os.popen("ls "+summary_count_path+"/*_reads_blastx_daa_summary_count.tsv").read().split('\n')
meta_contig_files = [x for x in meta_contig_files_tmp if x]
read_files = [x for x in read_files_tmp if x]
        
if "-m" not in list(sum(opts,())) and "--md-ref" not in list(sum(opts,())):
	main_frame_save = pd.read_csv(meta_contig_files[0], sep='\t', header=None, names=["TaxID","Taxa","tmp"],usecols=[0,1])
	num_rows = len(main_frame_save)
	print(num_rows)
	new_col = ["NT_count"]*num_rows
	main_frame_save["variable"] = new_col
else:
	main_frame_save = pd.read_csv(main_ref_file)

main_frame = main_frame_save

for file in meta_contig_files:
	sample = os.path.basename(file).replace("_metaspades_contigs_blastx_daa_summary_count.tsv","")
	filter_tempa=pd.DataFrame()
	filter_tempb=pd.DataFrame()
	try:
		temp_frame1 = pd.read_csv(summary_count_path+sample+"_metaspades_contigs_blastx_daa_summary_count.tsv", sep='\t', header=None, names=["TaxID","Taxa", sample])
		temp_framea = temp_frame1[["TaxID","Taxa",sample]]
		filter_tempa = temp_framea.loc[temp_framea[sample] > nfilt]
		num_rows = len(filter_tempa)
		new_col = ["NT_count"]*num_rows
		filter_tempa["variable"] = new_col
	except FileNotFoundError:
		print("File Not Found: " +summary_count_path+sample+"_metaspades_contigs_blastx_daa_summary_count.tsv")
		continue
	try:
		temp_frame2 = pd.read_csv(summary_count_path+sample+"_short_reads_blastx_daa_summary_count.tsv", sep='\t', header=None, names=["TaxID","Taxa", sample])
		temp_frameb = temp_frame2[["TaxID","Taxa",sample]]
		filter_tempb = temp_frameb.loc[temp_frameb[sample] > rfilt]
		num_rows = len(filter_tempb)
		#print(num_rows)
		new_col = ["Read_count"]*num_rows
		filter_tempb["variable"] = new_col
	except FileNotFoundError:
		print("File Not Found: "+summary_count_path+sample+"_short_reads_blastx_daa_summary_count.tsv")

	if not filter_tempa.empty:
		temp_frame = filter_tempa
		if not filter_tempb.empty:
			temp_frame = pd.concat([filter_tempa, filter_tempb])
	elif not filter_tempb.empty:
		temp_frame = filter_tempb
	curr_main_frame = pd.merge(main_frame, temp_frame, how="outer", on=["TaxID", "Taxa", "variable"])
	main_frame = curr_main_frame


if "-p" not in list(sum(opts,())) and "--project" not in list(sum(opts,())):
	main_frame.to_csv(summary_count_path+"MD_counts_merged.csv")
else:
	main_frame.to_csv(summary_count_path+"/"+project+"_MD_counts_merged.csv")
