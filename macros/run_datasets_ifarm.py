import os

def write_config(config_dict,filepath):
    outfile = open(filepath,'w')
    lines = []
    for item in config_dict:
	lines.append("%s %s\n" % (item,config_dict[item]))
    outfile.writelines(lines)
    outfile.close()


base_dir = '/w/halld-scifs17exp/home/yunjiey/run_FastDIRC_track_selection'
#treefiles_dir   = '/lustre/expphy/volatile/halld/home/jrsteven/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass05/merged'
treefiles_dir   = '/lustre/expphy/volatile/halld/home/jrsteven/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass09/merged'
#treefiles_dir   = '/lustre/expphy/volatile/halld/home/jrsteven/RunPeriod-2019-01/recon/analysis_REST/ver01_pass00/merged'



#outfiles_dir  = base_dir + '/result_dir/output_rootfiles'
outfiles_dir  = base_dir + '/result_dir_ver08_pass09/output_rootfiles'
#outfiles_dir  = base_dir + '/result_dir_ver01_pass00/output_rootfiles'

geometry_file = base_dir + '/bayes_calib/FastDIRC_dev/geometry_files/FastDIRC_HDDS_Nominal.csv'

config_dict = {}

dataset_dict = {}
dataset_dict["LI1"] = [60700,60738]
dataset_dict["LI2"] = [60739,60763]
dataset_dict["LI3"] = [60770,60788]
dataset_dict["MI1"] = [60789,60799]
dataset_dict["MI2_1"] = [60830,60837]
dataset_dict["MI2_2"] = [60842,60848]
dataset_dict["HI1"] = [60811,60829]
dataset_dict["HI2"] = [60838,60841]

#datasets_to_run = ["LI1"]
#datasets_to_run = ["LI2","LI3"]
#datasets_to_run = ["MI1","MI2_1","MI2_2"]
#datasets_to_run = ["HI1","HI2"]
datasets_to_run = ["LI1","LI2","LI3","MI1","MI2_1","MI2_2","HI1","HI2"]

cmd = 'cp ' + base_dir +'/bayes_calib/FastDIRC_dev/fastdirc_exe .'
print cmd
os.system(cmd)

for dataset in datasets_to_run:
    dataset_dir = '%s/%s'%(outfiles_dir,dataset)
    cmd = 'mkdir %s'%dataset_dir
    if not os.path.exists(dataset_dir):
    	print cmd
        os.system(cmd)
    start_RunNumber = dataset_dict[dataset][0]
    end_RunNumber   = dataset_dict[dataset][1]
    tree_counter = 0
    runs_analyzed = []
    for run_i in range(start_RunNumber,end_RunNumber+1):
	treename = treefiles_dir + "/hd_root_0%d.root"%run_i
	if not os.path.exists(treename):
	    continue

	tree_counter += 1
	outfiles_dir_run = '%s/0%d'%(dataset_dir,run_i)

	if not os.path.exists(outfiles_dir_run):
	    cmd = 'mkdir %s' % (outfiles_dir_run)
	    print cmd
	    os.system(cmd)

	config_dict["DIRCTREE_INFILE"] = "'%s'"%(treename)
	config_dict["GEOMETRY_INFILE"] = "'%s'"%(geometry_file)
	config_dict["OUTFILE"] = "'%s/output_hists_0%d.root'"%(outfiles_dir_run,run_i)
	config_dict["n_phi_phots"] = "12500"
	config_dict["Ntracks"] = "-1"
	config_dict["Ntracks_pion"] = "1000"
	config_dict["Ntracks_kaon"] = "-1"

	config_file = "%s/config_0%d.in"%(outfiles_dir_run,run_i)
	write_config(config_dict,config_file)

	cmd = './fastdirc_exe %s'%(config_file)
	print 'begin running...'
	os.system(cmd)
	print 'run finished.'
	runs_analyzed.append(run_i)

    if len(runs_analyzed) == 0:
	continue

    cmd_hadd = 'hadd -f %s/output_hists_0%d_0%d.root '%(dataset_dir,runs_analyzed[0],runs_analyzed[-1])
    for run_j in runs_analyzed:
	cmd_hadd += ' %s/0%d/output_hists_0%d.root'% (dataset_dir,run_j,run_j)
    print cmd_hadd
    os.system(cmd_hadd)

    print "Finished running dataset %s" % (dataset)
