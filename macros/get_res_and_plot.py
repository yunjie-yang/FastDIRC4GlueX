import os,csv
def read_csv(infile):
    dict_to_return = {}
    with open(infile,'r') as f:
        reader = csv.reader(f, delimiter=',')
        for line in reader:
            dict_to_return.update({line[0].strip():float(line[1])})
    return dict_to_return


hist_type = "tree"
#hist_type = "sim"

box_label = "both"
run_type = "KinBins_all"

label = "RunLabel"
rootfile = "../output_hists.root"


momentum = 3.


plot_label = label + "_" + hist_type

plotdir  = "../plots"

resfile  = "../output_res_%s.csv"%(label)


#------------- Compute resolution --------------------#
cmd = "root -l -b -q 'graphicHistos.C(\"%s\",\"%s\",0,%.01f,0,\"%s\",\"%s\",\"%s\",\"%s\",\"%s\")'"%(rootfile,resfile,momentum,label,plotdir,hist_type,run_type,box_label)
print cmd
os.system(cmd)

#------------- Plot DLL -------------------------------#
res_dict = read_csv(resfile)
locRes = res_dict["matching_resolution"]

cmd = "root -l -b -q 'draw_DLLs_arg.C(\"%s\",\"%s\",\"%s\",%f,\"%s\",\"%s\",\"%s\")'"%(rootfile,plotdir,plot_label,locRes,hist_type,run_type,box_label)
print cmd
os.system(cmd)
