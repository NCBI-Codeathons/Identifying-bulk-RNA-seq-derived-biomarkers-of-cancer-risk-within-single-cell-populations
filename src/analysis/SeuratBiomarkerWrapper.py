import sys, getopt, subprocess

from gsea_analysis import get_gsea_result
from color_genes_on_phate import color_genes_on_phate

#use Rscript to call the Seurat clustering script
def seurat_clusters_script(clusteringres):
    #call the Seurat command
    
    seuratcommandcall = "/usr/local/bin/Rscript SINGLE_CLUSTERING.R"
    subprocess.call(seuratcommandcall, shell=True)


def main(argv):
    
    #parse command line arguments
    inputgenelist=''
    clusteringres=''
    
    try:
        opts, args = getopt.getopt(argv, "hg:r:",["ifile1=","ifile2="])
    except getopt.GetoptError:
        print('SeuratBiomarkerWrapper.py -g <genes> -r <resolution>')
        sys.exit(3)
        
    for opt, arg in opts:
        if opt == '-h':
            print('SeuratBiomarkerWrapper.py -g <genes> -r <resolution')
            sys.exit()
        elif opt in ("-g", "--genes"):
            inputgenelist = arg
        elif opt in ("-r", "--resolution"):
            clusteringres = arg
            
    #error conditions: missing input and confliciting command line arguments
    if inputgenelist=='':
        print("Error: Missing biomarker gene list input")
        sys.exit(3)
        
    #split the input of gene names
    sgenes = inputgenelist.split(",")
    print(sgenes)
    
    seurat_clusters_script(clusteringres)
    
if __name__== "__main__":
    main(sys.argv[1:])