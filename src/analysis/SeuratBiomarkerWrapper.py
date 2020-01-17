import sys, getopt, subprocess

#use Rscript to call the Seurat clustering script
def seurat_clusters_script():
    #call the Seurat command
    seuratcommandcall = "/usr/local/bin/Rscript SINGLE_CLUSTERING.R"
    subprocess.call(seuratcommandcall, shell=True)


def main(argv):
    
    #parse command line arguments
    inputgenelist=''
    
    try:
        opts, args = getopt.getopt(argv, "hx:g:",["ifile1=","ifile2="])
    except getopt.GetoptError:
        print('SeuratBiomarkerWrapper.py -g <genes>')
        sys.exit(3)
        
    for opt, arg in opts:
        if opt == '-h':
            print('SeuratBiomarkerWrapper.py -g <genes>')
            sys.exit()
        elif opt in ("-g", "--genes"):
            inputgenelist = arg
            
    #error conditions: missing input and confliciting command line arguments
    if inputgenelist=='':
        print("Error: Missing input")
        sys.exit(3)
        
    #convert input sequences into bytes
    sgenes = inputgenelist.split(",")
    print(sgenes)
    
    seurat_clusters_script()
    
if __name__== "__main__":
    main(sys.argv[1:])