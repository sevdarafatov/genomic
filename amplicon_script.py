import argparse
import pathlib
import os

class WES: # Defining class WES
    def __init__(self, ref, kit, folder, outputdir):
        self.ref = ref
        self.kit = kit
        self.folder = folder
        self.outputdir = outputdir
        self.kitpath = self.check_kit()
        self.samples = self.run_wes()
        self.runcommand = self.runCommand()
        self.arrangeoutputs = self.arrange_outputs()
        self.checknofnormallyfinished = self.check_nof_normally_finished()
        
        
    def check_kit(self): # Defining a function to check kits
        kits = {
            "cftr_hg38" : "/staging/01.Clean_fastq_reads/01.General_fastq_files/01.WES_files/tek_gen/CFTR/CFTR_hg38.bed",
            "mefv_hg38": "/staging/01.Clean_fastq_reads/01.General_fastq_files/01.WES_files/tek_gen/MEFV/MEFV_hg38.bed",
            "chrM": "/staging/01.Clean_fastq_reads/01.General_fastq_files/01.WES_files/tek_gen/mtDNA/Twist_MitoPanel_chrM_all_hg38_target.bed",
            "paragon_mtDNA": "/staging/01.Clean_fastq_reads/01.General_fastq_files/01.WES_files/tek_gen/mtDNA/Paragon_mito_kit.bed"
        }
        
        try:
            return kits[self.kit]
        except Exception as e:
            print(e)
            exit()

    def run_wes(self):
        pathlib.Path(self.outputdir).mkdir(parents=True, exist_ok=True)
        listoffasts = []
        for file in os.listdir(self.folder):
            if file.endswith("_1.fq.gz"):
                full_path_r1 = os.path.join(self.folder, file)
                r2 = file.replace("_1.fq.gz","_2.fq.gz")
                fullpath_r2 = os.path.join(self.folder, r2)
                rgsm = file.split(".")[0]
                if os.path.exists(fullpath_r2):
                    listoffasts.append([full_path_r1, fullpath_r2, rgsm])
            
            elif file.endswith("R1_001.fastq.gz"):
                r2 = file.replace("R1_001.fastq.gz","R2_001.fastq.gz")
                full_path_forward = os.path.join(self.folder, file)
                fullpath_reverse = os.path.join(self.folder, r2)
                rgsm = file.split(".")[0]
                if os.path.exists(fullpath_reverse):
                    listoffasts.append([full_path_forward, fullpath_reverse, rgsm])
        return listoffasts
            
    def runCommand(self):
        # loop samples
        for sample in self.samples:
            cmd = f'dragen -f -r {self.ref} --enable-variant-caller true --enable-map-align true --enable-map-align-output true --enable-duplicate-marking=false  --enable-cnv true --cnv-target-bed {self.kitpath} --remove-duplicates=false --output-directory {self.outputdir} --output-file-prefix {sample[2]} --vc-target-bed {self.kitpath} --vc-target-bed-padding 75 --qc-coverage-region-1 {self.kitpath} --RGID {sample[2]} --RGSM {sample[2]} -1 {sample[0]} -2 {sample[1]}'
            os.system(cmd)
    
    def arrange_outputs(self):
        os.system(f'mv -fn {self.outputdir}/*usage* /staging/04.Output_Vcf_Bam/Usage/')
        pathlib.Path(self.outputdir + '/others').mkdir(parents=True, exist_ok=True)
        os.system(f'mv -fn {self.outputdir}/* {self.outputdir}/others')
        pathlib.Path(self.outputdir + '/vcfs').mkdir(parents=True, exist_ok=True)
        os.system(f'mv -fn {self.outputdir}/others/*vcf* {self.outputdir}/vcfs')
        pathlib.Path(self.outputdir + '/bams').mkdir(parents=True, exist_ok=True)
        os.system(f'mv -fn {self.outputdir}/others/*bam* {self.outputdir}/bams')

    def check_nof_normally_finished(self):
        file = open(self.folder + 'nohup.out', 'r')
        nohup_out = file.read()
        n_of_normally_finished = nohup_out.count("DRAGEN finished normally")
        print(f'DRAGEN finished {n_of_normally_finished} files normally')

if __name__== "__main__":
        parser = argparse.ArgumentParser(description ='Dragen WES')
        parser.add_argument('-r', type = str, required = True, help = 'Reference Genome')
        parser.add_argument('-k', type = str, required = True, help = 'kit type')
        parser.add_argument('-f', type = str, required = True, help = 'Fastq file folder')
        parser.add_argument('-o', type = str, required = True, help = 'Output directory')
        args = parser.parse_args()
        
        WES(args.r, args.k, args.f, args.o)
