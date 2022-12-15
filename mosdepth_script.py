import pandas as pd
import os

class mosdepth:
    def __init__(self,bam_directory,threads,kit):
        self.bamD = bam_directory
        self.outputFolder = "/home/gnks/Desktop/sevda/mosdepth_trials/genomic/output"
        self.kit = kit
        self.threads = threads
        self.bams = self.get_bams()
        self.commandMaker()
        self.bed_wrangler()

    def get_bams(self): 
        bams_list=[]
        for i in os.listdir(self.bamD):
            #gets list of bams:
            if i.endswith('bam'):
                id = os.path.basename(i)[:-4]
                path = os.path.join(self.bamD, i)
                output_dir = os.path.join(self.outputFolder, id )
                bams_list.append([id, path, output_dir])
        return bams_list

    def commandMaker(self):
        for sample in self.bams:
            if not os.path.isdir(sample[2]):
                os.mkdir(sample[2])
            out = os.path.join(sample[2], sample[0])
            mosdepthCommand = f"mosdepth {out} -q 9 -t {self.threads} --by {self.kit} -n -T 1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100 {sample[1]}"
            print(mosdepthCommand)
            os.system(mosdepthCommand)

    def bed_wrangler(self):
        for sample in self.bams:
            for i in os.listdir(sample[2]):
                if i.endswith('thresholds.bed.gz'):
                    i = os.path.join(sample[2], i)
                    print(i)

                    df = pd.read_csv(i, compression='gzip', sep = '\t', dtype = str)
                    header = ['#chrom','start','end','region','1X','5X','10X','15X','20X','25X','30X','35X','40X','45X','50X','55X','60X','65X','70X','75X','80X','85X','90X','95X','100X']
                    
                    df['range'] = df['end'].astype(int) - df['start'].astype(int)
                    df = df[['#chrom','start','end','range','region','1X','5X','10X','15X','20X','25X','30X','35X','40X','45X','50X','55X','60X','65X','70X','75X','80X','85X','90X','95X','100X']]
                    for column in df.columns[5:]:
                        df[column] = (df[column].astype(int)/df['range']) * 100
                    df = df.round(2)

                    thresholds_bed_output = os.path.join(sample[2], f"{sample[0]}.thresholds_PCT.json")
                    df.to_json(thresholds_bed_output, orient='records', indent=1)

        
mosdepth("/home/gnks/Desktop/sevda/mosdepth_trials/genomic",4,"/home/gnks/Desktop/sevda/mosdepth_trials/genomic/twist2_hg38.bed")