import pandas as pd
import os
import json

class mosdepth:
    def __init__(self,bam,output,threads,kit):
        self.bam = bam
        self.outputFolder = output
        self.kit = kit
        self.threads = threads
        self.id = os.path.basename(self.bam).split(".")[0]
        self.out = os.path.join(self.outputFolder, self.id)
        self.commandMaker()
        self.bed_wrangler()


    def commandMaker(self):

        if not os.path.isdir(self.outputFolder):
            os.mkdir(self.outputFolder)
        print("Running the mosdepth command ...")
        mosdepthCommand = f"mosdepth {self.out} -q 9 -t {self.threads} --by {self.kit} -n -T 1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100 {self.bam}"
        os.system(mosdepthCommand)

    def bed_wrangler(self):
        print("start processing the result ...")
        file = self.out + ".thresholds.bed.gz"
        if os.path.isfile(file):
            df = pd.read_csv(file, compression='gzip', sep = '\t', dtype = str)
            header = ['#chrom','start','end','region','1X','5X','10X','15X','20X','25X','30X','35X','40X','45X','50X','55X','60X','65X','70X','75X','80X','85X','90X','95X','100X']
            
            df['range'] = df['end'].astype(int) - df['start'].astype(int) 
            df = df[['#chrom','start','end','range','region','1X','5X','10X','15X','20X','25X','30X','35X','40X','45X','50X','55X','60X','65X','70X','75X','80X','85X','90X','95X','100X']]

            #gene coverage
            #group by region
            result = {}
            df_grouped = df.groupby(["region"])
            for key in df_grouped.groups.keys():
                subset = df_grouped.get_group(key)
                #print(subset.head)
                X1_sum = subset['1X'].astype(int).sum()
                X5_sum = subset['5X'].astype(int).sum()
                X10_sum = subset['10X'].astype(int).sum()
                X15_sum = subset['15X'].astype(int).sum()
                X20_sum = subset['20X'].astype(int).sum()
                X25_sum = subset['25X'].astype(int).sum()
                X30_sum = subset['30X'].astype(int).sum()
                X35_sum = subset['35X'].astype(int).sum()
                X40_sum = subset['40X'].astype(int).sum()
                X45_sum = subset['45X'].astype(int).sum()
                X50_sum = subset['50X'].astype(int).sum()
                X55_sum = subset['55X'].astype(int).sum()
                X60_sum = subset['60X'].astype(int).sum()
                X65_sum = subset['65X'].astype(int).sum()
                X70_sum = subset['70X'].astype(int).sum()
                X75_sum = subset['75X'].astype(int).sum()
                X80_sum = subset['80X'].astype(int).sum()
                X85_sum = subset['85X'].astype(int).sum()
                X90_sum = subset['90X'].astype(int).sum()
                X95_sum = subset['95X'].astype(int).sum()
                X100_sum = subset['100X'].astype(int).sum()
                range_sum = subset['range'].astype(int).sum()

                result[key] = {
                        "1x": X1_sum/ range_sum,
                        "5x": X5_sum/ range_sum,
                        "10X": X10_sum/ range_sum,
                        "15X": X15_sum/ range_sum,
                        "20X": X20_sum/ range_sum,
                        "25X": X25_sum/ range_sum,
                        "30X": X30_sum/ range_sum,
                        "35X": X35_sum/ range_sum,
                        "40X": X40_sum/ range_sum,
                        "45X": X45_sum/ range_sum,
                        "50X": X50_sum/ range_sum,
                        "55X": X55_sum/ range_sum,
                        "60X": X60_sum/ range_sum,
                        "65X": X65_sum/ range_sum,
                        "70X": X70_sum/ range_sum,
                        "75X": X75_sum/ range_sum,
                        "80X": X80_sum/ range_sum,
                        "85X": X85_sum/ range_sum,
                        "90X": X90_sum/ range_sum,
                        "95X": X95_sum/ range_sum,
                        "100X": X100_sum/ range_sum
                    }

            with open("output_gene.json", "w") as T:
                json.dump(result, T, indent=1)       

            df2 = df
            for column in df2.columns[5:]:
                df2[column] = (df2[column].astype(int)/df2['range']) * 100
            df2 = df2.round(2)

            thresholds_bed_output = os.path.join(self.outputFolder, f"{self.id}.thresholds_PCT.json")
            df2.to_json(thresholds_bed_output, orient='records', indent=1)


mosdepth("/home/gnks/Desktop/sevda/mosdepth_trials/deneme/bams/W01-102_hg38.bam","/home/gnks/Desktop/sevda/mosdepth_trials/deneme/output", 4,"/home/gnks/Desktop/sevda/mosdepth_trials/genomic/twist2_hg38.bed")
