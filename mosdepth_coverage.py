import pandas as pd
import os
import json

class mosdepth:
    def __init__(self,bam_directory,threads,kit):
        self.bamD = bam_directory
        self.outputFolder = "/home/gnks/Desktop/sevda/27122022_run/mitochondrial/mitochondrial_bams/output"
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

                    thresholds_bed_output = os.path.join(sample[2], f"{sample[0]}.thresholds_PCT.json")
                    df2.to_json(thresholds_bed_output, orient='records', indent=1)




        
mosdepth("/home/gnks/Desktop/sevda/27122022_run/mitochondrial/mitochondrial_bams",4,"/home/gnks/Desktop/sevda/27122022_run/mitochondrial/mitochondrial_bams/bed_file/mitochondrial.bed")
