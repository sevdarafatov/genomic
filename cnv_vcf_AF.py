import vcf
import os

class cnv_AF:
    def __init__(self,cnv_directory):
        self.cnvD = cnv_directory
        self.outputFolder = "/home/sevda/Desktop/sevda/allele_frequency"
        self.cnvs = self.get_cnvs()
        self.commandMaker()
        self.count_frequency()

    def get_cnvs(self): 
        cnvs_list=[]  # creating an empty list to get cnv vcfs
        for i in os.listdir(self.cnvD): # iterate through the cnv_directory
            #gets list of cnvs: 
            if i.endswith('cnv.vcf.gz'):
                #id = os.path.basename(i)[:-11]
                path = os.path.join(self.cnvD, i)
                cnvs_list.append(path)
        return cnvs_list
    
    def commandMaker(self):
        mergeCommand = "bcftools merge " + " ".join(self.cnvs) + " > output_new.txt"
        os.system(mergeCommand)
        
    def count_frequency(self):
        vcf_reader = vcf.Reader(open('output_new.txt', 'r'))
        num_samples = len(vcf_reader.samples)
        #print(num_samples)

        # Create a dictionary to store the counts of mutations for each variant
        #variant_counts = {}

        # Create a new VCF writer
        vcf_writer = vcf.Writer(open('output_with_info_FINAL.txt', 'w'), vcf_reader)
        print(vcf_writer)
 
        # Iterate through each record in the VCF file
        for record in vcf_reader:
            # Get the variant key (chromosome and position)
            #variant_key = f"{record.CHROM}:{record.POS}"
            # Initialize the count for this variant to 0
            #variant_counts[variant_key] = 0
            count = 0
            # Iterate through each sample in the record
            for sample in record.samples:
                # Get the GT value for the sample
                gt_value = sample['GT']
                # If the GT value is not "./." or "." (i.e. not a reference genotype or missing), increment the variant count
                if gt_value != './.' and gt_value != '.':
                    count += 1
            # Add the new INFO fields to the record
            record.add_info("TM", count)
            record.add_info("SN", num_samples)
            record.add_info("AF", round(count / num_samples, 6))
            # Write the updated record to the new VCF file
            vcf_writer.write_record(record)
            
        vcf_writer.close()

cnv_AF("/home/sevda/Desktop/sevda/allele_frequency")
