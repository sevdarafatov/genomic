import os
import vcf

class trio:
    def __init__(self, vcf_proband, vcf_mother, vcf_father):
        self.vcfproband = vcf_proband
        self.vcfmother = vcf_mother
        self.vcffather = vcf_father
        self.outputFolder = "/home/gnks/Desktop/sevda/trio_class/output"
        self.commandMaker()
        self.counter()
        
    
    def commandMaker(self):
        mergeCommand = f'bcftools merge {self.vcfproband} {self.vcfmother} {self.vcffather} > output_new.txt'
        os.system(mergeCommand)
        
    def counter(self):
        # Open merged VCF file for reading and new VCF file for writing
        vcf_reader = vcf.Reader(open('output_new.txt', 'r'))


        # Add new 'DeNovo' INFO field to VCF headerla
        vcf_reader.infos['DeNovo'] = vcf.parser._Info('DeNovo', 0, 'Flag', 'Variant is de novo in proband', source=None, version=None)
        vcf_reader.infos['mutant_het_couple'] = vcf.parser._Info('mutant_het_couple', 0, 'Flag', 'couple is heterozygote for the mutation', source=None, version=None)
        vcf_reader.infos['segregated_homo'] = vcf.parser._Info('segregated_homo', 0, 'Flag', 'Proband is Homozygote and parents are heterozygote', source=None, version=None)
        vcf_reader.infos['PR_GT'] = vcf.parser._Info('PR_GT', 0, 'Flag', 'GT of proband', source=None, version=None)
        vcf_reader.infos['M_GT'] = vcf.parser._Info('MT_GT', 0, 'Flag', 'GT of mother', source=None, version=None)
        vcf_reader.infos['F_GT'] = vcf.parser._Info('FT_GT', 0, 'Flag', 'FT of father', source=None, version=None)

        
        vcf_writer = vcf.Writer(open('output_with_info_FINAL.txt', 'w'), vcf_reader)
        
        # Iterate through each record in the VCF file
        for record in vcf_reader:
            # Remove the AN field from the VCF file
            del record.INFO['AN']
            # Get the genotypes for the proband, mother, and father
            gt_value_proband = record.samples[0]['GT']
            gt_value_mother = record.samples[1]['GT']
            gt_value_father = record.samples[2]['GT']

            # Check if proband is homozygous and parents are homozygous wildtype
            
            if (gt_value_proband in ('0/1', '0|1') and gt_value_mother in ('0/0', '0|0')) and gt_value_father in ('0/0', '0|0'):
                # If above condition is true, mark variant as de novo
                record.INFO['DeNovo'] = 'T'
            
            # Check if mutant hetero couple
            if  gt_value_mother in ('0/1', '0|1') and gt_value_father in ('0/1', '0|1'):
                
                record.INFO['mutant_het_couple'] = 'T'
                
                if  (gt_value_proband in ('1/1', '1|1')):
                    #Check if segregated
                    record.INFO['segregated_homo'] = 'T'
            AN = 6
            AC = 0
            # Defining GT values for proband       
            if  gt_value_proband in ('0/1', '0|1'):
                record.INFO['PR_GT'] = 'Het'
                AC += 1
            if  gt_value_proband in ('1/1', '1|1') :
                record.INFO['PR_GT'] = 'Homo'
                AC += 2
            if  gt_value_proband in ('0/0', '0|0') :
                record.INFO['PR_GT'] = 'Hom-ref'
                
            # Defining GT values for mother  
            if  gt_value_mother in ('0/1', '0|1'):
                record.INFO['M_GT'] = 'Het'
                AC += 1
            if  gt_value_mother in ('1/1', '1|1') :
                record.INFO['M_GT'] = 'Homo'
                AC += 2
            if  gt_value_mother in ('0/0', '0|0') :
                record.INFO['M_GT'] = 'Hom-ref'
                
            # Defining GT values for father
            if  gt_value_father in ('0/1', '0|1'):
                record.INFO['F_GT'] = 'Het'
                AC += 1 
            if  gt_value_father in ('1/1', '1|1') :
                record.INFO['F_GT'] = 'Homo'
                AC += 2
            if  gt_value_father in ('0/0', '0|0') :
                record.INFO['F_GT'] = 'Hom-ref'
            # Write AC and AN to info
            record.add_info("AC", AC)
            record.add_info("AN", AN)
                
            # Write the modified record to the new VCF file
            vcf_writer.write_record(record)

        vcf_writer.close()

                                
trio("/home/gnks/Desktop/sevda/trio_class/input/vcfs/4-G-EnDem_hg38.hard-filtered.vcf.gz", "/home/gnks/Desktop/sevda/trio_class/input/vcfs/147-G-CeDem_hg38.hard-filtered.vcf.gz", "/home/gnks/Desktop/sevda/trio_class/input/vcfs/148-G-IsDem_hg38.hard-filtered.vcf.gz" )
