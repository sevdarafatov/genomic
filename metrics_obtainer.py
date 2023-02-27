with open("sample.mapping_metrics.csv", 'r') as T:
    my_list = T.readlines()
    result = {}
    for el in my_list:
        if el.startswith("MAPPING/ALIGNING PER RG"):
            el = el.split(",")[-3:]
            result[el[0]] = el[1]

    #print(result)
    newResult = {}
    
    for key in result:
        if key == "Total reads in RG":
            newResult["total_reads"] = result[key]
        if key == "'Mapped reads adjusted for filtered mapping":
            newResult["pf_reads"] = result[key]
        if key == "'Mapped reads adjusted for filtered mapping":
            newResult["pf_reads_aligned"] = result[key]
        if key == "Mapped reads R1":
            newResult["pf_reads_aligned"] = result[key]
        if key == "Paired reads (itself & mate mapped)":
            newResult["reads_aligned_in_pairs"] = result[key]
        if  key == "Not properly paired reads (discordant)":
            newResult["pf_reads_improper_pairs"] = result[key]
        if  key == "Q30 bases":
            newResult["pf_aligned_bases" and "pf_hq_aligned_bases"] = result[key]
        if  key == "Reads with MAPQ [40:inf)":
            newResult["pf_hq_aligned_reads"] = result[key]

    print(newResult)
