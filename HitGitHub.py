
""" This script is used to parse hmmscan --domtblout output to return all profile hits for proteins that have a better e value than e value set at
e_Vb variable. This script is written for Hit.txt example file, for your data you might want to check the float regions if your query name length differs
or if results seems odd.
The output is in format: protein name "\t" protein length "\t" short profile annotation "\t" profile accession "\t" hit e value "\t"
hit region from - to (tells the region where on the query protein the profile hit was alligned).
This output can be used by OverlapGitHub.py script to select only nonoverlapping profile hits with the best e values"""

fr = open("/home/paulius/BIOINFORMATICS/GitHub/Hit.txt")
lines = fr.readlines()
f_out = open("/home/paulius/BIOINFORMATICS/GitHub/HitS.txt","w")
sw1 = "off"
e_Vb = float(1e-3) # can be changed to set different e value limit
e_Vt = float(0)
for line in lines:
    if sw1 == "on" and line.startswith("#"): 
        sw1 = "off"
    if sw1 == "on":
        query_name = line[38:46]
        accession = line[21:31]
        accession = accession.replace(" ","")
        region_from = float(line[141:146])
        region_to = float(line[148:153])
        region_from = str(region_from)
        region_to = str(region_to)
        region_from = region_from.replace(".0","")
        region_to = region_to.replace(".0","")
        length = float(line[59:64])
        length = str(length)
        length = length.replace(".0","")
        length = length.replace(" ","")
        short_anot = line[169:-1]
        E_value = float(line[66:73])
        if e_Vt <= E_value and e_Vb >= E_value:
            f_out.write(query_name + "\t" + length + "\t" + short_anot + "\t" + accession + "\t" + str(E_value) + "\t" + str(region_from) + "-" + str(region_to)
                        + "\n")
        else:
            pass
    if sw1 == "off" and line.startswith("#-------------------"):
        sw1 = "on"
fr.close()
f_out.close()


