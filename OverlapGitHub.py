""" This script is used to select good nonoverlapping profile hits obtained from hmmscan --domtblout search. If profile hits overlap the profile hit with better
e value is selected (the overlap limit can be set at overlap_limit variable). This script accepts file format:

protein name "\t" protein length "\t" short profile annotation "\t" profile accession "\t" hit e value "\t"
hit region from - to (tells the region where on the query protein the profile hit was alligned)
which can be obtained from raw hmmscan --domtblout file by using script HitGitHub.py

This script results are given in output format:
protein "\t" protein length "\t" profile short annotation "\t" profile accession "\t" e value "\t" region of the protein |@| next good profile......"""

overlap_limit = 50 # enter value between 1 and 100. A percentage number that is tolerable when two profiles overlap.
# If the overlap is higher than the number entered the profile with better e value will be selected
ove = 100 / overlap_limit


from collections import defaultdict
import re
prot_hits = defaultdict(list) # key protein ID, value - all profile hits for the protein
prot_len = {} # key protein ID, value it's length
AccShortAnot = {} # key profile accession, value profile short annotation
fa = open("/home/paulius/BIOINFORMATICS/GitHub/HitS.txt")
lines = fa.readlines()
for line in lines:
    ID = line.split("\t")[0]
    HIT = line.split("\t")[3:]
    prot_hits.setdefault(ID, []).append(HIT)
    prot_l = line.split("\t")[1]
    prot_len[ID] = (prot_l)
    acce = line.split("\t")[3]
    shrtAn = line.split("\t")[2]
    AccShortAnot[acce] = (shrtAn)
fa.close()


f_out = open("/home/paulius/BIOINFORMATICS/GitHub/NoOverlappingHits.txt","w")
for pre in prot_hits.keys():
    hits = prot_hits[pre]
    profile_hits = defaultdict(list) # profile hits, some proteins had several hits from the same profile for its different regions,
    # so we joined them up to make one larger hit. E value for same profile hits usually were identical or similar so we picked up the longest hit e value.
    profile_joined_hits = {} # key profile accession, value hit start - hit end
    profile_hits_e_values = defaultdict(list) # key profile accession, value hit length|e value
    for hit in hits:
        accession = hit[0]
        hit_region = hit[2]
        hit_region = hit_region.replace("\n","")
        e_val = float(hit[1])
        hit_start = int(hit_region.split("-")[0])
        hit_end = int(hit_region.split("-")[1])
        hit_length = hit_end - hit_start
        profile_hits_e_values.setdefault(accession, []).append(str(hit_length) + "|" + str(e_val))
        profile_hits.setdefault(accession, []).append(hit_start)
        profile_hits.setdefault(accession, []).append(hit_end)
        
    e_value = {} # key accession, value - e value (taken from the longest profile hit!)

    for accee in profile_hits_e_values:
        values = profile_hits_e_values[accee]
        values_for_max_selection = []
        for va in values:
            lengg = va.split("|")[0]
            ev = va.split("|")[1]
            zlengg = lengg.zfill(4)
            val = zlengg + "|" + ev
            values_for_max_selection.append(val)
        good_ev = max(values_for_max_selection)
        good_eval = good_ev.split("|")[1]
        good_eval = float(good_eval)
        e_value[accee] = (good_eval)
        
            
    for acc in profile_hits:
        values = profile_hits[acc]
        min_val = min(values) # earliest profile hit start
        max_val = max(values) # latest profile hit end
        profile_joined_hits[acc] = (str(min_val) + "-" + str(max_val))

    prof_len = {} # key profile accession, value profile hit length in the protein
    prof_range = {} # key profile accession, value profile hit range in the protein
    startProtProfAcc = [] # this list is not sorted and in format min_val(4 digit number, if it lacks numbers zeros are added to the beggining)|pozicija/Accesion
    # we will sort it to get startProtProfAccSort list where all protein profile hits are ordered by the start position.
    for profile in profile_joined_hits:
        profile_region = profile_joined_hits[profile]
        prof_start = int(profile_region.split("-")[0])
        prof_end = int(profile_region.split("-")[1])
        pradziaKformatu = str(prof_start).zfill(4)
        startProtProfAcc.append(pradziaKformatu + "|" + profile_region + "/" + profile)
        rangeas = range(prof_start, prof_end +1)
        prof_leng = len(rangeas)
        ranga = set(rangeas)
        prof_len[profile] = (prof_leng)
        prof_range[profile] = (ranga)
    startProtProfAccSort = sorted(startProtProfAcc) # profiles sorted by the hit start point in the protein. This is needed
    # to be able to iterate profile hits from the start of the protein till the end of the protein from left to right and select non overlapping
    # or overlapping at most the amount allowed by the overlap_limit variable.
    
    selected_protein_profiles = []
    
    for x in startProtProfAccSort:
            
        if startProtProfAccSort.index(x) + 1 == 1 and len(startProtProfAccSort) == 1: # the protein had only one profile hit, so there is no need to check overlaps!
            selected_protein_profiles.append(x)
            
            
        if startProtProfAccSort.index(x) == 0 and len(startProtProfAccSort) != 1: # This is the first profile in protein that has multiple profile hits
            # so we will check its overlap with the 2nd profile hit!
            prfacc1 = x.split("/")[1]
            range1 = prof_range[prfacc1]
            leng1 = prof_len[prfacc1]
            prfacc2Index = startProtProfAccSort.index(x) + 1
            prfaccDu = startProtProfAccSort[prfacc2Index]
            prfacc2 = prfaccDu.split("/")[1]
            range2 = prof_range[prfacc2]
            leng2 = prof_len[prfacc2]
            Intersec = range1.intersection(range2)
            IntersecLen = len(Intersec)
            evalu1 = e_value[prfacc1]
            evalu2 = e_value[prfacc2]
            evaluu1 = float(evalu1)
            evaluu2 = float(evalu2)
            if IntersecLen > leng1 / ove or IntersecLen > leng2 / ove: # only one profile is good because overlap is higher than indicated in overlap_limit!
                if evaluu1 < evaluu2:
                    good_profile = x
                    selected_protein_profiles.append(good_profile)
                if evaluu2 < evaluu1:
                    good_profile = prfaccDu
                    selected_protein_profiles.append(good_profile)
                if evaluu1 == evaluu2:
                    if leng1 < leng2:
                        good_profile = prfaccDu
                    if leng2 < leng1:
                        good_profile = x
                    if leng2 == leng1:
                        good_profile = x
                    selected_protein_profiles.append(good_profile)
            if IntersecLen <= leng1 / ove and IntersecLen <= leng2 / ove: # both profiles are good because overlap is lower than indicated in overlap_limit!
                good_profile1 = x
                good_profile2 = prfaccDu
                selected_protein_profiles.append(good_profile1)
                selected_protein_profiles.append(good_profile2)
                
        if startProtProfAccSort.index(x) + 1 != len(startProtProfAccSort) and startProtProfAccSort.index(x) != 0:
            # This profile is neither the first nor the last profile hit in the protein so we will check it with the last selected profile hit
            # To see whether they overlap and if they do select the profile hit with better e value!
            prfacc1 = x.split("/")[1]
            range1 = prof_range[prfacc1]
            leng1 = prof_len[prfacc1]
            prfaccDu = selected_protein_profiles[-1]
            prfacc2 = prfaccDu.split("/")[1]
            range2 = prof_range[prfacc2]
            leng2 = prof_len[prfacc2]
            Intersec = range1.intersection(range2)
            IntersecLen = len(Intersec)
            evalu1 = e_value[prfacc1]
            evalu2 = e_value[prfacc2]
            evaluu1 = float(evalu1)
            evaluu2 = float(evalu2)
            if IntersecLen > leng1 / ove or IntersecLen > leng2 / ove: # only one profile is good because overlap is higher than indicated in overlap_limit!
                if evaluu1 < evaluu2:
                    good_profile = x
                    selected_protein_profiles.pop()
                    selected_protein_profiles.append(good_profile)
                if evaluu2 < evaluu1:
                    good_profile = prfaccDu
                if evaluu1 == evaluu2:
                    if leng1 < leng2:
                        good_profile = prfaccDu
                    if leng2 < leng1:
                        good_profile = x
                    if leng2 == leng1:
                        good_profile = x
                    selected_protein_profiles.pop()
                    selected_protein_profiles.append(good_profile)
            if IntersecLen <= leng1 / ove and IntersecLen <= leng2 / ove: # both profiles are good because overlap is lower than indicated in overlap_limit!
                good_profile1 = x
                good_profile2 = prfaccDu
                selected_protein_profiles.append(good_profile1) # we add only good_profile1 because good_profile2 is already in the list!

        if (startProtProfAccSort.index(x) + 1 == len(startProtProfAccSort) and startProtProfAccSort.index(x) + 2 > len(startProtProfAccSort) and
        len(startProtProfAccSort) != 1):
            # This profile hit is the last in the protein, so we will check it with the last selected profile hit to check the overlap!
            prfacc1 = x.split("/")[1]
            range1 = prof_range[prfacc1]
            leng1 = prof_len[prfacc1]
            prfaccDu = selected_protein_profiles[-1]
            prfacc2 = prfaccDu.split("/")[1]
            range2 = prof_range[prfacc2]
            leng2 = prof_len[prfacc2]
            Intersec = range1.intersection(range2)
            IntersecLen = len(Intersec)
            evalu1 = e_value[prfacc1]
            evalu2 = e_value[prfacc2]
            evaluu1 = float(evalu1)
            evaluu2 = float(evalu2)
            if IntersecLen > leng1 / ove or IntersecLen > leng2 / ove: # only one profile is good because overlap is higher than indicated in overlap_limit!
                if evaluu1 < evaluu2:
                    good_profile = x
                    selected_protein_profiles.pop()
                    selected_protein_profiles.append(good_profile)
                if evaluu2 < evaluu1:
                    good_profile = prfaccDu
                if evaluu1 == evaluu2:
                    if leng1 < leng2:
                        good_profile = prfaccDu
                    if leng2 < leng1:
                        good_profile = x
                    if leng1 == leng2:
                        good_profile = x
                    selected_protein_profiles.pop()
                    selected_protein_profiles.append(good_profile)
            if IntersecLen <= leng1 / ove and IntersecLen <= leng2 / ove: # both profiles are good because overlap is lower than indicated in overlap_limit!
                good_profile1 = x
                good_profile2 = prfaccDu
                selected_protein_profiles.append(good_profile1) # we add only good_profile1 because good_profile2 is already in the list!


    good_profile_addition = [] # good selected profiles listed in format profile short annotation "\t" profile accession "\t" e value "\t" region of the protein
    # later on we will join them to make output like this:
    # protein "\t" protein length "\t" profile short annotation "\t" profile accession "\t" e value "\t" region of the protein |@| next good profile......
    
    for selected in selected_protein_profiles:
        accen1 = selected.split("/")[1]
        annott = AccShortAnot[accen1]
        reegion = selected.split("|")[1]
        reegion = reegion.split("/")[0]
        reegion = str(reegion)
        e_vvaal = e_value[accen1]
        e_vvaal = str(e_vvaal)
        good_prof_hit = annott + "\t" + accen1 + "\t" + e_vvaal + "\t" + reegion
        good_profile_addition.append(good_prof_hit)
    added = "|@|".join(good_profile_addition)
    outPt = pre + "\t" + prot_len[pre] + "\t" + added
    f_out.write(outPt + "\n")
    selected_protein_profiles.clear()
f_out.close()





