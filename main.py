input_file_handle = open(input())
unparsed_input = input_file_handle.readlines()
input_file_handle.close()

plasmid = unparsed_input[0].rstrip()
plasmid_restriction_site = unparsed_input[1].rstrip()
GFP = unparsed_input[2].rstrip()
restriction_sites = unparsed_input[3].split(" ")
GFP_site_1 = restriction_sites[0]
GFP_site_2 = restriction_sites[1].rstrip()


def getComplementaryStrand(strand):
    complementarity = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complementary = ""
    for base in strand:
        complementary = complementary + complementarity[base]
    return complementary


def findCutPoint(strand, site, reverse):
    if reverse:
        return strand.rfind(site)
    return strand.find(site)


def cutStrand(strand, search_start, cutpoint, i):
    cut_strand = ""
    for base in strand:
        if search_start < cutpoint + i:
            cut_strand = cut_strand + strand[search_start]
            search_start = search_start + 1
    return cut_strand


def cutGeneStrand(strand, cutpoint, site_index):
    cut_strand = ""
    i = cutpoint + site_index
    while i < len(strand):
        cut_strand = cut_strand + strand[i]
        i = i + 1
    return cut_strand


def completeLigationPrint(plasmid_1, plasmid_2, plasmid_3, plasmid_4, orig_GFP, comp_GFP):
    print(plasmid_1 + orig_GFP + plasmid_2)
    print(plasmid_3 + comp_GFP + plasmid_4)


# first two fragments of plasmid
plasmid_1 = cutStrand(plasmid, 0, findCutPoint(plasmid, plasmid_restriction_site, False), 1)
plasmid_2 = cutStrand(plasmid, len(plasmid_1), len(plasmid), 0)

# plasmid complements
plasmid_complement = getComplementaryStrand(plasmid)
plasmid_restriction_site_complement = getComplementaryStrand(plasmid_restriction_site)

# last two fragments of plasmid
plasmid_3 = cutStrand(plasmid_complement,0,findCutPoint(plasmid_complement,plasmid_restriction_site_complement,False),5)
plasmid_4 = cutStrand(plasmid_complement, len(plasmid_3), len(plasmid), 0)

# complements of GFP
GFP_complement = getComplementaryStrand(GFP)
GFP_site_1_complement = getComplementaryStrand(GFP_site_1)
GFP_site_2_complement = getComplementaryStrand(GFP_site_2)

# cutting GFP at the second restriction site
GFP_1 = cutStrand(GFP, 0, findCutPoint(GFP, GFP_site_2, True), 1)

# first sticky complete strand of cut GFP
GFP_2 = cutGeneStrand(GFP_1, findCutPoint(GFP_1, GFP_site_1, False), 1)

# cutting GFP complement at the second site
GFP_3 = cutStrand(GFP_complement, 0, findCutPoint(GFP_complement, GFP_site_2_complement, True), 5)

# cutting GFP complement at the second site, second sticky complete strand
GFP_4 = cutGeneStrand(GFP_3, findCutPoint(GFP_3, GFP_site_1_complement, False), 5)

# printing the ligated plasmid
completeLigationPrint(plasmid_1, plasmid_2, plasmid_3, plasmid_4, GFP_2, GFP_4)
