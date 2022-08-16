from bs4 import BeautifulSoup
from pprint import pprint
import xml.etree.cElementTree as etree
import sys
import argparse
import pandas as pd
import numpy as np
import re

pars = argparse.ArgumentParser(description='GET INFO OF DRUGBANK')
pars.add_argument("--file", help="Input file", required=True)
pars.add_argument("-o", "--out", help="Output file")
pars.add_argument("-og", "--out_gene", help="Output file of gene")
pars.add_argument("-os", "--out_sif", help="Create a sif file from ddis")
args = pars.parse_args()

infile = open(args.file, "r")

parser = BeautifulSoup(infile.read(), "xml")

print("Fichero leido")

drugs = parser.find_all("drug", created=True)  # read entire drug entries
# drugs is now a list with each entry

print(len(drugs))
print("Fichero parseado");

if args.out is not None:
    output = open(args.out, "w")
if args.out_gene is not None:
    output_gene = open(args.out_gene, "w")

# check if sif mode:
check = False
if args.out_sif is not None:
    check = True
    print(check)


def check_if_approved(drug_id, table):
    if drug_id in list(table.index):
        group = table.loc[drug_id, 'groups']  # retrieve groups

        # if re.search("approved", group):  # alternative
        if re.search(r"\bapproved\b", group, flags=re.IGNORECASE):  # to check
            # for variations of upper/lower case of approved word
            return True  # if drug1 is approved
        else:
            return False


if check:  # if only interested in obtaining sif

    rel = {}  # to have a controls of relations

    # read drugbank info file selecting ID and group columns
    df = pd.read_table('/data/network/nas01/fivibio'
                       '-data/data/drugbank/2020/2020_07_02'
                       '/drugbank_info_2020_07.txt', index_col=0)
    # #ID col is selected as rownames

    # open output of sif-file
    sif = open(args.out_sif, "w")

    sif.write("node_1\tnode_2\tdescription\n")

    for drug in drugs:  # access xml and iterates over every drug item

        name = drug.find("name").string  # gets drug name  ".string" att is
        # used when working with beautiful soup package to extract the string
        # from the entry itself (if not it will return the entire object)
        print(name)

        # drug being studied. Needs to pass att primary as it appears in xml
        id = drug.find("drugbank-id", primary=True).string  # gets drug id
        print(id)
        node1 = str(id)

        if check_if_approved(node1, df):
            ddis = drug.find_all("drug-interactions")  # extract object that
            # contains all ddis from xml

            for i in ddis:  # makes ddis accessible to use find all
                # obtain list of all ddis
                ddi = i.find_all("drug-interaction")  # list of ddi

                for j in ddi:
                    nodes_to = j.find_all("drugbank-id")  # extract all nodes 2
                    description = j.find_all("description")  # extract info

                    for k in range(len(nodes_to)):
                        node2 = nodes_to[k]
                        if node2.string is not None:  # if node2 ID exists;
                            # print(node2.string)
                            # create control points and add to dictionary
                            # If B to A in dict:
                            if str(node2.string + "_" + node1) in rel:
                                # print(node2.string + "_" + node1)
                                continue  # stop evaluating this relation

                            else:
                                # check if node2 is approved
                                if check_if_approved(node2.string, df):
                                    # print("pass_check")
                                    des = description[k].string
                                    rel[node1 + "_" + node2.string] = 1
                                    rel[node2.string + "_" + node1] = 1
                                    sif.write(node1 + "\t" + node2.string +
                                              "\t" +
                                              des + "\n")
    # pprint(rel)
    sif.close()
    quit()  # terminate code

else:
    print("normal mode")

output.write(
    "#ID\tname\talias\tatcCodes\tsmiles\ttransporters\ttargets\tenzymes\tcarriers\tgroups\n");
output_gene.write("TYPE\tGENE\tORGANISM\n")

atcCodes = {}

for drug in drugs:

    name = drug.find("name").string
    pprint(name)

    id = drug.find("drugbank-id", primary=True)

    smiles = ""
    transporters = []
    targets = []
    carriers = []
    enzymes = []
    syns = []
    atc = []
    groups = []
    organisms = []
    row = []

    cp = drug.select("calculated-properties > property")

    for p in cp:
        if p.kind.string == "SMILES":
            smiles = p.value.string
            break

    # Transporters

    transporters_elems = drug.find_all("transporter")

    for t_elem in transporters_elems:
        aux = t_elem.find("gene-name")
        org = t_elem.find("organism")

        if aux and aux.string != None:
            transporters.append(aux.string)

        # if org and org.string != None: # To save the organism
        # organisms.append(org.string)

        if org != None and aux != None:
            row.append("transporter\t{}\t{}\n".format(aux.string, org.string))

        if org == None and aux == None:
            row.append("transporter\t{}\t{}\n".format("None", "None"))

        if org == None and aux != None:
            row.append("transporter\t{}\t{}\n".format(aux.string, "None"))
        if aux == None and org != None:
            row.append("transporter\t{}\t{}\n".format("None", org.string))

    # Targets

    targets_elems = drug.find_all("target")

    for t_elem in targets_elems:
        aux = t_elem.find("gene-name")
        org = t_elem.find("organism")

        if aux and aux.string != None:
            targets.append(aux.string)

        if org != None and aux != None:
            row.append("target\t{}\t{}\n".format(aux.string, org.string))

        if org == None and aux == None:
            row.append("target\t{}\t{}\n".format("None", "None"))

        if org == None and aux != None:
            row.append("target\t{}\t{}\n".format(aux.string, "None"))
        if aux == None and org != None:
            row.append("target\t{}\t{}\n".format("None", org.string))

    # Group

    groups_elems = drug.find_all("group")

    for group_elem in groups_elems:
        groups.append(group_elem.string)

    for t_elem in targets_elems:
        aux = t_elem.find("gene-name")
        if aux and aux.string != None:
            targets.append(aux.string)

    # Enzymes

    enzymes_elems = drug.find_all("enzyme")

    for e_elem in enzymes_elems:
        aux = e_elem.find("gene-name")
        org = e_elem.find("organism")

        # print(aux,org)

        if aux and aux.string != None:
            enzymes.append(aux.string)

        if org != None and aux != None:
            row.append("enzyme\t{}\t{}\n".format(aux.string, org.string))

        if org == None and aux == None:
            row.append("enzyme\t{}\t{}\n".format("None", "None"))

        if org == None and aux != None:
            row.append("enzyme\t{}\t{}\n".format(aux.string, "None"))

        if aux == None and org != None:
            row.append("enzyme\t{}\t{}\n".format("None", org.string))

    # Carriers

    carriers_elems = drug.find_all("carrier")

    for c_elem in carriers_elems:
        aux = c_elem.find("gene-name")
        org = c_elem.find("organism")

        if aux and aux.string != None:
            carriers.append(aux.string)

        if org != None and aux != None:
            row.append("carrier\t{}\t{}\n".format(aux.string, org.string))

        if org == None and aux == None:
            row.append("carrier\t{}\t{}\n".format("None", "None"))

        if org == None and aux != None:
            row.append("carrier\t{}\t{}\n".format(aux.string, "None"))

        if aux == None and org != None:
            row.append("carrier\t{}\t{}\n".format("None", org.string))

    # Synonyms

    synonyms_elems = drug.find_all("synonym")

    for syn_elem in synonyms_elems:
        syns.append(syn_elem.string)

    # atc_code

    atc_codes_elems = drug.find_all("atc-code")

    for atc_elem in atc_codes_elems:
        code = atc_elem["code"]
        level_elems = atc_elem.find_all("level")
        atcCodes[code] = "-"

        atc_aux = [code]

        for level_elem in level_elems:
            code = level_elem["code"]
            value = level_elem.string
            atcCodes[code] = value
            atc_aux.append(code)

        atc.append(",".join(atc_aux))

    # print(row)
    for element in row:
        output_gene.write(element)

    try:

        output.write(
            id.string.replace("\t", "").strip() + "\t" +
            name.replace("\t", "").strip() + "\t" +
            ";".join(set(syns)).replace("\t", "").strip() + "\t" +
            ";".join(atc).replace("\t", "").strip() + "\t" +
            smiles.replace("\t", "").strip() + "\t" +
            ";".join(set(transporters)).replace("\t", "").strip() + "\t" +
            ";".join(set(targets)).replace("\t", "").strip() + "\t" +
            ";".join(set(enzymes)).replace("\t", "").strip() + "\t" +
            ";".join(set(carriers)).replace("\t", "").strip() + "\t" +
            ";".join(set(groups)).strip() +
            "\n"
        )

    except TypeError as e:
        # print drug
        pprint(id)
        pprint(smiles)
        pprint(transporters)
        pprint(targets)
        pprint(groups)
        raise

# print("-----------")

atc_file = open("atc_codes.tsv", "w")

for key, value in atcCodes.items():
    atc_file.write(key.replace("\t", "").strip() + "\t" + value.replace("\t",
                                                                        "").strip() + "\n")

atc_file.close()
infile.close()
output_gene.close()
