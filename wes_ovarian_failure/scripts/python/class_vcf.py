# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 13:07:12 2019

@author: ihenarejos
"""
import os
import re
from collections import OrderedDict
from pprint import *  # for code revision

from class_variant import Variant

######################
# Description        #
######################

"""
class vcf
"""


# dir to class file
# os.chdir("/home/ihenarejos/workspace/projects/pof/scripts/python")
class Vcf:
    """_summary_"""

    def __init__(
        self,
        vcf_file,
        list_variants=None,
        header=None,
        patients_list=None,
        variants_id=None,
    ):

        """Creates an VCF-file like object"""

        self.vcf = vcf_file  # absPath saved to use in functions that
        # create new VCF from original
        self.header = header  # Header of VCF
        self.list_variants = []  # List of variants in VCF
        self.variants_id = ""  # Individual Id of the variant, from
        # combination of (chr,pos,ref,alt)
        self.study_population = []  # List of patients by experimental
        # design
        self.variants_id = variants_id  # List of individual IDs for each
        # variant
        self.sample_list = []  # List of samples analysed during Variant
        # Calling written to the VCF
        self.directory = []  # must be empty, only filled when its a file in
        # a dir
        self.parse()  # Calls parse method

    def parse(self):
        """_summary_"""
        if type(self.vcf) is str:  # if str is passed, assume is a vcf file
            # provided in path
            self.list_variants = [
                linea.rstrip("\n")
                for linea in open(self.vcf)
                if not linea.startswith("#")
            ]
            self.header = [
                linea.rstrip("\n") for linea in open(self.vcf) if linea.startswith("#")
            ]
            self.directory = os.path.splitext(self.vcf)[
                0
            ]  # to save child-VCF in same directory but different name

        elif type(self.vcf) is list:  # several functions returns list to
            # build a VCF
            self.list_variants = [
                linea.rstrip("\n") for linea in self.vcf if not linea.startswith("#")
            ]
            self.header = [
                linea.rstrip("\n") for linea in self.vcf if linea.startswith("#")
            ]

        self.simple_header = "".join([i for i in self.header if i.startswith("#CHROM")])

        # Make a list of the same length as list of variants where an
        # specific ID using several fields is assigned
        ids = []
        for variant_selected in self.list_variants:
            fields = variant_selected.split("\t")
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            variant_id = chrom + "_" + pos + "_" + ref + "_" + alt
            ids.append(variant_id)
        self.variants_id = ids

        # Make two lists of the same length where classification of samples
        # and samples names will be stored
        study_population = []
        sample_list = []
        # if type(self.vcf) is str:  # Get tab delimited fields of VCF
        fields_list = [
            line.rstrip("\t") for line in self.header if line.startswith("#CHROM")
        ]
        study_population = [element.split("\t") for element in fields_list]
        study_population = study_population[0][9:]

        for i in study_population:
            sample = i.strip("\n")
            patient = i.split("_")[0]
            if sample == "UNKNOWN_150":  # Specific case to rename a
                # sample wrongly named
                sample = "CONTROL_150"
            self.sample_list.append(sample)
            self.study_population.append(patient)

    def num_variants(self):
        return print(f"Number of variants in vcf file: {len(self.list_variants)}")

    def save_vcf(self, outfile):
        """Saves vcf object to file, writing line to line"""

        with open(outfile, "w") as out:
            for element in self.header:  # writes header lines
                out.write(str(element) + "\n")
            for element in self.list_variants:  # writes variants lines
                out.write(str(element) + "\n")

        out.close()

    @staticmethod
    def dict_changes_by_allele(variant_to_check):
        """Returns a dictionary where for each allele,
        the values are the changes to protein level annotated by snpEff.

        """

        #  initialize variables
        fields = variant_to_check.split("\t")
        info = fields[7]
        alt = fields[4]
        alleles = {}

        alleles_values = alt.split(",")
        for allele in alleles_values:
            # we add alleles as keys and values
            alleles[allele] = set()

        info_fields = info.split(";")
        for info_elem in info_fields:
            info_elem_splits = info_elem.split("=")
            info_key = info_elem_splits[0]

            if info_key == "ANN":
                ann = info_elem_splits[1]

        snp_eff_ann = ann.split(",")
        for snpEff_elem in snp_eff_ann:
            al = snpEff_elem.split("|")[0]
            change = snpEff_elem.split("|")[1]

            alleles[al].add(change)

        return alleles

    def multiallelic_transform(self) -> None:

        new_list_variants = []

        for line in self.header:
            new_list_variants.append(line)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            fields = variant_selected.split("\t")

            if re.findall(string=fields[4], pattern=","):
                variant = Variant(variant_selected, self.sample_list)
                multi_variants = variant.split_multiallelic_variant()

                for j in multi_variants:
                    new_list_variants.append(j)
            else:
                new_list_variants.append(variant_selected)

        # print(new_list_variants)
        dir = self.directory + "_multialellic_fix.vcf"
        Vcf(new_list_variants).save_vcf(dir)

        return new_list_variants

    def multi_allele_variants_filter(self):
        """ Creates a tuple with two vcf objects from a vcf file where the first one has no multi-allele variants \
            and the other is only the multi-allele variants. Also creates two VCF in the same directory as the father.

        :Returns:
            tuple vcf objects : no multi-variants, multi-variants
        """

        variants_no_multi = []
        variants_only_multi = []

        for linea in self.header:  # Rebuild header for child-VCF
            variants_no_multi.append(linea)
            variants_only_multi.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            fields = variant_selected.split("\t")

            if re.findall(string=fields[4], pattern=","):
                # If the delimiter for multiple alleles in same pos is found:
                variants_only_multi.append(
                    variant_selected
                )  # add those to another list
            else:  # if it only finds an allele on alt field:
                variants_no_multi.append(variant_selected)

        tuple_vcf = Vcf(variants_no_multi), Vcf(variants_only_multi)
        print(self.directory)
        save1 = self.directory + "_only_multi.vcf"
        save2 = self.directory + "_no_multi.vcf"
        (Vcf(variants_no_multi).save_vcf(save1))
        (Vcf(variants_only_multi).save_vcf(save2))

        return tuple_vcf

    def synonymous_variants_filter(self, change_list):
        """ Creates a tuple with two vcf objects from a vcf file where the first one has no multi-allele variants \
            and the other is only the multi-allele variants. Also creates two VCF in the same directory as the father.

        :Returns:
            tuple vcf objects : no multi-variants, multi-variants
        """

        variants_syn = []

        for line in self.header:  # Rebuild header for child-VCF
            variants_syn.append(line)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            fields = variant_selected.split("\t")
            syn = True

            if re.findall(string=fields[7], pattern="synonymous_variant"):
                # If synonymous_variant is is found:
                for change in change_list:
                    if re.findall(string=fields[7], pattern=change):
                        syn = False
                        break

                if syn:  # if there's a syn change and no moderate effect change
                    variants_syn.append(variant_selected)
                    continue

        print("File saved to", self.directory)
        save1 = self.directory + "_synonymous.vcf"
        (Vcf(variants_syn).save_vcf(save1))

        return Vcf(save1).num_variants()

    def filter_variants_by_gene_list(self, gene_list):
        """Filter variants if genes in list of genes are found as affected

        Parameters
        ----------
        :param gene_list:
            list of genes to check
        ...

        Returns
        -------
        object
            vcf object with variants that pass the filter
            child-vcf with the variants affected by said genes

        """
        if gene_list is list:
            pass
        else:
            try:
                gene_list = list(gene_list)
            except:
                print(
                    "Bad gene list input. Don't forget adding [] and "
                    + " separate by , "
                )
                exit()

        variants_info = []
        # add a way to get the list from an external file
        for linea in self.header:
            variants_info.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            genes = Variant(variant_selected, self.sample_list).get_genes_in_variant()

            # Check if genes in variant intersect with genelist
            matches = set(genes).intersection(set(gene_list))
            if len(matches) > 0:
                variants_info.append(variant_selected)

            # for gen in genes:
            #     if gen in gene_list:
            #         variants_info.append(variant_selected)
            #         break  # if condition is meet for a gene affected by the
            #         # variant, get out of the loop

        vcf_object = Vcf(variants_info)
        vcf_object.save_vcf(str(self.directory + "_filtered.vcf"))
        return vcf_object

    def filter_variants_by_gene_list_clinvar(self, gene_list):
        """Filter variants if genes in list of genes are found as affected

        Parameters
        ----------
        :param gene_list:
            list of genes to check
        ...

        Returns
        -------
        object
            vcf object with variants that pass the filter
            child-vcf with the variants affected by said genes

        """
        if gene_list is list:
            pass
        else:
            try:
                gene_list = list(gene_list)
            except:
                print(
                    "Bad gene list input. Don't forget adding [] and "
                    + " separate by , "
                )
                exit()

        variants_info = []
        # add a way to get the list from an external file
        for line in self.header:
            variants_info.append(line)

        for i in range(len(self.list_variants)):
            # print(self.variants_id[i])
            variant_selected = self.list_variants[i]
            # separate into fields
            genes = Variant(
                variant_selected, self.sample_list
            ).get_genes_in_variant_clinvar()

            # Check if genes in variant intersect with genelist
            matches = set(genes).intersection(set(gene_list))
            if len(matches) > 0:
                variants_info.append(variant_selected)

        vcf_object = Vcf(variants_info)
        # vcf_object.save_vcf(str(self.directory + "_filtered.vcf"))
        return vcf_object

    def filter_by_rs(self):
        """Creates two children-VCF files with the variants that have rs and the other that don't have rs

        Returns
        -------
        tuple :
            vcf with no rs, vcf with rs

        """

        variants_rs = []
        variants_no_rs = []

        for linea in self.header:  # add father header
            variants_rs.append(linea)
            variants_no_rs.append(linea)

        for i in range(len(self.list_variants)):
            v = self.list_variants[i]
            ids = (Variant(v, self.sample_list)).ids

            if ids != ["."]:
                variants_rs.append(v)
            else:
                variants_no_rs.append(v)

        Vcf(variants_no_rs).save_vcf(self.directory + "_no_rs.vcf")
        Vcf(variants_rs).save_vcf(self.directory + "_rs.vcf")

        return Vcf(variants_no_rs), Vcf(variants_rs)

    def filter_for_variants_igsr(self, annotated_term_phase3):
        """Will return a list of of lists, depending on how many annotations
        from different phases of 1000g
        it founds in the file provided
        Compatibility: phase3 in annotated form of provided and phase1 from
        snpSift and snpEff annotations \
        For phase 3 and 1

        Parameters
        ..........
        annotated_term_phase3 : str
            How phase3 is annotated in variants

        Returns ....... tuple A tuple of vcf objects in this order : "variants
        only in phase3,variants only in
        phase1, variants not in 1000g, variants with phase 3 annotations (that
        could have annotations of phase1)"

        ...
        """

        variants_only_phase3 = []
        variants_only_phase1 = []
        variants_not_in_1000g = []
        variants_1000g_phase3 = []

        for linea in self.header:  # add mother header
            variants_only_phase3.append(linea)
            variants_only_phase1.append(linea)
            variants_not_in_1000g.append(linea)
            variants_1000g_phase3.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            v_id = self.variants_id[i]
            # separate into fields
            fields = variant_selected.split("\t")
            info = fields[7]

            # phase3
            af = re.search("(%s)+\D+\d(\.\d*)*" % annotated_term_phase3, info)
            # phase1
            af2 = re.search("(dbNSFP_1000Gp1_AF)+\D+\d(\.\d*)*", info)

            if af is not None:
                variants_1000g_phase3.append(variant_selected)
            else:
                variants_not_in_phase3.append(variant_selected)

            if af is not None and af2 is None:
                variants_only_phase3.append(variant_selected)

            if af2 is not None and af is None:
                variants_only_phase1.append(variant_selected)

            if af is None and af2 is None:
                variants_not_in_1000g.append(variant_selected)

        vcf_phase3 = Vcf(variants_only_phase3)  # only phase3
        vcf_phase1 = Vcf(variants_only_phase1)  # only phase1
        vcf_not_in = Vcf(variants_not_in_1000g)  # no annotation in any phase
        vcf_1000g = Vcf(variants_1000g_phase3)  # with phase3 (and could
        # have phase1)
        vcf_not_phase3 = Vcf(variants_not_in_phase3)  # not in phase3

        return vcf_phase3, vcf_phase1, vcf_not_in, vcf_1000g, vcf_not_phase3

    def calculate_and_annotate_maf_igsr(self, to_add_term, annotated_term_phase3):

        variant_own_maf = []  # save calculated maf for every variant

        for linea in self.header:  # add mother header
            variant_own_maf.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            fields = variant_selected.split("\t")
            chrom = fields[0]
            pos = fields[1]
            ids = fields[2]
            ref = fields[3]
            alt = fields[4]
            quality = fields[5]
            filtered = fields[6]
            info = fields[7]
            gt = fields[8]
            genotypes = fields[9:]

            # phase 1 annotation
            maf_phase1 = 0
            af = re.search("(dbNSFP_1000Gp1_AF)+\D+\d(\.\d*)*", info)
            if af is not None:  # if af is found
                look_for_af1000 = af.group()
                af1000g = look_for_af1000.split("=")[1]  # get float value
                af = float(af1000g)
                if af <= 0.5:
                    maf_phase1 = af
                else:
                    maf_phase1 = 1 - af

            # phase 3 annotation
            maf_phase3 = ""
            af = re.search("(%s)+\D+\d(\.\d*)*" % annotated_term_phase3, info)
            if af is not None:  # if af is found
                look_for_af1000 = af.group()
                af1000g = look_for_af1000.split("=")[1]  # get float value
                # now, there could be multiple af values when studying multiple
                # -allele variants.
                # we choose to keep the lower AF value possible
                af1000g = af1000g.split(",")
                list_maf = []
                for af in af1000g:
                    if af != ".":
                        af = float(af)
                        if af <= 0.5:
                            maf = af
                        else:
                            maf = 1 - af
                        list_maf.append(maf)

                if len(list_maf) > 0:
                    maf_phase3 = min(list_maf)

            info = info + ";" + str(to_add_term) + "=" + str(maf_phase3)
            # modify and save variant
            variant_selected_to_append = self.__variant_to_rebuild(
                chrom, pos, ids, ref, alt, quality, filtered, info, gt, genotypes
            )
            variant_own_maf.append(variant_selected_to_append)

        variant_object = Vcf(variant_own_maf)
        return variant_object

    def maf_filter_igsr(self, maf_term, threshold=0.05):
        """filters vcf objects depending of maf_threshold. Maf must be
        annotated in info field. if not, it will assume that the variants is
        not in 1000g and add it to filtered list

        ...

        """
        variant_under_filter = []  # to save variants that pass the filter

        for linea in self.header:  # add mother header
            variant_under_filter.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            fields = variant_selected.split("\t")
            info = fields[7]
            maf_search = re.search("(%s)+\D+\d(\.\d*)*" % maf_term, info)

            if maf_search is not None:

                look_for_maf = maf_search.group()
                maf_value = look_for_maf.split("=")[1]
                maf = float(maf_value)

                if float(maf) < threshold:
                    variant_under_filter.append(variant_selected)

            elif maf_search is None:
                variant_under_filter.append(variant_selected)

        vcf_object = Vcf(variant_under_filter)
        return vcf_object

    def contingency_table_criteria_genes(
        self,
        outfile,
        target_change_list,
        allele_depth=4,
        depth_of_variant=10,
        genotype_quality=90,
        mutated_reads_over_total=0.2,
    ):
        """Creates a contingency table where the mutation needs to pass
        several filters to count the individual as affected, and saves it
        to a txt file, counting by genes

        ...
        ...

        Parameters ---------- genotype_quality : int Desired genotype
        quality from 1 to 99 (credibility of haplotype calculated by GATK)
        allele_depth : int Desired allele depth (number of reads that
        contain the variant) depth_of_variant : int Desired depth in the
        position of the variant (number of total reads in that position)
        mutated_reads_over_total : float Desired percentage of mutated reads
        over the total (alleDepth/depth_of_variant) outfile : txt file
        target_change_list: list Variant effects annotation to take into
        account when counting sample

        """

        # genes super dict; each value is going to be a gene and each
        # gene is going to be a dict where keys are samples and values
        # are number of affections
        genes_global = {}

        out = open(outfile, "w")
        string = "gene\t"
        for key in self.sample_list:
            # if key == "FOP9":
            string = string + key + "\t"
        out.writelines(string + "\n")

        depth_of_variant = depth_of_variant
        allele_depth = allele_depth
        genotype_quality = genotype_quality
        mutated_reads_over_total = mutated_reads_over_total

        # initialize variables
        variants = []
        for line in self.header:  # gets every line from header and put it
            # on a new list to recreate child VCF
            variants.append(line)

        for i in range(len(self.list_variants)):
            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            genotype = fields[9:]

            # genes in variant
            genes_info = Variant(
                variant_selected, self.sample_list
            ).set_genes_in_variant_changelist(target_change_list)
            # check genes affected by variant, add to global dict
            genes_in = genes_info[0]
            genes_changes = genes_info[1]
            for gene in genes_in:
                if gene not in genes_global:
                    genes_global[gene] = {sample: 0 for sample in self.sample_list}

            alleles = Vcf.dict_changes_by_allele(variant_selected)
            # put alternatives alleles in a new list that we
            # will use to access previous dictionary
            alleles_values = [fields[3]]
            list_alts = alternates.split(",")
            for alt in list_alts:
                alleles_values.append(alt)

            # Create a dict that holds the specific FORMAT for that variant,
            # saving each possible format as keys
            # in a dictionary with the values being an "index" that
            # we will need when checking every SAMPLE
            format_dict = {}
            format_value = 0

            for elem in fields[8].split(":"):
                format_dict[elem] = format_value
                format_value = format_value + 1

            if len(genotype) == len(self.study_population):

                for j in range(len(genotype)):
                    # initialize variables
                    alleles_splits = []
                    index_value_1 = 0
                    index_value_2 = 0
                    depth_of_variant_sample = ""
                    genes_affected_1 = ""
                    genes_affected_2 = ""
                    allele_depth_sample = ""

                    # Extract to what study group from
                    # the population the sample
                    sample_id = self.sample_list[j]

                    # extract the SAMPLE field for
                    # respective sample as well as genotype
                    sample = genotype[j]
                    sample_gt = sample.split(":")
                    haplotype = sample_gt[0]

                    # print(sample_id)

                    if haplotype == "0/0" or haplotype == "./.":
                        continue  # if haplotype is either
                        # same as genome of reference or data is missing, skip
                        # to next sample
                    else:  # if for that sample, there's presence of an
                        # alternative allele;
                        # check annotations on FORMAT field
                        if "AD" in format_dict:  # For certain variants,
                            # allele depth is omitted
                            ad = sample_gt[
                                format_dict["AD"]
                            ]  # and its the same than the number
                            # of reads obtained
                        else:
                            ad = sample_gt[format_dict["DP"]]  # so AD = DP

                        dp = sample_gt[format_dict["DP"]]
                        gq = sample_gt[format_dict["GQ"]]

                        # Extract alleles from SAMPLE
                        alleles_splits = haplotype.split(sep="/")
                        allele1 = alleles_splits[0]
                        allele2 = alleles_splits[1]
                        index_value_1 = int(allele1)
                        index_value_2 = int(allele2)
                        allele1_change_is_deleterious = False
                        allele2_change_is_deleterious = False

                        if index_value_1 != 0:  # integer assigned
                            # to alleles are the indexes for
                            if alleles_values[index_value_1] in genes_changes:
                                genes_affected_1 = genes_changes[
                                    alleles_values[index_value_1]
                                ]
                                allele1_change_is_deleterious = True

                        if index_value_2 != 0:
                            if alleles_values[index_value_2] in genes_changes:
                                genes_affected_2 = genes_changes[
                                    alleles_values[index_value_2]
                                ]
                                allele2_change_is_deleterious = True

                        if (
                            allele1_change_is_deleterious is False
                            and allele2_change_is_deleterious is False
                        ):
                            continue  # skip to next sample

                        if allele1_change_is_deleterious:  # check first allele
                            depth_values = str.split(ad, ",")
                            if (
                                len(depth_values) == 1 and depth_values[0] == "."
                            ):  # Total reads for a given allele
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_1]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                #  if all criteria is passed, count the sample
                                #  as affected for that gene in the group it
                                #  belongs;
                                # genes_affected_1 could be 1 gene or more (,)

                                for key in genes_affected_1:
                                    if key in genes_global:
                                        genes_global[key][sample_id] += 1

                                continue  # After counting this sample,
                                # don't evaluate second allele and jump to
                                # next sample

                        if allele2_change_is_deleterious:  # check second allele
                            depth_values = str.split(ad, ",")
                            if len(depth_values) == 1 and depth_values[0] == ".":
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_2]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):

                                for key in genes_affected_2:
                                    if key in genes_global:
                                        genes_global[key][sample_id] += 1

        for gene in genes_global:
            final_res = ""
            for sample_id in self.sample_list:
                res = genes_global[gene][sample_id]
                final_res = final_res + "\t" + str(res)

            out.writelines(gene + final_res + "\n")

        out.close()

        return genes_global

    def get_all_vep_consequences(self, consequence_list):

        consequences = set()
        for i in range(len(self.list_variants)):
            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            info = fields[7]

            info_fields = info.split(";")
            for info_elem in info_fields:
                info_elem_splits = info_elem.split("=")

                if "ANN" in info_elem_splits:

                    # save info without CSQ
                    indices = [i for i, s in enumerate(info_fields) if "ANN=" in s]

                    indx = indices[0]
                    data = info_fields[indx]

                    alleles = data.split(",")

                    for alelle in alleles:
                        changes = alelle.split("|")
                        changes = changes[1].split("&")
                        for change in changes:
                            if (
                                change not in consequences
                                and change in consequence_list
                            ):
                                consequences.add(change)
        return consequences

    def contingency_table_criteria(
        self,
        outfile,
        target_change_list,
        allele_depth=4,
        depth_of_variant=10,
        genotype_quality=90,
        mutated_reads_over_total=0.2,
    ):
        """Creates a contingency table where the mutation needs to pass
        several filters to count the individual as affected, and saves it
        to a txt file

        ...
        ...

        Parameters ---------- genotype_quality : int Desired genotype
        quality from 1 to 99 (credibility of haplotype calculated by GATK)
        allele_depth : int Desired allele depth (number of reads that
        contain the variant) depth_of_variant : int Desired depth in the
        position of the variant (number of total reads in that position)
        mutated_reads_over_total : float Desired percentage of mutated reads
        over the total (alleDepth/depth_of_variant) outfile : txt file
        target_change_list: list Variant effects annotation to take into
        account when counting sample

        """

        depth_of_variant = depth_of_variant
        allele_depth = allele_depth
        genotype_quality = genotype_quality
        mutated_reads_over_total = mutated_reads_over_total

        out = open(outfile, "w")

        out.write(
            "variant" + "\t" + "control" + "\t" + "foo" + "\t" + "fop" + "\t" + "cases"
        )

        # initialize variables
        variants = []
        for line in self.header:  # gets every line from header and put it
            # on a new list to recreate child VCF
            variants.append(line)

        for i in range(len(self.list_variants)):
            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            genotype = fields[9:]

            # dictionary with changes of protein consequence for each allele
            alleles = dict_changes_by_allele(variant_selected)
            # put alternatives alleles in a new list that we
            # will use to access previous dictionary
            alleles_values = [fields[3]]
            list_alts = alternates.split(",")
            for alt in list_alts:
                alleles_values.append(alt)

            # Create a dict that holds the specific FORMAT for that variant,
            # saving each possible format as keys
            # in a dictionary with the values being an "index" that
            # we will need when checking every SAMPLE
            format_dict = {}
            format_value = 0

            # Initialize where number of samples for each group
            # will be allocated
            controls = 0
            foo = 0
            fop = 0

            for elem in fields[8].split(":"):
                format_dict[elem] = format_value
                format_value = format_value + 1

            if len(genotype) == len(self.study_population):

                for j in range(len(genotype)):
                    # initialize variables
                    alleles_splits = []
                    index_value_1 = 0
                    index_value_2 = 0
                    depth_of_variant_sample = ""
                    allele_depth_sample = ""

                    # Extract to what study group from
                    # the population the sample
                    sample_group = self.study_population[j]

                    # extract the SAMPLE field for
                    # respective sample as well as genotype
                    sample = genotype[j]
                    sample_gt = sample.split(":")
                    haplotype = sample_gt[0]

                    if haplotype == "0/0" or haplotype == "./.":
                        continue  # if haplotype is either
                        # same as genome of reference or data is missing, skip
                        # to next sample
                    else:  # if for that sample, there's presence of an
                        # alternative allele;
                        # check annotations on FORMAT field
                        if "AD" in format_dict:  # For certain variants,
                            # allele depth is omitted
                            ad = sample_gt[
                                format_dict["AD"]
                            ]  # and its the same than the number
                            # of reads obtained
                        else:
                            ad = sample_gt[format_dict["DP"]]  # so AD = DP
                        dp = sample_gt[format_dict["DP"]]
                        gq = sample_gt[format_dict["GQ"]]

                        # Extract alleles from SAMPLE
                        alleles_splits = haplotype.split(sep="/")
                        allele1 = alleles_splits[0]
                        allele2 = alleles_splits[1]
                        index_value_1 = int(allele1)
                        index_value_2 = int(allele2)
                        allele1_change_is_deleterious = False
                        allele2_change_is_deleterious = False

                        if index_value_1 != 0:  # integer assigned
                            # to alleles are the indexes for
                            list_of_changes1 = alleles[alleles_values[index_value_1]]

                            for change_elem in list_of_changes1:
                                change_splits = change_elem.split(
                                    "&"
                                )  # there are combinations
                                # of changes with &,
                                # we have to split them to basic naming
                                # convention

                                for change in change_splits:
                                    if change in target_change_list:
                                        allele1_change_is_deleterious = True

                        if index_value_2 != 0:
                            list_of_changes2 = alleles[alleles_values[index_value_2]]
                            for change_elem in list_of_changes2:
                                change_splits = change_elem.split("&")

                                for change in change_splits:
                                    if change in target_change_list:
                                        allele2_change_is_deleterious = True

                        # Modification of SAMPLE field
                        if (
                            allele1_change_is_deleterious is False
                            and allele2_change_is_deleterious is False
                        ):
                            # if both alleles don't cause a targeted
                            # change in the protein:
                            for key in format_dict:

                                if key == "GT":
                                    sample_gt[
                                        format_dict[key]
                                    ] = "./."  # change genotype to
                                    # corresponding missing
                                else:
                                    sample_gt[
                                        format_dict[key]
                                    ] = "."  # change the rest of
                                    # fields to missing

                            genotype[j] = ":".join(
                                sample_gt
                            )  # Rebuild SAMPLE field for that
                            # specific sample
                            continue  # skip to next sample

                        if allele1_change_is_deleterious:  # check first allele
                            depth_values = str.split(ad, ",")
                            if (
                                len(depth_values) == 1 and depth_values[0] == "."
                            ):  # Total reads for a given allele
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_1]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                # if all criteria is passed, count the
                                # sample as affected on the group it belongs;
                                if sample_group == "CONTROL":
                                    controls += 1
                                elif sample_group == "FOO":
                                    foo += 1
                                elif sample_group == "FOP":
                                    fop += 1

                                continue  # After counting this sample,
                                # don't evaluate second allele and jump to
                                # next sample

                        if allele2_change_is_deleterious:  # check second
                            # allele
                            depth_values = str.split(ad, ",")
                            if len(depth_values) == 1 and depth_values[0] == ".":
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_2]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):

                                if sample_group == "CONTROL":
                                    controls += 1
                                elif sample_group == "FOO":
                                    foo += 1
                                elif sample_group == "FOP":
                                    fop += 1

            controls_in_variant = controls
            foo_in_variant = foo
            fop_in_variant = fop
            count_cases = foo_in_variant + fop_in_variant

            if controls_in_variant == 0 and count_cases == 0:

                continue  # skip to next variant, and don't write
                # information of current one on txt/tsv, since no sample
                # passed criteria
            else:

                out.write(
                    "\n"
                    + str(self.variants_id[i])
                    + "\t"
                    + str(controls_in_variant)
                    + "\t"
                    + str(foo_in_variant)
                    + "\t"
                    + str(fop_in_variant)
                    + "\t"
                    + str(count_cases)
                )
        out.close()

        return out

    def modify_sample_field_to_missing(
        self,
        target_change_list,
        allele_depth=4,
        depth_of_variant=10,
        genotype_quality=90,
        mutated_reads_over_total=0.2,
    ):
        """Creates a new VCF from the one loaded where SAMPLE field is
        changed to missing if FORMAT thresholds are not meet. New VCF must
        be saved using the corresponding function. ... ...

        :parameter: ---------- genotype_quality : int Desired genotype
        quality from 1 to 99 (credibility of genotype calculated by GATK)
        allele_depth : int Desired allele depth (number of reads that
        contain the variant) depth_of_variant : int Desired depth in the
        position of the variant (number of total reads in that position)
        mutated_reads_over_total : float Desired percentage of mutated reads
        over the total (alleDepth//depth_of_variant)

        :returns:
                vcf object
        """

        # initialize variables
        variants = []

        for line in self.header:  # gets every line from header and put it
            # on a new list to recreate child VCF
            variants.append(line)

        for i in range(len(self.list_variants)):

            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            genotype = fields[9:]

            # dictionary with changes of protein consequence for each allele
            alleles = dict_changes_by_allele(variant_selected)
            # put alternatives alleles in a new list that we will use to
            # access previous dictionary
            alleles_values = [fields[3]]
            list_alts = alternates.split(",")
            for alt in list_alts:
                alleles_values.append(alt)

            # Create a dict that holds the specific FORMAT for that variant,
            # saving each possible format as keys in a dictionary with the
            # values being an "index" that we will need when checking every
            # SAMPLE
            format_dict = {}
            format_value = 0

            for elem in fields[8].split(":"):
                format_dict[elem] = format_value
                format_value = format_value + 1

            if len(genotype) == len(self.study_population):

                for j in range(len(genotype)):

                    # initialize variables
                    alleles_splits = []
                    index_value_1 = 0
                    index_value_2 = 0
                    depth_of_variant_sample = ""
                    allele_depth_sample = ""

                    # extract the SAMPLE field for respective sample as well
                    # as genotype
                    sample = genotype[j]
                    sample_gt = sample.split(":")
                    haplotype = sample_gt[0]

                    if haplotype == "0/0" or haplotype == "./.":
                        continue  # if haplotype is either same as genome of
                        # reference or data is missing, skip
                        # to next sample
                    else:  # if for that sample, there's presence of an
                        # alternative allele;
                        # check annotations on FORMAT field
                        if "AD" in format_dict:  # For certain variants,
                            # allele depth is not written in FORMAT or SAMPLE
                            ad = sample_gt[
                                format_dict["AD"]
                            ]  # and its the same than the number
                            # of reads obtained
                        else:
                            ad = sample_gt[format_dict["DP"]]  # so AD = DP
                        dp = sample_gt[format_dict["DP"]]
                        gq = sample_gt[format_dict["GQ"]]

                        # Extract alleles from SAMPLE
                        alleles_splits = haplotype.split(sep="/")
                        allele1 = alleles_splits[0]
                        allele2 = alleles_splits[1]
                        index_value_1 = int(allele1)
                        index_value_2 = int(allele2)
                        allele1_change_is_deleterious = False
                        allele2_change_is_deleterious = False

                        if index_value_1 != 0:  # integer assigned to
                            # alleles are the indexes for
                            list_of_changes1 = alleles[alleles_values[index_value_1]]

                            for change_elem in list_of_changes1:
                                change_splits = change_elem.split(
                                    "&"
                                )  # there are combinations of changes

                                for change in change_splits:
                                    if change in target_change_list:
                                        allele1_change_is_deleterious = True

                        if index_value_2 != 0:
                            list_of_changes2 = alleles[alleles_values[index_value_2]]
                            for change_elem in list_of_changes2:
                                change_splits = change_elem.split("&")

                                for change in change_splits:
                                    if change in target_change_list:
                                        allele2_change_is_deleterious = True

                        # Modification of SAMPLE field
                        if (
                            allele1_change_is_deleterious is False
                            and allele2_change_is_deleterious is False
                        ):
                            # if both alleles don't cause a targeted change
                            # in the protein:
                            for key in format_dict:

                                if key == "GT":
                                    sample_gt[
                                        format_dict[key]
                                    ] = "./."  # change genotype to
                                    # corresponding missing
                                else:
                                    sample_gt[
                                        format_dict[key]
                                    ] = "."  # change the rest of
                                    # fields to missing

                            genotype[j] = ":".join(
                                sample_gt
                            )  # Rebuild SAMPLE field for that
                            # specific sample
                            continue  # skip to next sample

                        if allele1_change_is_deleterious:  # check first allele
                            depth_values = str.split(ad, ",")
                            if (
                                len(depth_values) == 1 and depth_values[0] == "."
                            ):  # Total reads for a given allele
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_1]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                continue  # skip to new sample
                            else:
                                if allele2_change_is_deleterious is False:
                                    # if allele 1 does not pass quality filters

                                    for key in format_dict:
                                        if key == "GT":
                                            sample_gt[format_dict[key]] = "./."
                                        else:
                                            sample_gt[format_dict[key]] = "."
                                    genotype[j] = ":".join(
                                        sample_gt
                                    )  # modify SAMPLE to missing

                        if allele2_change_is_deleterious:  # check second
                            # allele
                            depth_values = str.split(ad, ",")
                            if len(depth_values) == 1 and depth_values[0] == ".":
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_2]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                pass

                            else:
                                for key in format_dict:
                                    if key == "GT":
                                        sample_gt[format_dict[key]] = "./."
                                    else:
                                        sample_gt[format_dict[key]] = "."
                                genotype[j] = ":".join(sample_gt)

                modified_variant = self.__variant_to_rebuild(
                    fields[0],
                    fields[1],
                    fields[2],
                    fields[3],
                    fields[4],
                    fields[5],
                    fields[6],
                    fields[7],
                    fields[8],
                    genotype,
                )
                variants.append(modified_variant)
        Vcf(variants).save_vcf(self.directory + "_mod_to_missings_qual.vcf")

        return Vcf(variants)

    def modify_sample_field_to_missing_qual(
        self,
        allele_depth=4,
        depth_of_variant=10,
        genotype_quality=90,
        mutated_reads_over_total=0.2,
    ):
        """Creates a new VCF from the one loaded where SAMPLE field is
        changed to missing if FORMAT thresholds are not meet. New VCF must
        be saved using the corresponding function. ... ...

        :parameter: ---------- genotype_quality : int Desired genotype
        quality from 1 to 99 (credibility of genotype calculated by GATK)
        allele_depth : int Desired allele depth (number of reads that
        contain the variant) depth_of_variant : int Desired depth in the
        position of the variant (number of total reads in that position)
        mutated_reads_over_total : float Desired percentage of mutated reads
        over the total (alleDepth//depth_of_variant)

        :returns:
                vcf object
        """

        # initialize variables
        variants = []

        for line in self.header:  # gets every line from header and put it
            # on a new list to recreate child VCF
            variants.append(line)

        for i in range(len(self.list_variants)):

            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            genotype = fields[9:]

            # dictionary with changes of protein consequence for each allele
            alleles = dict_changes_by_allele(variant_selected)
            # put alternatives alleles in a new list that we will use to
            # access previous dictionary
            alleles_values = [fields[3]]
            list_alts = alternates.split(",")
            for alt in list_alts:
                alleles_values.append(alt)

            # Create a dict that holds the specific FORMAT for that variant,
            # saving each possible format as keys in a dictionary with the
            # values being an "index" that we will need when checking every
            # SAMPLE
            format_dict = {}
            format_value = 0

            for elem in fields[8].split(":"):
                format_dict[elem] = format_value
                format_value = format_value + 1

            if len(genotype) == len(self.study_population):

                for j in range(len(genotype)):

                    # initialize variables
                    alleles_splits = []
                    index_value_1 = 0
                    index_value_2 = 0
                    depth_of_variant_sample = ""
                    allele_depth_sample = ""

                    # extract the SAMPLE field for respective sample as well
                    # as genotype
                    sample = genotype[j]
                    sample_gt = sample.split(":")
                    haplotype = sample_gt[0]

                    if haplotype == "0/0" or haplotype == "./.":
                        continue  # if haplotype is either same as genome of
                        # reference or data is missing, skip
                        # to next sample
                    else:  # if for that sample, there's presence of an
                        # alternative allele;
                        # check annotations on FORMAT field
                        if "AD" in format_dict:  # For certain variants,
                            # allele depth is not written in FORMAT or SAMPLE
                            ad = sample_gt[
                                format_dict["AD"]
                            ]  # and its the same than the number
                            # of reads obtained
                        else:
                            ad = sample_gt[format_dict["DP"]]  # so AD = DP
                        dp = sample_gt[format_dict["DP"]]
                        gq = sample_gt[format_dict["GQ"]]

                        # Extract alleles from SAMPLE
                        alleles_splits = haplotype.split(sep="/")
                        allele1 = alleles_splits[0]
                        allele2 = alleles_splits[1]
                        index_value_1 = int(allele1)
                        index_value_2 = int(allele2)
                        allele1_change_is_deleterious = False
                        allele2_change_is_deleterious = False

                        if index_value_1 != 0:  # integer assigned to
                            # alleles are the indexes for

                            allele1_change_is_deleterious = True

                        if index_value_2 != 0:

                            allele2_change_is_deleterious = True

                        # Modification of SAMPLE field
                        if (
                            allele1_change_is_deleterious is False
                            and allele2_change_is_deleterious is False
                        ):
                            # if both alleles 0

                            for key in format_dict:

                                if key == "GT":
                                    sample_gt[
                                        format_dict[key]
                                    ] = "./."  # change genotype to
                                    # corresponding missing
                                else:
                                    sample_gt[
                                        format_dict[key]
                                    ] = "."  # change the rest of
                                    # fields to missing

                            genotype[j] = ":".join(
                                sample_gt
                            )  # Rebuild SAMPLE field for that
                            # specific sample
                            continue  # skip to next sample

                        if allele1_change_is_deleterious != 0:  # check first allele
                            depth_values = str.split(ad, ",")
                            if (
                                len(depth_values) == 1 and depth_values[0] == "."
                            ):  # Total reads for a given allele
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_1]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                continue  # skip to new sample
                            else:
                                if allele2_change_is_deleterious is False:
                                    # if allele 1 does not pass quality filters

                                    for key in format_dict:
                                        if key == "GT":
                                            sample_gt[format_dict[key]] = "./."
                                        else:
                                            sample_gt[format_dict[key]] = "."
                                    genotype[j] = ":".join(
                                        sample_gt
                                    )  # modify SAMPLE to missing

                        if allele2_change_is_deleterious:  # check second
                            # allele
                            depth_values = str.split(ad, ",")
                            if len(depth_values) == 1 and depth_values[0] == ".":
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_2]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                pass

                            else:
                                for key in format_dict:
                                    if key == "GT":
                                        sample_gt[format_dict[key]] = "./."
                                    else:
                                        sample_gt[format_dict[key]] = "."
                                genotype[j] = ":".join(sample_gt)

                modified_variant = self.__variant_to_rebuild(
                    fields[0],
                    fields[1],
                    fields[2],
                    fields[3],
                    fields[4],
                    fields[5],
                    fields[6],
                    fields[7],
                    fields[8],
                    genotype,
                )
                variants.append(modified_variant)
        Vcf(variants).save_vcf(self.directory + "_mod_to_missings")

        return Vcf(variants)

    def annotate_counts(
        self,
        target_change_list,
        allele_depth=4,
        depth_of_variant=10,
        genotype_quality=90,
        mutated_reads_over_total=0.2,
    ):
        """Annotates to info field counts of individuals. \
        For each individual that have at least one mutated allele, sum +1.
        It also filters by quality, so it omits counting samples below
        thresholds
        Parameters
        ----------

        Returns
        -------
        object
            vcf file object with variants annotated
        ...

        Parameters
        ----------
        :param genotype_quality:
        :param allele_depth:
        :param mutated_reads_over_total:
        :param depth_of_variant:
        :param target_change_list: effect at protein level
        :parameter: ---------- genotype_quality : int Desired genotype
        quality from 1 to 99 (credibility of genotype calculated by GATK)
        allele_depth : int Desired allele depth (number of reads that
        contain the variant) depth_of_variant : int Desired depth in the
        position of the variant (number of total reads in that position)
        mutated_reads_over_total : float Desired percentage of mutated reads
        over the total (alleDepth//depth_of_variant)

        """

        to_vcf = []  # list where header and variants will be saved
        for line in self.header:
            to_vcf.append(line)
        #  initialize variables
        depth_of_variant = depth_of_variant
        allele_depth = allele_depth
        genotype_quality = genotype_quality
        mutated_reads_over_total = mutated_reads_over_total

        for i in range(len(self.list_variants)):
            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            genotype = fields[9:]

            # dictionary with changes of protein consequence for each allele
            alleles = dict_changes_by_allele(variant_selected)
            # put alternatives alleles in a new list that we
            # will use to access previous dictionary
            alleles_values = [fields[3]]
            list_alts = alternates.split(",")
            for alt in list_alts:
                alleles_values.append(alt)

            # Create a dict that holds the specific FORMAT for that variant,
            # saving each possible format as keys
            # in a dictionary with the values being an "index" that
            # we will need when checking every SAMPLE
            format_dict = {}
            format_value = 0

            # Initialize where number of samples for each group
            # will be allocated
            controls = 0
            foo = 0
            fop = 0

            for elem in fields[8].split(":"):
                format_dict[elem] = format_value
                format_value = format_value + 1

            if len(genotype) == len(self.study_population):

                for j in range(len(genotype)):
                    # initialize variables
                    alleles_splits = []
                    index_value_1 = 0
                    index_value_2 = 0
                    depth_of_variant_sample = ""
                    allele_depth_sample = ""

                    # Extract to what study group from
                    # the population the sample
                    sample_group = self.study_population[j]

                    # extract the SAMPLE field for
                    # respective sample as well as genotype
                    sample = genotype[j]
                    sample_gt = sample.split(":")
                    haplotype = sample_gt[0]

                    if haplotype == "0/0" or haplotype == "./.":
                        continue  # if haplotype is either
                        # same as genome of reference or data is missing, skip
                        # to next sample
                    else:  # if for that sample, there's presence of an
                        # alternative allele;
                        # check annotations on FORMAT field
                        if "AD" in format_dict:  # For certain variants,
                            # allele depth is omitted
                            ad = sample_gt[
                                format_dict["AD"]
                            ]  # and its the same than the number
                            # of reads obtained
                        else:
                            ad = sample_gt[format_dict["DP"]]  # so AD = DP
                        dp = sample_gt[format_dict["DP"]]
                        gq = sample_gt[format_dict["GQ"]]

                        # Extract alleles from SAMPLE
                        alleles_splits = haplotype.split(sep="/")
                        allele1 = alleles_splits[0]
                        allele2 = alleles_splits[1]
                        index_value_1 = int(allele1)
                        index_value_2 = int(allele2)
                        allele1_change_is_deleterious = False
                        allele2_change_is_deleterious = False

                        if index_value_1 != 0:  # integer assigned
                            # to alleles are the indexes for
                            list_of_changes1 = alleles[alleles_values[index_value_1]]

                            for change_elem in list_of_changes1:
                                change_splits = change_elem.split(
                                    "&"
                                )  # there are combinations
                                # of changes with &,
                                # we have to split them to basic naming
                                # convention

                                for change in change_splits:
                                    if change in target_change_list:
                                        allele1_change_is_deleterious = True

                        if index_value_2 != 0:
                            list_of_changes2 = alleles[alleles_values[index_value_2]]
                            for change_elem in list_of_changes2:
                                change_splits = change_elem.split("&")

                                for change in change_splits:
                                    if change in target_change_list:
                                        allele2_change_is_deleterious = True

                        # Modification of SAMPLE field
                        if (
                            allele1_change_is_deleterious is False
                            and allele2_change_is_deleterious is False
                        ):
                            # if both alleles don't cause a targeted
                            # change in the protein:
                            for key in format_dict:

                                if key == "GT":
                                    sample_gt[
                                        format_dict[key]
                                    ] = "./."  # change genotype to
                                    # corresponding missing
                                else:
                                    sample_gt[
                                        format_dict[key]
                                    ] = "."  # change the rest of
                                    # fields to missing

                            genotype[j] = ":".join(
                                sample_gt
                            )  # Rebuild SAMPLE field for that
                            # specific sample
                            continue  # skip to next sample

                        if allele1_change_is_deleterious:  # check first allele
                            depth_values = str.split(ad, ",")
                            if (
                                len(depth_values) == 1 and depth_values[0] == "."
                            ):  # Total reads for a given allele
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_1]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                # if all criteria is passed, count the
                                # sample as affected on the group it belongs;
                                if sample_group == "CONTROL":
                                    controls += 1
                                elif sample_group == "FOO":
                                    foo += 1
                                elif sample_group == "FOP":
                                    fop += 1

                                continue  # After counting this sample,
                                # don't evaluate second allele and jump to
                                # next sample

                        if allele2_change_is_deleterious:  # check second
                            # allele
                            depth_values = str.split(ad, ",")
                            if len(depth_values) == 1 and depth_values[0] == ".":
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_2]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):

                                if sample_group == "CONTROL":
                                    controls += 1
                                elif sample_group == "FOO":
                                    foo += 1
                                elif sample_group == "FOP":
                                    fop += 1

            controls_in_variant = controls
            foo_in_variant = foo
            fop_in_variant = fop
            count_cases = foo_in_variant + fop_in_variant

            # Modify and save new variant
            chrom = fields[0]
            pos = fields[1]
            ids = fields[2]
            ref = fields[3]
            alt = fields[4]
            quality = fields[5]
            filtered = fields[6]
            info = fields[7]
            gt = fields[8]
            genotypes = fields[9:]

            new_info = (
                info
                + ";control_f="
                + str(controls_in_variant)
                + ";foo_f="
                + str(foo_in_variant)
                + ";fop_f="
                + str(fop_in_variant)
                + ";cases_f="
                + str(count_cases)
            )

            # __variant_to_rebuild
            variant_to_append = self.__variant_to_rebuild(
                chrom, pos, ids, ref, alt, quality, filtered, new_info, gt, genotypes
            )

            to_vcf.append(variant_to_append)

        vcf_object = Vcf(to_vcf)
        return vcf_object

    def annotate_counts_quality(
        self,
        allele_depth=4,
        depth_of_variant=10,
        genotype_quality=90,
        mutated_reads_over_total=0.2,
    ):
        """Annotates to info field counts of individuals. \
        For each individual that have at least one mutated allele, sum +1.
        Only filters by quality
        Parameters
        ----------

        Returns
        -------
        object
            vcf file object with variants annotated
        ...

        Parameters
        ----------
        :param genotype_quality:
        :param allele_depth:
        :param mutated_reads_over_total:
        :param depth_of_variant:
        :param target_change_list: effect at protein level
        :parameter: ---------- genotype_quality : int Desired genotype
        quality from 1 to 99 (credibility of genotype calculated by GATK)
        allele_depth : int Desired allele depth (number of reads that
        contain the variant) depth_of_variant : int Desired depth in the
        position of the variant (number of total reads in that position)
        mutated_reads_over_total : float Desired percentage of mutated reads
        over the total (alleDepth//depth_of_variant)

        """

        to_vcf = []  # list where header and variants will be saved
        for line in self.header:
            to_vcf.append(line)
        #  initialize variables
        depth_of_variant = depth_of_variant
        allele_depth = allele_depth
        genotype_quality = genotype_quality
        mutated_reads_over_total = mutated_reads_over_total

        for i in range(len(self.list_variants)):
            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            genotype = fields[9:]

            # dictionary with changes of protein consequence for each allele
            alleles = dict_changes_by_allele(variant_selected)
            # put alternatives alleles in a new list that we
            # will use to access previous dictionary
            alleles_values = [fields[3]]
            list_alts = alternates.split(",")
            for alt in list_alts:
                alleles_values.append(alt)

            # Create a dict that holds the specific FORMAT for that variant,
            # saving each possible format as keys
            # in a dictionary with the values being an "index" that
            # we will need when checking every SAMPLE
            format_dict = {}
            format_value = 0

            # Initialize where number of samples for each group
            # will be allocated
            controls = 0
            foo = 0
            fop = 0

            for elem in fields[8].split(":"):
                format_dict[elem] = format_value
                format_value = format_value + 1

            if len(genotype) == len(self.study_population):

                for j in range(len(genotype)):
                    # initialize variables
                    alleles_splits = []
                    index_value_1 = 0
                    index_value_2 = 0
                    depth_of_variant_sample = ""
                    allele_depth_sample = ""

                    # Extract to what study group from
                    # the population the sample
                    sample_group = self.study_population[j]

                    # extract the SAMPLE field for
                    # respective sample as well as genotype
                    sample = genotype[j]
                    sample_gt = sample.split(":")
                    haplotype = sample_gt[0]

                    if haplotype == "0/0" or haplotype == "./.":
                        continue  # if haplotype is either
                        # same as genome of reference or data is missing, skip
                        # to next sample
                    else:  # if for that sample, there's presence of an
                        # alternative allele;
                        # check annotations on FORMAT field
                        if "AD" in format_dict:  # For certain variants,
                            # allele depth is omitted
                            ad = sample_gt[
                                format_dict["AD"]
                            ]  # and its the same than the number
                            # of reads obtained
                        else:
                            ad = sample_gt[format_dict["DP"]]  # so AD = DP
                        dp = sample_gt[format_dict["DP"]]
                        gq = sample_gt[format_dict["GQ"]]

                        # Extract alleles from SAMPLE
                        alleles_splits = haplotype.split(sep="/")
                        allele1 = alleles_splits[0]
                        allele2 = alleles_splits[1]
                        index_value_1 = int(allele1)
                        index_value_2 = int(allele2)
                        allele1_is_valid = False
                        allele2_is_valid = False

                        if index_value_1 != 0:  # integer assigned
                            # to alleles are the indexes for
                            allele1_is_valid = True

                        if index_value_2 != 0:
                            allele2_is_valid = True

                        if allele1_is_valid:  # check first allele
                            depth_values = str.split(ad, ",")
                            if (
                                len(depth_values) == 1 and depth_values[0] == "."
                            ):  # Total reads for a given allele
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_1]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):
                                # if all criteria is passed, count the
                                # sample as affected on the group it belongs;
                                if sample_group == "CONTROL":
                                    controls += 1
                                elif sample_group == "FOO":
                                    foo += 1
                                elif sample_group == "FOP":
                                    fop += 1

                                continue  # After counting this sample,
                                # don't evaluate second allele and jump to
                                # next sample

                        if allele2_is_valid:  # check second
                            # allele
                            depth_values = str.split(ad, ",")
                            if len(depth_values) == 1 and depth_values[0] == ".":
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_2]
                            depth_of_variant_sample = dp  # DP

                            if (
                                int(gq) >= genotype_quality
                                and int(allele_depth_sample) >= allele_depth
                                and int(depth_of_variant_sample) >= depth_of_variant
                                and (int(allele_depth_sample))
                                / (int(depth_of_variant_sample))
                                >= mutated_reads_over_total
                            ):

                                if sample_group == "CONTROL":
                                    controls += 1
                                elif sample_group == "FOO":
                                    foo += 1
                                elif sample_group == "FOP":
                                    fop += 1

            controls_in_variant = controls
            foo_in_variant = foo
            fop_in_variant = fop
            count_cases = foo_in_variant + fop_in_variant

            # Modify and save new variant
            chrom = fields[0]
            pos = fields[1]
            ids = fields[2]
            ref = fields[3]
            alt = fields[4]
            quality = fields[5]
            filtered = fields[6]
            info = fields[7]
            gt = fields[8]
            genotypes = fields[9:]

            new_info = (
                info
                + ";control_f="
                + str(controls_in_variant)
                + ";foo_f="
                + str(foo_in_variant)
                + ";fop_f="
                + str(fop_in_variant)
                + ";cases_f="
                + str(count_cases)
            )

            # __variant_to_rebuild
            variant_to_append = self.__variant_to_rebuild(
                chrom, pos, ids, ref, alt, quality, filtered, new_info, gt, genotypes
            )

            to_vcf.append(variant_to_append)

        vcf_object = Vcf(to_vcf)
        return vcf_object

    def filter_by_counts(self, control_threshold, foo_threshold, fop_threshold):
        """filters vcf objects depending of samples threshold %.
        counts must be annotated in info field

        Parameters
        ----------
        control_threshold : int
            Number of control samples you want to find variants below
        foo_threshold : int
            Number of control samples you want to find variants over
        fop_threshold : int
            Number of control samples you want to find variants over
        ...

        """
        variant_under_filter = []  # to save variants that pass the filter

        for linea in self.header:  # add mother header
            variant_under_filter.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            fields = variant_selected.split("\t")
            info = fields[7]
            control_value = 0
            foo_value = 0
            fop_value = 0

            control_search = re.search("(control_f=)+\d*", info)
            look_for_control = control_search.group()
            control_value = float(look_for_control.split("=")[1])

            foo_search = re.search("(foo_f=)+\d*", info)
            look_for_foo = foo_search.group()
            foo_value = float(look_for_foo.split("=")[1])

            fop_search = re.search("(fop_f=)+\d*", info)
            look_for_fop = fop_search.group()
            fop_value = float(look_for_fop.split("=")[1])

            if (
                control_value <= control_threshold
                and foo_value >= foo_threshold
                and fop_value >= fop_threshold
            ):
                variant_under_filter.append(variant_selected)

        vcf_object = Vcf(variant_under_filter)
        return vcf_object

    # TODO DO BETTER
    def calculate_annotate_own_maf(self, annotation_term):

        """Creates a list of variants where for each variant population size
         maf is calculated and annotated at the end of info field

        ...

        Parameters
        ----------
        annotation_term : str
            String to annotate the calculated maf with
        """
        list_variants = []  # were we will save maf results
        variant_own_maf = []  # save calculated maf for every variant
        for linea in self.header:  # add mother header
            variant_own_maf.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            v_id = self.variants_id[i]
            # separate into fields
            fields = variant_selected.split("\t")
            chrom = fields[0]
            pos = fields[1]
            ids = fields[2]
            ref = fields[3]
            alt = fields[4]
            quality = fields[5]
            filtered = fields[6]
            info = fields[7]
            gt = fields[8]
            genotypes = fields[9:]

            # variables to calculate own maf
            maf = 0
            af = 0
            missing = 0
            allele_count = 0

            for j in genotypes:
                sample_gt = j.split(":")
                haplotype = sample_gt[0]
                allele_1 = haplotype.split(sep="/")[0]
                allele_2 = haplotype.split(sep="/")[1]

                # Check both alleles
                if allele_1 == ".":
                    missing += 1
                else:
                    if int(allele_1) >= 1:
                        allele_count += 1

                if allele_2 == ".":
                    missing += 1
                else:
                    if int(allele_2) >= 1:
                        allele_count += 1

            #  Calculate af and maf
            af = allele_count / (len(self.study_population) * 2)
            if af >= 0.5:
                maf = 1 - af
            else:
                maf = af

            info = info + ";" + str(annotation_term) + "=" + str(maf)

            variant_selected = self.__variant_to_rebuild(
                chrom, pos, ids, ref, alt, quality, filtered, info, gt, genotypes
            )

            variant_selected_to_append = Variant.build_variant(variant_selected)
            variant_own_maf.append(variant_selected_to_append)

        variant_object = Vcf(variant_own_maf)
        return variant_object

    def evaluate_variants(self, missing_threshold):
        """Evaluate variants by filtering depending of threshold provided
        for missing values

        ...

        """

        variant_under_filter = []

        for linea in self.header:
            variant_under_filter.append(linea)

        for i in range(len(self.list_variants)):
            variant_selected = Variant(self.list_variants[i], self.sample_list)

            if variant_selected.analyze_missing(missing_threshold):
                # Will return True if variant passes threshold
                variant_under_filter.append(variant_selected)

        variant_object = Vcf(variant_under_filter)
        return variant_object

    def get_snpsift_case_control(self, output):
        """Gets scores from each model of SnpSift

        :param output: file to write results
        :return: list
        """
        snpsift_annotations = []

        out = open(output, "w")
        out.write(
            "variant\t"
            + "CC_DOM"
            + "\t"
            + "CC_REC"
            + "\t"
            + "CC_ALL"
            + "\t"
            + "CC_GENO"
            + "\t"
            + "CC_TREND"
        )

        for i in range(len(self.list_variants)):
            variant_selected = Variant(self.list_variants[i], self.sample_list)
            list_scores = variant_selected.get_snpsift_case_control()
            res = ""
            res = "\t".join(list_scores)  # this method doesnt work with floats
            res = "\t".join(map(str, list_scores))
            out.write(res)

        out.close()

        return snpsift_annotations

    def get_groups_proportions(self, outfile):
        """Makes use of calculate_proportions from class variant to get
        proportions for each group.
        """

        results = ""
        with open(outfile, "w") as out:
            out.writelines(
                "v\t" + "ctr_ref\tctr_het\tctr_hom\tctr_miss\tfoo_ref"
                "\tfoo_het\tfoo_hom\tfoo_miss\tfop_ref\tfop_het"
                "\tfop_hom\tfop_miss"
            )
            for i in range(len(self.list_variants)):
                results = Variant(
                    self.list_variants[i], self.sample_list
                ).calculate_proportions()

                out.writelines("\n" + str(self.variants_id[i]) + "\t" + results)

        out.close()

    def get_cc_counts(self, outfile):
        """Returns count for each genotype for cases vs controls"""

        results = ""

        with open(outfile, "w") as out:
            out.writelines(
                "v\t" + "ctr_ref\tctr_het\tctr_hom\tcase_ref\tcase_het" "\tcase_hom"
            )
            for i in range(len(self.list_variants)):
                results = Variant(
                    self.list_variants[i], self.sample_list
                ).calculate_counts_cc()

                out.writelines("\n" + str(self.variants_id[i]) + "\t" + results)

        out.close()

    def get_haplotype_counts(self, outfile):
        """Returns count for each genotype for cases vs controls"""

        results = ""

        with open(outfile, "w") as out:
            out.writelines("v\t" + "h00\th01\th11")
            for i in range(len(self.list_variants)):
                results = Variant(
                    self.list_variants[i], self.sample_list
                ).calculate_counts_haplotypes()

                out.writelines("\n" + str(self.variants_id[i]) + "\t" + results)

        out.close()

    def mean_coverage_file(self):
        """
        Mean DP value across all VCF
        """

        result = 0

        for i in range(len(self.list_variants)):
            result = result + int(
                (Variant(self.list_variants[i], self.sample_list)).get_dp()
            )

        return result / len(self.list_variants)

    def mean_coverage_sd_samples(self, outfile):
        """
        based on DP of all variants selected
        mean is computed for all samples where variant was sequenced
        meaning, missing values are not considered, since it could mean that
        sample wasn't sequenced in that position
        """

        result = 0
        with open(outfile, "w") as out:
            out.writelines("v\tcoverage\n")

            for i in range(len(self.list_variants)):
                result = int(
                    (
                        Variant(self.list_variants[i], self.sample_list)
                    ).get_mean_dp_samples()
                )
                sd = str(
                    (
                        Variant(self.list_variants[i], self.sample_list)
                    ).get_sd_dp_samples()
                )
                out.writelines(
                    (self.variants_id[i] + "\t" + str(result) + "\t" + sd + "\n")
                )

        out.close()

    def check_all_samples_presence_absence(self, outfile):
        """Creatures a discrete feature table for WEKA

        :param outfile: dir and name of file to write results
        """

        results_dict = {}

        for i in self.sample_list:
            results_dict[i] = 0
        # Force and save order of keys in dictionary
        sorted(results_dict.items(), key=lambda x: x[0])
        results_dict = OrderedDict(sorted(results_dict.items()))

        with open(outfile, "w") as out:
            string = "v\t"
            for key in results_dict:
                string = string + key + "\t"
            out.writelines(string)

            for i in range(len(self.list_variants)):
                results = Variant(
                    self.list_variants[i], self.sample_list
                ).genotype_samples_variant()
                row = ""

                for key in results:
                    value = results[key]  # get the value from variant dict
                    results_dict[key] = str(results_dict[key]) + value
                    row = row + str(value) + "\t"

                out.writelines("\n" + str(self.variants_id[i]) + "\t" + row)

        out.close()

    def samples_affected(self, outfile):
        """Extract samples affected in each variant

        :param outfile: Dir and name of file to write results
        """

        with open(outfile, "w") as out:
            out.writelines("id\tsample")

            for i in range(len(self.list_variants)):
                genes = ""
                genes = Variant(
                    self.list_variants[i], self.sample_list
                ).set_samples_affected()

                for gene in genes:
                    out.writelines("\n" + str(self.variants_id[i]) + "\t" + str(gene))

            out.close()

    def calculate_normalized_positions_mutated(
        self, bed_file, synonym_file, gene_list, changelist=None, sum_all=False
    ):
        """
        check how many variants are affecting each gene; needs a bed file
        with the size of each gene can work with a list of biological
        annotations of variants provided by user.

        :param sum_all indicates if user wants to store results individually
        for each gene or just get the cumulative sum of all
        """

        log_num = 1
        synonym_dict = {}
        synonym = open(synonym_file, "r")
        lines = synonym.readlines()

        # create a dictionary where each key has the original genename and
        # correct synonym is stored as value
        for line in lines:
            line = line.strip("\n")
            splits = line.split("\t")
            synonym_dict[splits[0]] = splits[1]

        print("LOG[", log_num, "]: synonyms file processed")
        log_num += 1
        print("LOG[", log_num, "]: Bed file being processed")
        log_num += 1

        c = 0  # to control found genes
        sizes = {}
        not_found = []
        bed = open(bed_file, "r")
        lines = bed.readlines()

        for line in lines:
            if lines[0]:
                pass
            splits = (line.strip("\n")).split("\t")
            gene = splits[0]
            size = splits[1]
            sizes[gene] = size

        print("LOG[", log_num, "]: Bed file processed")
        log_num += 1
        print(
            "LOG[",
            log_num,
            "]: Processing VCF file to check genes affected in each variant",
        )
        log_num += 1

        dict_genes_variant = {}
        for i in range(len(self.list_variants)):
            set_of_genes_in = Variant(
                self.list_variants[i], self.sample_list
            ).set_genes_in_variant_changelist(changelist)

            for gene in set_of_genes_in:
                if gene in synonym_dict:
                    # if gene is found as key in dictionary of synonyms,
                    # its value is modified so its the correct synonym to
                    # use for the next step
                    gene = synonym_dict[gene]

                if gene not in dict_genes_variant:
                    dict_genes_variant[gene] = 0
                    dict_genes_variant[gene] += 1

                else:
                    dict_genes_variant[gene] += 1

        print("LOG[", log_num, "]: VCF file processed")
        log_num += 1
        print(
            "LOG[", log_num, "]: Calculating normalized rate of sites mutated by gene"
        )
        log_num += 1

        if not sum_all:
            # Dict of affection of genes by affections/gene_size
            dict_results_affection = {}

            for key in dict_genes_variant:  # dict obtained by combining
                # snpEff data and correct synonyms

                if key in sizes:
                    gene_size = int(sizes[key])
                    dict_results_affection[key] = (
                        int(dict_genes_variant[key]) / gene_size
                    )
                else:
                    c += 1  # check how many are not found
                    not_found.append(key)

            return (
                dict_results_affection,
                not_found,
            )  # returns dict of found and list of not found
            # genes

        else:
            dict_results_affection = {}
            sum_size = 0
            all_sites_mutated = 0
            for key in dict_genes_variant:
                if key in sizes:
                    if key in gene_list:
                        gene_size = int(sizes[key])
                        sum_size += gene_size
                        dict_results_affection[key] = int(dict_genes_variant[key])

            for key in dict_results_affection:
                sites_mutated_gene = dict_results_affection[key]
                all_sites_mutated += sites_mutated_gene

            return all_sites_mutated, sum_size

    def return_snpeff_field(self, outfile, field):
        """Returns snpeff annotation desired for each variant

        :param field: str annotation from snpeff to return
        :param outfile: dir + file where file with results will be written
        :return: txt file
        """

        file = open(outfile, "w")
        file.writelines("id" + "\t" + str(field))

        for i in range(len(self.list_variants)):
            variant_id = self.variants_id[i]
            bio_type = (
                Variant(self.list_variants[i], self.sample_list)
            ).get_snpeff_field(str(field))

            result = str(variant_id) + "\t" + str(bio_type)

            file.writelines("\n" + result)
        file.close()

    def samples_affected_by_variant(self, variant_id):

        for i in range(len(self.list_variants)):

            if self.variants_id[i] == variant_id:
                print("[LOG #] : Found variant ")
                v = self.list_variants[i]

                return (Variant(v, self.sample_list)).set_samples_affected()

    def which_variant_affect_same_samples(self, pattern_length):
        """Returns number of samples affecting a sample

        :param pattern_length: number of samples that must be affected by
        variant :return:
        """

        d_results = {}
        import collections

        print("# LOG: Beginning processing variants")
        for i in range(len(self.list_variants)):

            v = self.list_variants[i]

            samples = list((Variant(v, self.sample_list)).set_samples_affected())

            if len(samples) == int(pattern_length):
                string_of_samples = ",".join(samples)

                id = self.variants_id[i]

                d_results[id] = string_of_samples

        print("# LOG: Finished processing variants")

        target = collections.defaultdict(tuple)

        for key in d_results:
            target[d_results[key]] += (key,)

        # pp.pprint(target)

        res = []

        for i, j in target.items():

            if len(j) > 1:
                print(i, j)

                res.append(j)

        return res

    # # TO REVIEW
    # def variants_of_gene(self, changelist, outfile):
    #     """ Extract variants affecting a gene specified.

    #     :return: dict where keys are the genes of the list and values are id
    #     of variants that affect said genes
    #     """

    #     variants_aff = {}
    #     for i in range(len(self.list_variants)):
    #         v = self.list_variants[i]
    #         vid = self.variants_id[i]
    #         counts = (
    #             Variant(v, self.sample_list)).parse_counts()  # returns a dict
    #         res = ""

    #         for key, value in counts.items():
    #             res = res + key + ":" + value + " "
    #         change_type = (
    #             Variant(v, self.sample_list)).get_snpeff_field(
    #             "Annotation")

    #         samples_affected = list(
    #             (Variant(v, self.sample_list)).set_samples_affected())
    #         string_samples = ""

    #         for element in samples_affected:
    #             string_samples = string_samples + "," + str(element)

    #         information = str(res) + str(change_type) + " " + str(
    #             string_samples)

    #         for gene in list(
    #                 (Variant(v,
    #                          self.sample_list)).set_genes_in_variant_changelist
    #                     (changelist)):

    #             tuple_data = (information, vid)
    #             if gene not in variants_aff:
    #                 variants_aff[gene] = []
    #                 variants_aff[gene].append(tuple_data)
    #             else:
    #                 variants_aff[gene].append(tuple_data)

    #     file = open(outfile, "w")
    #     file.write("gene\tfop\tcontrol\tfoo\tchange\tsample\tvariant\n")

    #     for key, value in variants_aff.items():
    #         for info in value:
    #             res = str(info).replace("\t", "")
    #             res = ((((
    #                 ((str(res)).replace("(", "")).replace(")", ""))).replace(
    #                 "'", "")).replace(" ", "\t"))

    #             file.write(str(key) + "\t" + res + "\n")
    #     file.close()

    #     return variants_aff

    def get_variants_using_id(self, id_file=None, id_list_set=None):
        """Uses id (chr_pos_ref_alt) of the variants to look for then in
        father vcf file.

        Parameters
        ----------
        id_file : file
            a file where each line is a variant id
        id_list : list
            list of id variants (chrom_pos_ref_alt)
        ...
        Returns
        -------
        object :
            vcf object with the variants found using the ids

        """

        # Check if args are default
        ext_ids = ""
        if id_file is not None:
            ext_ids = [linea.rstrip("\n") for linea in open(id_file, "r")]
        if id_list_set is not None:
            ext_ids = id_list_set

        # initialize variable
        header = [i for i in self.header]

        selection = [i for i in self.return_id_variant_all() if i in ext_ids]

        filtered = [
            i
            for i in self.list_variants
            if Variant(i, self.sample_list).return_id_variant(sep="_") in selection
        ]

        vcf_object = Vcf([*header, *filtered])
        return vcf_object

    def return_rs_all(self, outfile):
        """Extract rs of each variants and save them to a new file

        :param outfile: Txt file to write results
        """
        results = ""

        with open(outfile, "w") as out:
            out.writelines("v\trs")

            for i in range(len(self.list_variants)):
                results = Variant(self.list_variants[i], self.sample_list).return_rs()
                out.writelines("\n" + str(self.variants_id[i]) + "\t" + results)

            out.close()

    def return_mean_ad_coverage(self, outfile):

        """Creates a contingency table where the mutation needs to pass
        several filters to count the individual as affected, and saves it
        to a txt file

        ...
        ...

        Parameters ---------- genotype_quality : int Desired genotype
        quality from 1 to 99 (credibility of haplotype calculated by GATK)
        allele_depth : int Desired allele depth (number of reads that
        contain the variant) depth_of_variant : int Desired depth in the
        position of the variant (number of total reads in that position)
        mutated_reads_over_total : float Desired percentage of mutated reads
        over the total (alleDepth/depth_of_variant) outfile : txt file
        target_change_list: list Variant effects annotation to take into
        account when counting sample

        """
        out = open(outfile, "w")

        out.write("variant" + "\t" + "mean_ad\n")

        # initialize variables
        variants = []
        for line in self.header:  # gets every line from header and put it
            # on a new list to recreate child VCF
            variants.append(line)

        for i in range(len(self.list_variants)):

            # results of variant ad
            res = 0
            # aff
            aff = 0

            # extract fields to work with
            variant_selected = self.list_variants[i]
            fields = variant_selected.split("\t")
            alternates = fields[4]
            genotype = fields[9:]

            # dictionary with changes of protein consequence for each allele
            alleles = dict_changes_by_allele(variant_selected)
            # put alternatives alleles in a new list that we
            # will use to access previous dictionary
            alleles_values = [fields[3]]
            list_alts = alternates.split(",")
            for alt in list_alts:
                alleles_values.append(alt)

            # Create a dict that holds the specific FORMAT for that variant,
            # saving each possible format as keys
            # in a dictionary with the values being an "index" that
            # we will need when checking every SAMPLE
            format_dict = {}
            format_value = 0

            # Initialize where number of samples for each group
            # will be allocated
            controls = 0
            foo = 0
            fop = 0

            for elem in fields[8].split(":"):
                format_dict[elem] = format_value
                format_value = format_value + 1

            if len(genotype) == len(self.study_population):

                for j in range(len(genotype)):
                    # initialize variables
                    alleles_splits = []
                    index_value_1 = 0
                    index_value_2 = 0
                    depth_of_variant_sample = ""
                    allele_depth_sample = ""

                    # Extract to what study group from
                    # the population the sample
                    sample_group = self.study_population[j]

                    # extract the SAMPLE field for
                    # respective sample as well as genotype
                    sample = genotype[j]
                    sample_gt = sample.split(":")
                    haplotype = sample_gt[0]

                    if haplotype == "0/0" or haplotype == "./.":
                        continue  # if haplotype is either
                        # same as genome of reference or data is missing, skip
                        # to next sample
                    else:  # if for that sample, there's presence of an
                        # alternative allele;
                        # check annotations on FORMAT field
                        if "AD" in format_dict:  # For certain variants,
                            # allele depth is omitted
                            ad = sample_gt[
                                format_dict["AD"]
                            ]  # and its the same than the number
                            # of reads obtained
                        else:
                            ad = sample_gt[format_dict["DP"]]  # so AD = DP
                        dp = sample_gt[format_dict["DP"]]
                        gq = sample_gt[format_dict["GQ"]]

                        # Extract alleles from SAMPLE
                        aff += 1  # patient has variant
                        alleles_splits = haplotype.split(sep="/")
                        allele1 = alleles_splits[0]
                        allele2 = alleles_splits[1]
                        index_value_1 = int(allele1)
                        index_value_2 = int(allele2)
                        allele1_change_is_deleterious = False
                        allele2_change_is_deleterious = False

                        if index_value_1 != 0:  # integer assigned
                            # to alleles are the indexes for
                            allele1_change_is_deleterious = True

                        if index_value_2 != 0:
                            allele2_change_is_deleterious = True

                        if allele1_change_is_deleterious:  # check first allele
                            depth_values = str.split(ad, ",")
                            if (
                                len(depth_values) == 1 and depth_values[0] == "."
                            ):  # Total reads for a given allele
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_1]
                            depth_of_variant_sample = dp  # DP

                            res = res + int(allele_depth_sample)

                        if allele2_change_is_deleterious:  # check second allele
                            depth_values = str.split(ad, ",")
                            if len(depth_values) == 1 and depth_values[0] == ".":
                                allele_depth_sample = dp
                            elif len(depth_values) == 1:
                                allele_depth_sample = dp
                            elif len(depth_values) == len(alleles_values):
                                allele_depth_sample = depth_values[index_value_2]
                            depth_of_variant_sample = dp  # DP

                            if index_value_2 != index_value_1:
                                res = res + int(allele_depth_sample)
            print(aff, res)
            out.write(self.variants_id[i] + "\t" + str(res // aff) + "\n")

        out.close()

        return out

    def return_genes(self, outfile="", list=False):
        """Returns genes affected by each variant

        :param outfile: Txt file to write results
        """

        if list == True:
            results = set()
            for i in range(len(self.list_variants)):
                genes = Variant(
                    self.list_variants[i], self.sample_list
                ).get_genes_in_variant()
                for j in genes:
                    results.add(j)

            return results

        else:

            out = open(outfile, "w")
            out.writelines("id\tgene")

            for i in range(len(self.list_variants)):
                results = ""

                results = (
                    str(self.variants_id[i])
                    + "\t"
                    + str(
                        ",".join(
                            list(
                                Variant(
                                    self.list_variants[i], self.sample_list
                                ).get_genes_in_variant()
                            )
                        )
                    )
                )

                out.writelines("\n" + results)
            out.close()

    def modify_genes(self, gene_list, syns_list):

        for i in range(len(self.list_variants)):
            variant_selected = self.list_variants[i]
            # separate into fields
            fields = variant_selected.split("\t")
            chrom = fields[0]
            pos = fields[1]
            ids = fields[2]
            ref = fields[3]
            alt = fields[4]
            quality = fields[5]
            filtered = fields[6]
            info = fields[7]
            gt = fields[8]
            genotypes = fields[9:]

            genes_in = Variant(
                self.list_variants[i], self.sample_list
            ).get_genes_in_variant()

            list_genes = genes_in.split(",")

            for i in list_genes:
                if i in gene_list:
                    index_value = gene_list.index(i)
                    syn = syns_list(index_value)
                    res = ""

                else:
                    continue

    def return_bioinfo(self, outfile):
        """Returns several information from each variant, like id, effect
        on protein, predictors scores...

        :param outfile: Txt file to write results
        """

        out = open(outfile, "w")
        out.writelines(
            "id\tchrom\tpos\tref\talt\trs\tgene\taa\tbioEff\tCADD"
            "\tMutTaster\tcounts\tnumGenesAff"
        )

        for i in range(len(self.list_variants)):
            results = ""
            results = (
                str(self.variants_id[i])
                + "\t"
                + str(
                    # Calls variant class to extract required information
                    Variant(self.list_variants[i], self.sample_list).return_id_variant()
                )
                + "\t"
                + ",".join(
                    list(
                        Variant(
                            self.list_variants[i], self.sample_list
                        ).get_genes_in_variant()
                    )
                )
                + "\t"
                + (
                    Variant(self.list_variants[i], self.sample_list).get_snpeff_field(
                        "HGVS.p"
                    )
                )
                + "\t"
                + (
                    Variant(self.list_variants[i], self.sample_list).get_snpeff_field(
                        "Annotation"
                    )
                )
                + "\t"
                + str(
                    Variant(self.list_variants[i], self.sample_list).return_cadd_score()
                )
                + "\t"
                + str(
                    Variant(
                        self.list_variants[i], self.sample_list
                    ).return_mutation_taster_score()
                )
                + "\t"
                + str(Variant(self.list_variants[i], self.sample_list).parse_counts())
                + "\t"
                + str(
                    len(
                        list(
                            Variant(
                                self.list_variants[i], self.sample_list
                            ).get_genes_in_variant()
                        )
                    )
                )
                + "\t"
                + str(Variant(self.list_variants[i], self.sample_list).counts_cc())
            )
            out.writelines("\n" + results)
        out.close()

    def return_bioinfo_alternative(self, outfile):
        """Returns a tsv file with several information from the variants selected from the
        VCF. Calls variant class to obtain the corresponding fields.
        """

        out = open(outfile, "w")
        out.writelines(
            "Chromosome & Position\tChange at sequence|aa level\tType of change\tGenes\tRs\t"
            + "Polarity change\tNumber of cases affected"
        )

        for i in range(len(self.list_variants)):
            results = ""
            results = (
                Variant(self.list_variants[i], self.sample_list).create_alternative_id()
                + "\t"
                + "".join(Variant(self.list_variants[i], self.sample_list).ids)
                + "\t"
                + Variant(self.list_variants[i], self.sample_list).aminoacid_review()
                + "\t"
                + str(Variant(self.list_variants[i], self.sample_list).counts_cc())
            )
            out.writelines("\n" + results)

        out.close()

    def return_bioinfo_clinvar(self, outfile):
        """Returns a txt file with several information from clinvar
        variants VCF.

        :returns: id from variant, rs, clinic diagnosis, id from clinvar,
        disease related, gene affected in clinvar
        """

        with open(outfile, "w") as out:
            out.writelines(
                "id\tRefSNP\tGene Symbol\tclinical significance\treview\teffects\tscore"
            )
            # out.writelines("id\trs\tclinic\tvariant\tdisease\tgene")

            for i in range(len(self.list_variants)):
                results = ""  # Calls variant class and get the info needed
                # from a clinvar variant

                id = Variant(self.list_variants[i], self.sample_list).vId
                rs = "".join(Variant(self.list_variants[i], self.sample_list).ids)
                variant = Variant(self.list_variants[i], self.sample_list)
                gene = variant.clinvar_d["clinvar_gene"]
                significance = variant.clinvar_d["clinvar_significance"]
                review = variant.clinvar_d["clinvar_review"]
                effects = variant.clinvar_d["effects"]
                results = "{}\t{}\t{}\t{}\t{}\t{}".format(
                    id, rs, gene, significance, review, effects
                )

                out.writelines("\n" + results)
            out.close()

    def return_id_variant_all(self) -> list:
        """Returns a list file with the id of variants (chr,pos,ref,alt)

        Returns:
            list: List with variant ids.
        """

        return [
            Variant(i, self.sample_list).return_id_variant(sep="_")
            for idx, i in enumerate(self.list_variants)
        ]

    def return_infoToburden(self):

        header = list([self.simple_header + "\n"])
        results = []
        for i in range(len(self.list_variants)):
            results.append(
                Variant(self.list_variants[i], self.sample_list).return_variant(True)
            )
        return [*header, *results]

    def __variant_to_rebuild(
        self, chromosome, pos, ids, ref, alt, quality, filtered, info, gt, genotypes
    ):
        """Rebuild a specific variant using it fields. These fields can be
        modified

        ...

        """

        # New fields attributes are assigned
        self.chromosome = chromosome
        self.pos = pos
        self.ids = ids
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.filtered = filtered
        self.info = info
        self.gt = gt
        self.genotypes = genotypes

        return (
            # Returns a line with the new fields of the variant tab delimited
            self.chromosome
            + "\t"
            + self.pos
            + "\t"
            + self.ids
            + "\t"
            + self.ref
            + "\t"
            + self.alt
            + "\t"
            + self.quality
            + "\t"
            + self.filtered
            + "\t"
            + self.info
            + "\t"
            + self.gt
            + "\t"
            + "\t".join(self.genotypes)
        )
