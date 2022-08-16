# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:07:12 2019

@author: ihenarejos
"""

import re
import statistics
import traceback
import warnings
from collections import OrderedDict, defaultdict
from pprint import *

# from sympy import false


class Variant:
    """_summary_"""

    def __init__(self, variant, sample_list=None):

        #  Fields from variant
        self.variant = variant
        self.chrom = []
        self.pos = 0
        self.ids = []
        self.ref = ""
        self.alt = []
        self.quality = 0
        self.filtered = []
        self.original_info = ""
        self.info = {}
        self.format = []
        self.genotypes = []

        # Specific data
        self.sample_data = {}
        self.clinvar_d = {}
        self.snpEff = {}
        self.sample_list = sample_list
        self.sample_types = {}
        self.sample_types_case_control = {}
        self.sample_outsiders = {}
        self.vId = ""
        self.predictors_scores = {}
        self.counts = ""
        self.csq = {}
        self.csq_alleles = {}

        self.parse()

    def parse(self):

        splits = self.variant.split("\t")
        self.chrom = splits[0]
        self.pos = int(splits[1])
        self.ids = splits[2].split(";")
        self.ref = splits[3]
        self.alt = splits[4].split(",")

        if splits[5] == ".":
            self.quality = "."
        else:
            self.quality = float(splits[5])

        self.filtered = splits[6].split(";")

        if len(splits) > 8:  # format field could be optional
            self.format = splits[8].split(":")
        if len(splits) > 9:  # Genotype field is optional
            self.genotypes = splits[9:]

        self.vId = (
            str(self.chrom)
            + "_"
            + str(self.pos)
            + "_"
            + str(self.ref)
            + "_"
            + str(",".join(self.alt))
        )

        # parse counts
        self.counts = self.parse_counts()

        # parse INFO field
        info_d = {}
        # csq_d = {}
        self.original_info = splits[7]
        info_fields = splits[7].split(";")
        for info_elem in info_fields:
            info_elem_splits = info_elem.split("=")

            if len(info_elem_splits) == 1:
                info_key = info_elem_splits[0]
                info_value = ""
            else:
                info_key = info_elem_splits[0]
                info_value = info_elem_splits[1]
            info_d[info_key] = info_value
        self.info = info_d

        #     if 'CSQ' in info_elem_splits:

        #         # save info without CSQ
        #         indices = [i for i, s in enumerate(info_fields) if 'CSQ=' in s]

        #         for indx in indices:
        #             info_fields.pop(indx)
        #         self.original_info = ';'.join(info_fields)

        #         alleles = info_elem_splits[1].split(',')

        #         for allele in alleles:
        #             al = allele[0]
        #             csq_d[al] = allele

        #         self.csq = csq_d
        # self.parse_CSQ_annotations()

        # parse clinvar annotations in info field
        # self.parse_clinvar_annotations()

        # parse snpeff annotations in INFO field
        self.parse_snpeff_annotations()

        # parse predictors scores in INFO field
        # self.parse_predictors_annotations()

        # parse genotype (sample_data) field
        sample_d = {}

        if len(self.sample_list) == 0:
            warnings.warn(f"proceeding without sample list.")

            for i in range(len(self.genotypes)):
                gt_d = {}
                sample = self.genotypes[i]
                sample_splits = sample.split(":")
                index = i + 1

                if len(sample_splits) == len(self.format):

                    for j in range(len(self.format)):
                        gt_d[self.format[j]] = sample_splits[j]
                        # sample_d[sample_id] = gt_d
                        sample_d = gt_d

            self.sample_data = sample_d

        else:
            if len(self.sample_list) > 0:
                if len(self.sample_list) == len(self.genotypes):

                    for i in range(len(self.sample_list)):
                        gt_d = {}
                        sample = self.genotypes[i]
                        sample_id = self.sample_list[i]
                        sample_splits = sample.split(":")
                        if len(sample_splits) == len(self.format):
                            for j in range(len(self.format)):
                                gt_d[self.format[j]] = sample_splits[j]
                                sample_d[sample_id] = gt_d

                    self.sample_data = sample_d

                    # Case control dictionary
                    for k in range(len(self.sample_list)):
                        class_patient = self.sample_list[k].split("_")[0]
                        if class_patient != "CONTROL":
                            class_patient = "CASE"

                        self.sample_types_case_control[class_patient] = 0

                    for k in range(len(self.sample_list)):
                        class_patient = self.sample_list[k].split("_")[0]
                        if class_patient != "CONTROL":
                            class_patient = "CASE"

                        self.sample_types_case_control[class_patient] = (
                            self.sample_types_case_control[class_patient] + 1
                        )

                    # Group dictionary

                    for k in range(len(self.sample_list)):
                        class_patient = self.sample_list[k].split("_")[0]
                        self.sample_types[class_patient] = 0

                    for k in range(len(self.sample_list)):
                        class_patient = self.sample_list[k].split("_")[0]

                        if class_patient == "UNKNOWN":
                            class_patient = "CONTROL"

                        self.sample_types[class_patient] = (
                            self.sample_types[class_patient] + 1
                        )

                    # Different populations dictionary
                    list_not_north = [
                        "FOP_28",
                        "FOP_29",
                        "FOP_30",
                        "FOP_31",
                        "FOP_32",
                        "FOP_33",
                        "FOP_34",
                        "FOP_35",
                    ]

                    for k in range(len(self.sample_list)):
                        if self.sample_list[k] in list_not_north:
                            class_patient = "EAST"
                        else:
                            class_patient = "NORTH"

                        self.sample_outsiders[class_patient] = 0

                    for k in range(len(self.sample_list)):
                        if self.sample_list[k] in list_not_north:
                            class_patient = "EAST"
                        else:
                            class_patient = "NORTH"

                        self.sample_outsiders[class_patient] = (
                            self.sample_outsiders[class_patient] + 1
                        )

    def analyze_missing(self, missing_threshold=0.5):
        """Analyze if variants has too many missing values.

        Parameters
        ----------
        missing_threshold : float
        ...

        """

        missing = 0

        for i in range(len(self.genotypes)):
            genotype_sample = self.genotypes[i]
            haplotype = genotype_sample.split(":")[0]
            if haplotype == "./.":
                missing += 1

        #  Compare against threshold
        if missing / (len(self.genotypes)) > missing_threshold:
            return False
        else:
            return True

    def get_snpsift_case_control(self):

        snpsift_annotations = [self.vId]

        if "CC_DOM" in self.info:
            snpsift_annotations.append(self.info["CC_DOM"])
        else:
            snpsift_annotations.append("")

        if "CC_REC" in self.info:
            snpsift_annotations.append(self.info["CC_REC"])
        else:
            snpsift_annotations.append("")

        if "CC_ALL" in self.info:
            snpsift_annotations.append(self.info["CC_ALL"])
        else:
            snpsift_annotations.append("")

        if "CC_GENO" in self.info:
            snpsift_annotations.append(self.info["CC_GENO"])
        else:
            snpsift_annotations.append("")

        if "CC_TREND" in self.info:
            snpsift_annotations.append(self.info["CC_TREND"])
        else:
            snpsift_annotations.append("")

        return snpsift_annotations

    def parse_clinvar_annotations(self):

        clinvar_d = {}

        if "CLNDN" in self.info:
            clinvar_d["clinical_tr"] = self.info["CLNDN"]
        else:
            clinvar_d["clinical_tr"] = "."

        if "CLNHGVS" in self.info:
            clinvar_d["clinical_v"] = self.info["CLNHGVS"]
        else:
            clinvar_d["clinical_v"] = "."

        if "MC" in self.info:
            clinvar_d["clinvar_eff"] = self.info["MC"]
        else:
            clinvar_d["clinvar_eff"] = "."

        if "GENEINFO" in self.info:
            clinvar_d["clinvar_gene"] = self.info["GENEINFO"].split(":")[0]
            clinvar_d["clinvar_full"] = self.info["GENEINFO"]
        else:
            clinvar_d["clinvar_gene"] = "."
            clinvar_d["clinvar_full"] = "."

        if "CLNREVSTAT" in self.info:
            clinvar_d["clinvar_review"] = self.info["CLNREVSTAT"]
        else:
            clinvar_d["clinvar_review"] = "."

        if "CLNSIG" in self.info:
            clinvar_d["clinvar_significance"] = self.info["CLNSIG"]
        else:
            clinvar_d["clinvar_significance"] = "."

        # Variant effects
        if "MC" in self.info:
            effects = self.info["MC"].split("|")[1:]

            # if there are more than one effect:
            if len(effects) > 1:
                res = []
                for j in effects:
                    head, sep, tail = j.partition(",")
                    res.append(head)
                effects = ";".join(res)
                clinvar_d["effects"] = effects
            else:
                effects = ";".join(effects)
                clinvar_d["effects"] = effects
        else:
            clinvar_d["effects"] = "."

        self.clinvar_d = clinvar_d
        return self.clinvar_d

    def parse_snpeff_annotations(self):
        """Returns snpEff annotations of variant in a dict form where each
        key is the unique allele and for each allele the value is a list of
        dicts for each value affected and where the keys are the snpEff
        fields with their respective values.

        Returns
        -------
        dict
            snpEff annotations
        ...

        """

        snpeff_dict = {}

        list_snpeff = (
            "Allele",
            "Annotation",
            "Annotation_Impact",
            "Gene_Name",
            "Gene_ID",
            "Feature_Type",
            "Feature_ID",
            "Transcript_BioType",
            "Rank",
            "HGVS.c",
            "HGVS.p",
            "cDNA.pos/cDNA.length",
            "CDS.pos/CDS.length",
            "AA.pos/AA.length",
            "Distance",
            "ERRORS/WARNINGS/INFO",
        )

        if "ANN" in self.info:
            annotations = self.info["ANN"]
            transcripts = annotations.split(",")
            for transcript in transcripts:

                d = {}
                transcripts_split = transcript.split("|")
                t_key = transcripts_split[0]

                if t_key not in snpeff_dict:
                    snpeff_dict[t_key] = []

                if len(transcripts_split) == len(list_snpeff):
                    for annotation, split in zip(list_snpeff, transcripts_split):
                        d[annotation] = split
                snpeff_dict[t_key].append(d)

        self.snpEff = snpeff_dict
        return self.snpEff

    def parse_CSQ_annotations(self):

        csq_list = [
            "Allele",
            "Consequence",
            "IMPACT",
            "SYMBOL",
            "Gene",
            "Feature_type",
            "Feature",
            "BIOTYPE",
            "EXON",
            "INTRON",
            "HGVSc",
            "HGVSp",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "Existing_variation",
            "DISTANCE",
            "STRAND",
            "FLAGS",
            "SYMBOL_SOURCE",
            "HGNC_ID",
            "SOURCE",
            "CADD_PHRED",
            "CADD_RAW",
            "gnomADe",
            "gnomADe_AF",
        ]

        csq_alleles = defaultdict(dict)

        for k in self.csq.keys():
            split = self.csq[k].split("|")

            if len(split) == len(csq_list):
                for i in range(0, len(csq_list)):
                    format = csq_list[i]
                    csq_alleles[k][format] = split[i]

        self.csq_alleles = csq_alleles

        return self.csq_alleles

    def get_dp(self):
        """
        Returns DP value for the variant

        """

        dp_val = 0

        if "DP" in self.info:
            dp_annotation = self.info["DP"]
            dp_val = dp_annotation.split("=")[0]

        return dp_val

    def get_mean_dp_samples(self):
        """
        based on DP of all variants selected
        mean is computed for all samples where variant was sequenced
        meaning, missing values are not considered, since it could mean that
        sample wasn't sequenced in that position
        """
        presence = 0
        results = 0
        for key in self.sample_data:
            dict_sample = self.sample_data[key]
            for k in dict_sample:
                if k == "DP":
                    value = dict_sample[k]
                    print(value)
                    if value != ".":
                        # value = value.split(",")[1]
                        presence = presence + 1
                        results = results + int(value)

        if presence == 0:

            return 0

        else:

            return results / presence

    def get_sd_dp_samples(self):
        """
        based on DP of all variants selected
        mean is computed for all samples where variant was sequenced
        meaning, missing values are not considered, since it could mean that
        sample wasn't sequenced in that position
        """
        presence = 0
        results = []
        for key in self.sample_data:
            dict_sample = self.sample_data[key]
            for k in dict_sample:
                if k == "DP":  # AD, DP
                    value = dict_sample[k]
                    if value != ".":
                        # value = value.split(",")[1]
                        presence = presence + 1
                        results.append(int(value))

        if len(results) >= 2:
            return round(statistics.stdev(results), 2)

        else:
            return 0

    def set_samples_affected(self):
        """Check samples affected in each variant
        :return: set of samples that are affected
        """

        results = set()
        for sample in self.sample_data:
            sample_dict = self.sample_data[sample]

            for key in sample_dict:
                if key == "GT":
                    haplotype = sample_dict[key]
                    alleles = haplotype.split("/")

                    for allele in alleles:
                        if allele != "." and int(allele) >= 1:
                            results.add(str(sample))

        return results

    def sample_format_field(self, sample_id, format_fields=[]):
        """
        sample_id: str
        format_fields: list

        Returns either a string or a dictionaries with the information from the sample
        field and format specified

        Returns
        -------
        str if only one field specified
        dict if more than one field specified

        """
        if type(format_fields) != type(list()):
            try:
                format_fields = [format_fields]
            except:
                print(f"Wrongly formatted format field. Exitting...")
                traceback.print_exc()
                exit()

        if len(format_fields) == 1:
            spec = format_fields[0]
            return self.sample_data[sample_id][spec]

        if len(format_fields) > 1:
            f_dict = {k: self.sample_data[sample_id][k] for k in format_fields}
            return f_dict

    def check_quality(self, sample_id):
        """Check if a sample specified passes the quality thresholds defined
        in study

        :param sample_id: str: sample id
        :return: int : add 1 to count
        """

        result = 0

        dict_qualities = self.sample_data[sample_id]
        gt = int(dict_qualities["GT"])
        ad = int(dict_qualities["AD"])
        dp = int(dict_qualities["DP"])

        if gt > 90 and ad > 4 and dp > 10 and ad / dp > 0.2:
            result = 1
        else:
            result = 0

        return result

    def check_haplotype_all_samples(self):
        """Creates a dict where each sample has its genotype value stored (
        ex CONTROL_118: REF) missing values are treated depending on the
        interests of the study
        """

        dict_samples = self.sample_data

        index = 0

        samples_haplotype = {}

        for i in range(len(self.genotypes)):

            if len(self.sample_list) > 0:
                index = self.sample_list[i]

            elif len(self.sample_list) == 0:
                index += 1
            haplotype = dict_samples[index]["GT"]

            if self.__check_if_homozygous_for_mut(haplotype):
                samples_haplotype[index] = "hom"

            elif self.__check_if_heterozygous_for_mut(haplotype):
                samples_haplotype[index] = "het"

            elif haplotype == "0/0":
                samples_haplotype[index] = "ref"

            elif self.__check_if_missing_for_mut(haplotype):  # change to
                # desired value
                pass
                # samples_haplotype[index]="ref" # OR "miss"

        return samples_haplotype

    @staticmethod
    def __check_if_reference(haplotype):
        """Returns boolean depending what conditions are meet"""

        if haplotype == "0/0":
            return True
        else:
            return False

    @staticmethod
    def __check_if_missing_for_mut(haplotype):
        """Returns boolean depending what conditions are meet"""

        if haplotype == "./.":

            return True
        else:
            return False

    @staticmethod
    def __check_if_homozygous_for_mut(haplotype):
        """Returns boolean depending what conditions are meet"""

        if haplotype != "0/0" and haplotype != "./.":

            split = haplotype.split("/")  # / / | \|

            allele_1 = split[0]

            allele_2 = split[1]

            if allele_1 == "0" or allele_2 == "0":
                return False
            else:
                return True

    @staticmethod
    def __check_if_heterozygous_for_mut(haplotype):
        """Returns boolean depending what conditions are meet"""
        split = haplotype.split("/")  # / / | \|
        allele_1 = split[0]
        allele_2 = split[1]

        if (allele_1 == "0" and allele_2 != "0") or (
            allele_2 == "0" and allele_1 != "0"
        ):
            return True
        else:
            return False

    def genotype_samples_variant(self):
        """
        Returns a dict where keys are the samples and they are given a
        certain value depending if alt is present.
        Missing values are treated as ref.
        """

        dict_values = {}

        for i in range(len(self.genotypes)):
            if len(self.sample_list) > 0:
                index = self.sample_list[i]
            haplotype = self.sample_data[index]["GT"]

            if self.__check_if_homozygous_for_mut(haplotype):
                dict_values[index] = "hom"
            elif self.__check_if_heterozygous_for_mut(haplotype):
                dict_values[index] = "het"
            elif self.__check_if_reference(haplotype):
                dict_values[index] = "ref"
            elif self.__check_if_missing_for_mut(haplotype):
                dict_values[index] = "ref"
        # Save and force order to keys in dictionary
        dict_values = OrderedDict(sorted(dict_values.items()))

        return dict_values

    def presence_absence(self, tab_delimited=False):
        """

        :return: str with 150 characters
        :att : tab_delimited T/F
        """

        res = ""

        for sample in self.sample_data:
            sample_dict = self.sample_data[sample]

            for key in sample_dict:
                if key == "GT":
                    haplotype = sample_dict[key]
                    alleles = haplotype.split("/")

                    for allele in alleles:
                        if allele != "." and int(allele) >= 1:

                            if not tab_delimited:
                                res = res + "1"
                            else:
                                res = res + "1" + "\t"

                        else:
                            if not tab_delimited:
                                res = res + "0"
                            else:
                                res = res + "0" + "\t"

        return res

    def calculate_proportions(self):
        """Using check_all_samples function to get which samples are hom/het
        or ref, calculate proportions
        for that variant
        """

        hom_c = 0
        het_c = 0
        ref_c = 0
        d_counts = {}
        samples_haplotype = self.check_haplotype_all_samples()

        for key in self.sample_types:  #
            """Creates a dict of the form CONTROL { HET: 0, HOM:0, REF: O,
            MISS:0} Add miss key at will.

            """
            d_counts[key] = {}
            d_counts[key]["hom"] = 0
            d_counts[key]["ref"] = 0
            d_counts[key]["het"] = 0

        for key in samples_haplotype:
            gt = samples_haplotype[key]
            key_type = key.split("_")[0]
            d_counts[key_type][gt] += 1

        for key in d_counts:
            d_counts[key]["hom"] = round(
                (d_counts[key]["hom"] / self.sample_types_case_control[key]) * 100
            )
            d_counts[key]["het"] = round(
                (d_counts[key]["het"] / self.sample_types_case_control[key]) * 100
            )
            d_counts[key]["ref"] = round(
                (d_counts[key]["ref"] / self.sample_types_case_control[key]) * 100
            )

        return (
            str(d_counts["CONTROL"]["ref"])
            + "\t"
            + str(d_counts["CONTROL"]["het"])
            + "\t"
            + str(d_counts["CONTROL"]["hom"])
            + "\t"
            + str(d_counts["FOO"]["ref"])
            + "\t"
            + str(d_counts["FOO"]["het"])
            + "\t"
            + str(d_counts["FOO"]["hom"])
            + "\t"
            + str(d_counts["FOP"]["ref"])
            + "\t"
            + str(d_counts["FOP"]["het"])
            + "\t"
            + str(d_counts["FOP"]["hom"])
        )

    def calculate_counts_cc(self):
        """Using check_all_samples function to get which samples are hom/het or
        ref, calculate counts for case vs control
        """

        hom_c = 0
        het_c = 0
        ref_c = 0
        d_counts = {}
        samples_haplotype = self.check_haplotype_all_samples()

        for key in self.sample_types_case_control:
            """Create a dict of the form CONTROL { HET: 0, HOM:0, REF: O,
            MISS:0}

            """
            d_counts[key] = {}
            d_counts[key]["hom"] = 0
            d_counts[key]["ref"] = 0
            d_counts[key]["het"] = 0

        for key in samples_haplotype:

            key_type = key.split("_")[0]
            if key_type == "FOO" or key_type == "FOP":
                key_type = "CASE"

            gt = samples_haplotype[key]
            d_counts[key_type][gt] += 1

        return (
            str(d_counts["CONTROL"]["ref"])
            + "\t"
            + str(d_counts["CONTROL"]["het"])
            + "\t"
            + str(d_counts["CONTROL"]["hom"])
            + "\t"
            + str(d_counts["CASE"]["ref"])
            + "\t"
            + str(d_counts["CASE"]["het"])
            + "\t"
            + str(d_counts["CASE"]["hom"])
        )

    def calculate_counts_haplotypes(self):
        """Using check_all_samples function to get which samples are hom/het or
        ref, calculate counts for case vs control
        """

        hom_c = 0
        het_c = 0
        ref_c = 0
        d_counts = {}
        samples_haplotype = self.check_haplotype_all_samples()

        for key in self.sample_types_case_control:
            """Create a dict of the form CONTROL { HET: 0, HOM:0, REF: O,
            MISS:0}

            """
            d_counts[key] = {}
            d_counts[key]["hom"] = 0
            d_counts[key]["ref"] = 0
            d_counts[key]["het"] = 0

        for key in samples_haplotype:

            key_type = key.split("_")[0]
            if key_type == "FOO" or key_type == "FOP":
                key_type = "CASE"

            gt = samples_haplotype[key]
            d_counts[key_type][gt] += 1

        h00 = int(d_counts["CONTROL"]["ref"]) + int(d_counts["CASE"]["ref"])
        h01 = int(d_counts["CONTROL"]["het"]) + int(d_counts["CASE"]["het"])
        h11 = int(d_counts["CONTROL"]["hom"]) + int(d_counts["CASE"]["hom"])

        return str(h00) + "\t" + str(h01) + "\t" + str(h11)

    def calculate_counts_diff_populations(self):
        """Count and differs from population specified"""

        hom_c = 0
        het_c = 0
        ref_c = 0
        d_counts = {}
        samples_haplotype = self.check_haplotype_all_samples()

        list_not_north = [
            "FOP_28",
            "FOP_29",
            "FOP_30",
            "FOP_31",
            "FOP_32",
            "FOP_33",
            "FOP_34",
            "FOP_35",
        ]

        for key in self.sample_outsiders:
            d_counts[key] = {}
            d_counts[key]["hom"] = 0
            d_counts[key]["ref"] = 0
            d_counts[key]["het"] = 0

        for key in samples_haplotype:
            key_type = ""
            if key in list_not_north:
                key_type = "EAST"
                gt = samples_haplotype[key]
                d_counts[key_type][gt] += 1
            else:
                key_type = "NORTH"
                gt = samples_haplotype[key]
                d_counts[key_type][gt] += 1

        return (
            str(d_counts["Y"]["ref"])
            + "\t"
            + str(d_counts["Y"]["het"])
            + "\t"
            + str(d_counts["Y"]["hom"])
            + "\t"
            + str(d_counts["X"]["ref"])
            + "\t"
            + str(d_counts["X"]["het"])
            + "\t"
            + str(d_counts["X"]["hom"])
        )

    def get_alt_alleles_changes(self):
        """Get all alts alleles in variant

        Returns
        ------
        dict :
            dict where each key is an alt allele and values are sets

        ...

        """

        alleles = {}

        alleles_values = self.alt

        for allele in alleles_values:
            alleles[allele] = set()

        for allele in self.snpEff:

            if allele in alleles:
                alleles[allele] = self.snpEff[allele]

        return alleles

    def which_sample_format(self, cond):
        """Returns which sample meets conditions in its format data"""

        results = []

        for key in self.sample_data:
            dic = self.sample_data[key]
            for k in dic:
                for element in cond:
                    if dic[k] == element:
                        results.append(key)

        return "Conditions meet %s: ", results

    def get_predictors_annotations(self):
        """Returns a string of predictors scores using dictionary
        comprehension

        """

        self.parse_predictors_annotations()

        return {k: v for (k, v) in self.predictors_scores.items()}

    def return_cadd_score(self):
        """Returns a string of CADD score"""
        self.parse_predictors_annotations()

        if "dbNSFP_CADD_phred" in self.predictors_scores:
            value = self.predictors_scores["dbNSFP_CADD_phred"]
        else:
            value = ""

        return str(value)

    def return_mutation_taster_score(self):
        """Returns a string of CADD score"""

        self.parse_predictors_annotations()

        if "dbNSFP_MutationTaster_pred" in self.predictors_scores:
            value = self.predictors_scores["dbNSFP_MutationTaster_pred"]
        else:
            value = ""

        return str(value)

    def parse_predictors_annotations(self):
        """Returns a dict of predictors scores"""

        predictors_d = {}

        list_predictors = [
            "dbNSFP_CADD_phred",
            "dbNSFP_MutationTaster_pred",
            "dbNSFP_Polyphen2_HDIV_pred",
            "dbNSFP_Polyphen2_HVAR_pred",
            "dbNSFP_SIFT_pred",
        ]

        for element in list_predictors:

            if element in self.info:
                values = self.info[element]
                predictors_d[element] = values

        self.predictors_scores = predictors_d
        return self.predictors_scores

    def parse_counts(self):
        """Returns counts as att"""
        counts_d = {}
        list_of_counts = ["control_f", "foo_f", "fop_f"]
        res = []

        for element in list_of_counts:
            if element in self.info:
                values = self.info[element]
                counts_d[element] = values
                res.append(values)
                # print(element, res)

        self.counts = counts_d

        # return self.counts
        return "\t".join(map(str, res))

    def counts_cc(self):
        """

        :return: str case str control
        """

        self.parse_counts()
        foo = 0
        fop = 0
        control = 0

        for key in self.counts:

            if key == "foo_f":
                foo = self.counts["foo_f"]

            elif key == "fop_f":
                fop = self.counts["fop_f"]

            elif key == "control_f":
                control = self.counts["control_f"]

        cases = int(foo) + int(fop)

        return str(cases)
        # return str(cases) + "\t" + str(control)

    def get_genes_in_variant(self):
        """
        Returns a list of unique genes in variant if transcript is protein
        coding
        """

        genes = list()

        if len(self.snpEff) > 0:

            for allele in self.snpEff:
                list_transcript = self.snpEff[allele]

                for transcript_annotated in list_transcript:
                    gene = transcript_annotated["Gene_Name"]
                    bio_type = transcript_annotated["Transcript_BioType"]

                    if bio_type == "protein_coding":
                        if gene not in genes:
                            genes.append(gene)

        else:
            genes = set()
            for allele in self.csq_alleles.keys():
                genes.add(self.csq_alleles[allele]["SYMBOL"])

        return genes

    def get_genes_in_variant_clinvar(self):
        """
        Returns a list of unique genes in variant if transcript is protein
        coding
        """

        genes = list()

        gene = self.clinvar_d["clinvar_gene"]
        gene = gene.split(":")[0]
        if gene not in genes:
            genes.append(gene)

        return genes

    def set_genes_in_variant_changelist(self, effect_list=None):
        """
        Returns a set of unique genes in variant if gene is protein coding
        and its affected by a change in provided list
        Also returns a dict of alleles and the genes they affect with their
        deleterious effect
        """

        genes = set()
        alleles = {}

        for allele in self.snpEff:
            list_transcript = self.snpEff[allele]

            if allele not in alleles:
                alleles[allele] = {}
            allele_dict = alleles[allele]  # only dict of allele of iteration

            for transcript_annotated in list_transcript:
                gene = transcript_annotated["Gene_Name"]
                bio_type = transcript_annotated["Transcript_BioType"]

                if bio_type == "protein_coding":
                    change_splits = transcript_annotated["Annotation"].split("&")

                    for change in change_splits:
                        if change in effect_list:
                            genes.add(gene)
                            if gene not in allele_dict:
                                allele_dict[gene] = set()
                            # if already in:
                            allele_dict[gene].add(change)

        return genes, alleles

    def get_snpeff_field(self, field):
        """
        Returns a string of the requested snpEff field
        'Allele','Annotation','Annotation_Impact','Gene_Name','Gene_ID',
        'Feature_Type','Feature_ID','Transcript_BioType',\
             'Rank','HGVS.c','HGVS.p','cDNA.pos/cDNA.length',
             'CDS.pos/CDS.length','AA.pos/AA.length','Distance',
             'ERRORS/WARNINGS/INFO')
        """

        value_list = []

        for allele in self.snpEff:
            list_transcripts = self.snpEff[allele]

            for tr in list_transcripts:
                value = tr[field]

                if value not in value_list:
                    value_list.append(value)

        return ",".join(value_list)

    def aminoacid_review(self):
        """
        Creates a dict of aminoacid properties from a file, reads a variant
        and annotates the change of said properties
        :return: str
        """

        in_aa = [
            line.rstrip("\n")
            for line in open(
                "/home/ihenarejos/workspace/projects/pof/data"
                "/txt_tsv"
                "/aminoacids_info.txt",
                "r",
            )
        ]
        aa_dict = {}

        for line in in_aa:
            aa_properties = line.split("\t")
            aa_dict[(aa_properties[1]).replace(" ", "")] = str(
                aa_properties[3] + " " + aa_properties[4] + " " + aa_properties[5]
            )

        allele = self.alt[0]
        first_transcript = {}
        first_transcript = self.snpEff[allele][0]

        aa = (first_transcript["HGVS.p"]).replace("p.", "")
        if len(aa) > 0:

            aa2 = re.sub("\d+", "_", aa)
            # pprint(aa2)

            aa_original = aa2.split("_")[0]
            aa_alternative = aa2.split("_")[1]

            p1 = ""
            p2 = ""

            if aa_original in aa_dict:
                p1 = aa_dict[aa_original]
            if aa_alternative in aa_dict:
                p2 = aa_dict[aa_alternative]

            res = str(p1) + " > " + str(p2)

            return res
        else:
            return ""

    def create_alternative_id(self):

        allele = self.alt[0]
        first_transcript = {}
        first_transcript = self.snpEff[allele][0]
        # print(first_transcript)
        gene = first_transcript["Gene_Name"]
        effect = first_transcript["Annotation"]
        dna = first_transcript["HGVS.c"]
        aa = first_transcript["HGVS.p"]

        cro_pos = (self.chrom).replace("chr", "") + ";" + str(self.pos)
        dna_aa = dna + ", " + aa

        return cro_pos + "\t" + dna_aa + "\t" + effect + "\t" + gene

    def allele_correspondence(self):
        """ """
        corr = {}

        allele_nuc = [self.ref] + [i for i in self.alt]

        for i in allele_nuc:
            idx = str(allele_nuc.index(i))
            corr[idx] = i

        return corr

    def get_snpeff_allele_genes(self):
        """Returns a dict where keys are alleles of variant and values are
        genes affected by that particular allele (based on SnpEff)
        """

        allele_gene = {}

        for allele in self.snpEff:
            list_transcripts = self.snpEff[allele]

            for tr in list_transcripts:
                gene = tr["Gene_Name"]

                if allele not in allele_gene:
                    allele_gene[allele] = set()
                    allele_gene[allele].add(gene)
                if allele in allele_gene:
                    allele_gene[allele].add(gene)

        return allele_gene

    def get_clinvar_in_variant(self, field):
        """Get values from clinvar dictionary"""
        value_list = []

        value = self.clinvar_d[field]

        if value not in value_list:
            value_list.append(value)

        return ",".join(value_list)

    def return_rs(self):

        rs = ",".join(self.ids)

        if rs == 0:
            rs = "."

        return rs

    def return_chr(self):

        return self.chrom

    def return_pos(self):

        return str(self.pos)

    def return_ref(self):

        return self.ref

    def return_alt(self, ignore_indels=False):

        if ignore_indels:

            return (",".join(self.alt)).replace(",*", "")

        else:
            if type(self.alt) is list:

                return ",".join(self.alt)

            else:

                return self.alt

    def return_genotype(self):

        return "\t".join(self.genotypes)

    def return_filter(self):

        return "".join(self.filtered)

    def return_format(self):

        return ":".join(self.format)

    def return_id_variant(self, sep):
        chr, pos = tuple([str(i) for i in [self.chrom, self.pos]])
        return sep.join(
            [chr.replace("chr", ""), pos, self.ref, self.return_alt(ignore_indels=True)]
        )

    def return_variant(self, reduce_info: bool = False) -> list:
        """Returns variant as a string

        Args:
            reduce_info ([bool]): If true, ID, QUAL, FILTER AND INFO fields are left as '.'. Format and genotype fields only return 'GT' string

        Returns:
            [list]: [a list with variants as str]
        """
        # if 'alt' in dict_changes.keys():
        #     return_alt = dict_changes['alt']
        # else:
        return_alt = self.return_alt()

        res = (
            f"{self.chrom}\t{self.pos}\t{self.return_rs()}\t{self.ref}"
            f"\t{return_alt}"
            f"\t{self.quality}\t{self.return_filter()}\t{self.info}"
            f"\t{self.return_format()}\t{self.return_genotype()}"
        )

        if reduce_info:

            format_d = {
                v: indx for indx, v in enumerate(self.return_format().split(":"))
            }
            genotypes_gt = "\t".join(
                [
                    i.split(":")[format_d["GT"]]
                    for i in self.return_genotype().split("\t")
                ]
            )
            # for i in self.return_genotype.split("\t"):
            # gt = i.split(":")[format_d['GT']]

            chrom = self.chrom.replace("chr", "")

            res = (
                f"{chrom}\t{self.pos}\t.\t{self.ref}"
                f"\t{return_alt}"
                f"\t.\t.\t."
                f"\tGT\t{(genotypes_gt)}\n"
            )
        return res

    def info_rebuild_csq(self, allele):

        info = []

        for catg in self.csq_alleles[allele].keys():
            info.append(self.csq_alleles[allele][catg])
        csq = ";CSQ=" + "|".join(info)

        # retrieve original info without csq and add it to rebuilt csq annotations:
        res = self.original_info + csq

        return res

    def split_multiallelic_variant(self):

        not_indel = False
        res = []
        dict_changes = {}

        def set_alt_and_transform(allele):
            # set current alelle as the only allele:
            dict_changes["alt"] = allele

            # get new info field
            new_info = self.info_rebuild_csq(allele)
            self.info = new_info
            geno = []

            # modify genotypes
            if len(self.genotypes) > 1:
                for j in self.genotypes:
                    genotype_split = "".join(j).split(":")
                    genotype_split[0] = re.sub(
                        pattern="[2-9]", repl="1", string=(genotype_split)[0]
                    )
                    genotype_split = ":".join(genotype_split)
                    geno.append(genotype_split)
                self.genotypes = geno

            else:
                genotype_split = "".join(self.genotypes).split(":")
                genotype_split[0] = re.sub(
                    pattern="[2-9]", repl="1", string=(genotype_split)[0]
                )
                genotype_split = ":".join(genotype_split)
                geno.append(genotype_split)
                self.genotypes = geno

            new_variant = self.return_variant(dict_changes)
            res.append(new_variant)

            return res

        # check if alleles are in info
        for k in self.csq_alleles.keys():
            if k not in self.alt:
                res = set_alt_and_transform(k)
            else:
                not_indel = True

        if not_indel == True:
            for i in self.alt:
                res = set_alt_and_transform(i)

        return res


class VariantSample(Variant):
    sample: str
    haplotype: str
    variant: Variant

    def __init__(self, variant, sample) -> None:
        self.variant = variant
        self.sample = sample
        self.haplotype = ""

        self.__get_haplotypes()

    def to_string(self):
        return f"{self.variant}\t{self.sample}\t{self.haplotype}"

    def __get_haplotypes(self) -> str:
        """
        Return a str with the haplotype associated with the variant
        """

        self.haplotype = self.variant.sample_format_field(self.sample, "GT").split(
            "/|\\|"
        )
        self.haplotype = "".join(self.haplotype)
        return self.haplotype

    # def get_changes_specific_gene(self, gene):
    #     '''
    #     Return changes (e.g. missense, frameshift) for the gene, sample provided
    #     '''
    #     changes: dict

    #     changes = {}
    #     alleles_and_genes = self.variant.get_snpeff_allele_genes()
    #     annotations = self.variant.snpEff

    #     # extract the haplotype of a given variant for a certain sample,
    #     # then retrieve its alleles and check if they are reference or not
    #     for allele in self.haplotype:
    #         if allele == '0' or '.':
    #                 gene_score = 1.0 # not really necessary
    #                 # don't include anything in changes, it will count as 1
    #         elif allele in alleles_and_genes.keys():
    #             # iterate over transcripts corresponding to alleles
    #             # and retrieve change
    #             for gene_snpeff in alleles_and_genes[allele]:
    #                 if gene_snpeff == gene:
    #                     changes = {j["Annotation"].split("&") for i in annotations[allele] for j in i}

    #     return changes
