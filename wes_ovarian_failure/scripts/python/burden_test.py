import os
from collections import defaultdict
from re import A

import pandas as pd

from class_vcf import *


def join_opencga_files(dir: str) -> None:
    """_summary_

    Args:
        dir (str): _description_
    """
    opencga_files = [
        f"{opencga_dir}/{f}"
        for f in os.listdir(opencga_dir)
        if os.path.isfile(os.path.join(opencga_dir, f))
    ]

    opencga_d = {}
    for file in opencga_files:
        fhand = open(file, "r")
        next(fhand)
        dhand = {
            "_".join(i): j
            for (i, j) in (
                (line.split("\t")[0:4], (line.split("\t")[4:]))
                for line in fhand.readlines()
            )
        }
        opencga_d.update(dhand)

    with open("/23_burden_data/opencga_variants_combined.tsv", "w") as fout:
        fout.write(
            "chr	pos	ref	alt	consequence	gene	cadd	revel	genotype	clinvar	ref_allele_gnomad	alt_allele_gnomad	alt_freq_gnomad	ref_allele_1000g	alt_allele_1000g	alt_freq_1000g\n"
        )

        for k, v in opencga_d.items():
            values = "\t".join([i for i in v])
            vid = k.replace("_", "\t")
            fout.write(f"{vid}\t{values}")


def generate_variant_id(*strings):
    """_summary_

    Returns:
        _type_: _description_
    """
    res = "_".join([str(i) for i in strings])
    return res


def get_vid_opencga(file: str, generate_vid=False) -> set:
    """_summary_

    Args:
        file (str): _description_

    Returns:
        set: _description_
    """

    df = pd.read_csv(file, sep="\t")

    if generate_variant_id:
        df["vid"] = df.apply(
            lambda row: generate_variant_id(
                row["chr"], row["pos"], row["ref"], row["alt"]
            ),
            axis=1,
        )

    return {i for i in df["vid"].tolist()}


def biomart_to_burden(file_path: str, outpath: str) -> None:
    """_summary_

    Args:
        file_path (str): _description_
        outpath (str): _description_
    """
    with open(file_path, "r") as fhand:
        next(fhand)
        lines = [i.strip() for i in fhand.readlines()]

    with open(outpath, "w") as fout:

        res = []
        for i in lines:
            if len(i.split("\t")) > 3:
                gene, chrom, start, end = i.split("\t")
                res.append(f"{gene}\t{chrom}:{start}-{end}\n")

        fout.writelines(res)


def get_genotype_using_variant_id(file: str, output: str) -> None:
    """_summary_

    Args:
        file (str): _description_
        output (str): _description_
    """
    set_vid = pd.read_csv(file, sep="\t")["vid"].to_list()
    output = open(output, "w")

    # vcf file
    file_vcf = (
        "/home/ihenarejos/workspace/projects/pof/data/merge_vcfs_v38/merged_v38.vcf"
    )

    vfile = Vcf(file_vcf)
    vfile.num_variants()
    header = vfile.header[-1].replace("UNKNOWN_150", "CONTROL_150")

    variants_filtered = vfile.get_variants_using_id(id_list_set=set_vid)
    variants_filtered.num_variants()
    variants_info = variants_filtered.return_infoToburden()
    output.writelines(variants_info)


def get_variants_1000g(dir_1000g: str, genotypes: str) -> None:
    """_summary_

    Args:
        dir_1000g (str): _description_
        genotypes (str): _description_
    """
    dir_1000g = "/home/ihenarejos/workspace/projects/pof/data/1000g_grch38"
    only_files = [
        (f"{dir_1000g}/{f}")
        for f in os.listdir(dir_1000g)
        if re.search(".vcf.gz$", f) and os.path.isfile(f"{dir_1000g}/{f}")
    ]
    with open(
        genotypes,
        "r",
    ) as fhand:
        next(fhand)
        lines = [i.strip() for i in fhand.readlines()]

        chr_dict = defaultdict()
        for i in lines:
            chrom, pos = i.split("\t")[0:2]
            if f"chr{chrom}" not in chr_dict:
                chr_dict[f"chr{chrom}"] = []
                chr_dict[f"chr{chrom}"].append(f"{chrom}\t{pos}\n")

            else:
                chr_dict[f"chr{chrom}"].append(f"{chrom}\t{pos}\n")

        with open(f"{dir_1000g}/bcftools_filter.sh", "w") as fbcf:

            for idx, f in enumerate(only_files):

                chrom = f.split("/")[-1].split(".")[1]

                if chrom == "wgs":

                    continue

                for k, v in chr_dict.items():
                    if k == chrom:
                        with open(f"{dir_1000g}/regions_{chrom}.tsv", "w") as fout:
                            fout.writelines(v)

                if idx == 0:
                    header = "-h"
                    arrow = ">"

                else:
                    header = "-H"
                    arrow = ">>"

                string = f"bcftools view -Ov {header} -S female_samples.tsv ALL.{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -R regions_{chrom}.tsv --force-samples {arrow}aux\n"

                fbcf.writelines(string)

            # at the end of the for loop:
            fbcf.writelines(f"cat aux | tail -n +22 > burden.vcf\nrm aux")


def retrieve_af_population(vid_list: list, output_file: str) -> None:
    """_summary_

    Args:
        vid_list (list): _description_
        gt_idx (list): _description_
    """

    def calculate_alt_freq(genotypes: list) -> float:
        def calculate_genotypes(samples: list) -> list:
            sol = 0
            for i in samples:
                alleles = [int(j) for j in i.split("/") if j != "."]
                if "." not in alleles:
                    sol += sum(alleles)

            return sol

        pop_sum = calculate_genotypes(genotypes)
        if pop_sum != 0:

            af = pop_sum / (len(genotypes) * 2)

            return af

        else:
            return 0

    with open(
        "/home/ihenarejos/workspace/projects/pof/results/24_burden_analysis/v38_opencga_burden_format_changes.tsv",
        "r",
    ) as fhand:

        lines = [i.strip() for i in fhand.readlines()]

    with open(
        output_file,
        "w",
    ) as fout:

        fout.write(f"vid\tcases\tcontrols\n")

        i = -1

        for line in lines:

            variant = line.split("\t")

            if i == -1:
                gt_fields = "\t".join(variant[9:])
                control_idx = [
                    idx
                    for idx, x in enumerate(gt_fields.split("\t"))
                    if x.split("_")[0] in ["CONTROL", "UNKNOWN"]
                ]

                cases_idx = [
                    idx
                    for idx, x in enumerate(gt_fields.split("\t"))
                    if x.split("_")[0] in ["FOO", "FOP"]
                ]

                i = 0
                continue

            else:

                chrom, pos, ref, alt, gt = (
                    variant[0],
                    variant[1],
                    variant[3],
                    variant[4].replace(
                        ",*", ""
                    ),  # alt can have this format alt,* when dealing with indels
                    variant[9:],
                )

                vid = "_".join([chrom, pos, ref, alt])

                if vid in vid_list:
                    maf_cases = calculate_alt_freq(gt[cases_idx[0] : cases_idx[-1]])
                    maf_controls = calculate_alt_freq(
                        gt[control_idx[0] : control_idx[-1]]
                    )

                    fout.write(f"{vid}\t{maf_cases}\t{maf_controls}\n")


if __name__ == "__main__":

    mode = "1000g"

    if mode == "join_variants":

        opencga_dir = "/home/ihenarejos/workspace/projects/uk_stage/data/06_tmp_files/opencga_variants"
        join_opencga_files(opencga_dir)
        exit()

    if mode == "mart":
        biomart_to_burden(
            "data/bed_v38/mart_export.txt",
            "results/24_burden_analysis/biomart_burden.tsv",
        )

    if mode == "pheno":

        sep = "\t"
        header = f"fid{sep}iid{sep}fatid{sep}matid{sep}sex{sep}of\n"

        with open("results/23_burden_data/ivi_samples.txt") as fhand:
            samples = [i.strip() for i in fhand.readlines()]

        with open("results/24_burden_analysis/samples_burden.tsv", "w") as fout:
            res = []
            res.append(header)
            for s in samples:
                if s.split("_")[0] in ["CONTROL", "UNKNOWN"]:
                    sample_class = "2"
                else:
                    sample_class = "1"

                res.append(f"{s}{sep}{s}{sep}{0}{sep}{0}{sep}{2}{sep}{sample_class}\n")

            fout.writelines(res)

    if mode == "maf":

        print(
            retrieve_af_population(
                ["1_942451_T_C", "1_69428_T_G", "1_69511_A_G"],
                output_file="results/24_burden_analysis/v38_variants_burden_format_filtered_changes_opencga_mafs.tsv",
            )
        )
