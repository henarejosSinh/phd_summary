# JAN 5 2023
# ihc.ca94@gmail.com
# circa gtf - vcf skill test

import re
import gzip
import argparse
from os import listdir
from collections import defaultdict
from bisect import bisect_right, bisect_left
from multiprocessing import Process

# For this test, we configure a subclass of multiprocessing.Process class to paralelize the vcf files so each core will handle several vcf files at the same time: 

class annotation_process(Process):
    
    # The object will take the full list of vcf files, the corresponding batch for the core, the number of base pairs to add to the variant position up/downstream and a dictionary in the form of (start position, end position): gene id.
    def __init__(self: list, batch: list, gtf_dict: dict, window_limit: int):
    # def __init__(self, vcf_files: list, batch: list, gtf_dict: dict, window_limit: int):
        Process.__init__(self)
        # self.vcf_files = vcf_files
        self.batch = batch
        self.gtf_dict = gtf_dict
        self.window_limit = window_limit
            
    def run(self):
        
        # the core will select its corresponding vcf according to the batch and launch the function to parse the variants in the file, and create a new compressed vcf with the corresponding annotations
        # for i in self.batch:
            # print(self.batch, self.vcf_files[i])
        for i in self.batch:
            print(i)
            
            output_name = i.replace('.gz', '_ann.gz')
            # output_name = self.vcf_files[i].replace('.gz', '_ann.gz')
            
            file_out = gzip.open(output_name, 'wt')
            file_out.writelines(read_vcf_and_annotate(i, self.window_limit, self.gtf_dict))
            # file_out.writelines(read_vcf_and_annotate(self.vcf_files[i], self.window_limit, self.gtf_dict))
            
            file_out.close()
            

def read_vcf_and_annotate(vcf_file: str, limit: int, gtf_dict: dict) -> str:
    """reads a vcf file and substitute info field with the requested annotations.

    Args:
        vcf_file (str): vcf file to open and parse
        limit (int): up/downstream base pairs to add to the variant to search for genes
        gtf_dict (dict): gtf dictionary (start, end): gene id

    Yields:
        _str_: variant with info field annotated
    """ 
    
    # here, we extract the corresponding chromosome from the name of the file
    chr_file = "".join(re.findall(r'chr\d+', vcf_file))
    
    # as we will use the bisect module afterwards, we need to sort the gtf positions according to the start position.
    sorted_positions = sorted(gtf_dict[chr_file].keys(), key=lambda tup: tup[0])
    
    # here, we iterate over every line in the vcf file. If the line is part of the header/metadata, is yielded to the new vcf without further modification. if the line is a variant, is processed to generate a modified string with the new fields in the info section. 
    # using yield instead of return allows us to avoid any issues with RAM as the new lines are not stored in memory but written directly to the new file.
    with gzip.open(vcf_file,'rt') as f:
            
            for line in f:    
                if line.startswith('##'):
                    yield line
                    
                elif line.startswith("#CHROM"):
                    header_ref = {field.strip("#").strip(): i for i, field in enumerate(line.split("\t"))}
                    yield line
                    
                elif len(line) == 0:
                    continue
                
                else:
                    res_variant = parse_variant_and_annotate(line, header_ref, gtf_dict, sorted_positions, limit )
                    
                    if res_variant is not None:
                        yield res_variant


def parse_variant_and_annotate(variant: str, header: dict, dict_gtf: dict, sorted_positions: list, limit: int) -> str:
    """ parses a variant (string) to analyze if it falls inside a gene region, close to it or calculate the nearest gene to the position, using information extracted from a gtf as reference.

    Args:
        variant (str): variant in string format
        header (dict): the format of the header, that will be used to extract target fields in the variant by position
        dict_gtf (dict): dictionary in the form of (start, end): gene id
        sorted_positions (list): sorted genomic positions for each gene in the gtf field, considering the start position
        limit (int): the up/downstream base pairs to add to the variant to study a wider window.

    Returns:
        str: the variant with the info field annotated.
    """    
    
    # we initialize the target fields
    genes_in = set()
    genes_window = set()
    nearest = {}
    
    # variant is separated into fields, and target position and chromosome is extracted
    fields = variant.rstrip().split('\t')
    pos = int(fields[header['POS']])
    gtf_d_chr = dict_gtf[fields[header['CHROM']]]
    
    # now, we perform a first filter, were we omit analyzing variants too far away from the known starting and end points from the gtf file
    # if pos + limit < int(sorted_positions[0][0]) - limit or pos - limit > int(sorted_positions[-1][1]) + limit:
        # return None
    
    # know we use the bisect module to get the closest point to the known starting points of the genes in the gtf file. this will allows us to analyze only the gene positions closest to the target variant position, skipping the rest of the gene positions too far away to take into account. This step requires of python 3.10 to work:
    
    bisect_idx_minus = bisect_left(sorted_positions, abs(pos - limit), key=lambda i: i[0])
    bisect_idx_plus = bisect_right(sorted_positions, abs(pos + limit), key=lambda i: i[0])
    
    # once we know the range of gene positions close to our variant, we calculate the distance from the variant to the gene star/end point. This is necessary for the variants that dont fall inside/close to a gene.
    for start, end in sorted_positions[bisect_idx_minus: bisect_idx_plus + 1]:
        
        nearest[gtf_d_chr[start, end]] = (
                min(
                    [abs(start - pos), abs(end - pos)]
                    )
            )
        
        # here we skip further analyzing the gene position if it is outside the window range:
        if pos + limit < start - limit or pos - limit > end + limit:
            continue
        
        # if the variant is inside the range, we analyze if is indeed inside a gene or at least said gene is in the window range:
        if start <= pos and pos <= end:
            genes_in.add(gtf_d_chr[start, end])
        
        if start <= pos + limit and pos - limit <= end:
            genes_window.add(gtf_d_chr[start, end])
    
    if len(genes_window) == 0:
        genes_window = "."
    
    # if the variant did not fall inside a gene, we obtain the closest gene to said variant:
    if len(genes_in) == 0:
        genes_in = '.'
        nearest_gene = min(nearest, key=nearest.get)
    else:
        nearest_gene = '.'
    
    # we prepare and return the new variant line
    ann_info = f'GENES_IN={",".join(genes_in)};GENES_200KB={",".join(genes_window)};GENE_NEAREST={nearest_gene}\n'
    fields[header['INFO']] = ann_info
    
    ann_variant = "\t".join(fields)
    
    return ann_variant


def parse_gtf_file(gtf_file: str) -> dict:
    """ Parses a gtf file to obtain a dictionary in the form of (start, end): gene id.

    Args:
        gtf_file (str): full route to gtf file

    Returns:
        dict: dictionary in the form of (start, end): gene id
    """
    
    gtf_dict = defaultdict(dict)
    
    with open(gtf_file, 'r') as gtf_f:
        
        for line in gtf_f:
            chr, start, end, ensembl_id = line.strip().split('\t')
            
            gtf_dict[chr][(int(start), int(end))] = ensembl_id.replace('"', '').replace(';', '')
            
    return gtf_dict

def annotation_parallel(vcf_files: list, num_process: int, window_limit: int, gtf_dict: dict):
    """launches each separate process with the corresponding vcf files

    Args:
        vcf_files (list): list of vcf files (full routes)
        num_process (int): number of process to use
        window_limit (int): the up/downstream base pairs to add to the variant to study a wider window.
        gtf_dict (dict):  dictionary in the form of (start, end)
    """    
    
    # here we divide the number of files by the number of processes
    divs = (len(vcf_files)) // num_process
    batch = []
    
    # we create sublists of files to pass to each core/process
    batch = [vcf_files[x: x + divs + 1] for x in range(0, len(vcf_files), divs)]
    process_list = []
    
    # we launch each process
    for i in range(0, num_process):
        
        print(f'Launching process nº: {i + 1}')
                
        process_list.append(annotation_process(batch[i], gtf_dict, window_limit))
        process_list[i].start()
    
    # and we check and stop them
    check = 0
    for i in process_list:
        check +=1
        i.join()
        print(f'Process number {check} stopped.')

if __name__ == '__main__':
    
    pars = argparse.ArgumentParser(description= 'Circa skill test: VCF-GTF annotation')
    pars.add_argument("--vcf", 
                    help="Directory with all the VCF.gz files",
                    required=True)
    pars.add_argument("--gtf", 
                      help='gtf file compressed (gz)', required=True)
    pars.add_argument('--num_process', help= "number of cores to dedicate", required=True)
    pars.add_argument('--pb', help= "pair of bases to add to the down/upstream gene search from the target genomic position", required=True)
    args = pars.parse_args()    
    
    
    # Read vcf files directory
    vcf_files = [(f'{args.vcf}/{f}') for f in listdir(args.vcf) if re.match('.*(vcf.gz)$', f)]
    
    # Read GTF file
    gtf_dict = parse_gtf_file(args.gtf)
    
    annotation_parallel(vcf_files, int(args.num_process), int(args.pb), gtf_dict)
