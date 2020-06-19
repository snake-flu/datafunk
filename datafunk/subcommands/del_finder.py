from datafunk.del_finder import *

def run(options):
    if options.append_snp and not options.fasta_out:
        sys.exit("Error: If you want to append deletions as SNPs, you need to specify a fasta file to write. Use --output-fasta or -o.")

    del_finder(fasta_in = options.fasta_in,
               fasta_out = options.fasta_out,
               del_file = options.del_file,
               genotypes_file = options.genotypes_file,
               append_snp = options.append_snp)
