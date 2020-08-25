from datafunk.sam_2_fasta import *

def run(options):

    samfile = pysam.AlignmentFile(options.sam, 'r')
    reference = SeqIO.read(options.reference, 'fasta')

    # The length of the reference sequence:
    RLEN = samfile.header['SQ'][0]['LN']
    # The name of the target sequence:
    REF = samfile.header['SQ'][0]['SN']

    if REF != reference.id:
        sys.exit('reference names differ!')
    if RLEN != len(reference.seq):
        sys.exit('reference lengths differ!')


    if options.stdout or not options.output_fasta:
        output = 'stdout'
    else:
        output = options.output_fasta

    trim = False
    if options.trim:
        trim = True
        if ':' in options.trim:
            s = options.trim.strip().replace("[","").replace("]","").split(':')
            if len(s) == 2:
                if len(s[0]) == 0:
                    trimstart = 0
                else:
                    if s[0].isdigit():
                        trimstart = int(s[0])
                    else:
                        trim = False
                if len(s[1]) == 0:
                    trimend = RLEN
                else:
                    if s[0].isdigit():
                        trimend = int(s[1])
                    else:
                        trim = False
            else:
                trim = False
        else:
            trim = False

    if options.trim and not trim:
        sys.stderr.write('Trim argument not formatted properly. Ignoring trimming.')
    if options.trim and trim:
        if trimstart > RLEN - 1 or trimend > RLEN:
            trim = False
            sys.stderr.write('Trim values are larger than length of reference. Ignoring trimming.')
        if trimstart > trimend:
            trim = False
            sys.stderr.write('Trim argument not formatted properly. Ignoring trimming.')

    if trim:
        sam_2_fasta(samfile = samfile,
                    reference = reference,
                    output = output,
                    prefix_ref = options.prefix_ref,
                    log_inserts = options.log_inserts,
                    log_all_inserts = options.log_all_inserts,
                    log_dels = options.log_dels,
                    log_all_dels = options.log_all_dels,
                    trim = True,
                    pad = options.pad,
                    trimstart = trimstart,
                    trimend = trimend)
    else:
        sam_2_fasta(samfile = samfile,
                    reference = reference,
                    output = output,
                    prefix_ref = options.prefix_ref,
                    log_inserts = options.log_inserts,
                    log_all_inserts = options.log_all_inserts,
                    log_dels = options.log_dels,
                    log_all_dels = options.log_all_dels)
