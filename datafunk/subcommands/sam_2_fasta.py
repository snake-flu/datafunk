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

    trim = False
    if options.trim:
        trim = True
        if ':' in options.trim:
            s = options.trim.strip().split(':')
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
        warnings.warn('Trim argument not formatted properly. Ignoring trimming.')
    if options.trim and trim:
        if trimstart > RLEN or trimend > RLEN:
            trim = False
            warnings.warn('Trim values are larger than length of reference. Ignoring trimming.')
        if trimstart > trimend:
            trim = False
            warnings.warn('Trim argument not formatted properly. Ignoring trimming.')


    if trim:
        if options.output and not options.stdout:
            f_out = open(options.output, 'w')
            if options.prefix_ref:
                f_out.write('>' + reference.id + '\n')
                f_out.write(str(reference.seq[trimstart:trimend]) + '\n')
        if options.stdout:
            if options.prefix_ref:
                print('>' + reference.id)
                print(reference.seq[trimstart:trimend])
    else:
        if options.output and not options.stdout:
            f_out = open(options.output, 'w')
            if options.prefix_ref:
                f_out.write('>' + reference.id + '\n')
                f_out.write(str(reference.seq) + '\n')
        if options.stdout:
            if options.prefix_ref:
                print('>' + reference.id)
                print(reference.seq)


    for query_seq_name, one_querys_alignment_lines in itertools.groupby(samfile, lambda x: parse_sam_line(x)['QNAME']):
        # one_querys_alignment_lines is an iterator corresponding to all the lines
        # in the SAM file for one query sequence
        seq = get_seq_from_block(sam_block = one_querys_alignment_lines, rlen = RLEN)

        if trim:
            if options.output and not options.stdout:
                f_out.write('>' + query_seq_name + '\n')
                f_out.write(seq[trimstart:trimend] + '\n')

            if options.stdout:
                print('>' + query_seq_name)
                print(seq[trimstart:trimend])
        else:
            if options.output and not options.stdout:
                f_out.write('>' + query_seq_name + '\n')
                f_out.write(seq + '\n')

            if options.stdout:
                print('>' + query_seq_name)
                print(seq)



    if options.output and not options.stdout:
        f_out.close()
