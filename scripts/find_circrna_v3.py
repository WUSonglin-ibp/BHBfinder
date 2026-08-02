from Bio import SeqIO
import re
import collections
import argparse
import gzip


# Global parameters
# maximal length of intron
maxintron = 5000
# minimal length of each hit in hybrid read
minarm = 15

crtseqs = collections.defaultdict(dict)
########################################################
# New method based on blastn table (fmt 6)
########################################################
def read_genome (genomefile):
    chr = collections.defaultdict(dict)
    with open(genomefile) as g:
        for k in SeqIO.parse(g, "fasta"):
            chr[k.id] = k.seq.upper()
    return chr


def read_bedgraph(file):
    """
    Read bedgraph file, coor is zero based
    """
    cov = {}
    for l in open(file):
        a = l.rstrip().split('\t')
        cov[(a[0], int(a[1])+1)] = int(float(a[2]))
    return cov


def read_junc_reads (fasta_file):
    """
    read fasta file
    """
    fa = {}
    with open(fasta_file) as f:
        for k in SeqIO.parse(f, "fasta"):
            fa[k.id] = k.seq.upper()
    return fa


def read_gff(file):
    global gff
    gff = collections.OrderedDict()
    gene_entry_list = ['gene', 'mobile_genetic_element','pseudogene']
    escaped_annotation_list = ['region', 'exon']

    read_lines = open(file).readlines()
    k = '' # k of gff
    need_annotation_flag = False
    for i, l in enumerate(read_lines):
        if l[0] == '#' or re.search(r'^\s$',l):
            continue
        a = l.rstrip().split('\t')
        if a[2] in escaped_annotation_list:
            # avoid the whole genome notation
            continue

        gene_entry_type = a[2]
        ## Meet new gene
        if gene_entry_type in gene_entry_list:
            need_annotation_flag = True
            chro = a[0]
            start = int(a[3])
            end = int(a[4])
            strand = a[6]
            k = (chro, start, end, strand)
            gff[k] = {}
            m = re.search('ID=(.+?);', a[8])
            gff[k]['gene_ID'] = m.group(1)
            gff[k]['product'] = ''
            gff[k]['gbkey'] = gene_entry_type
            continue
        ## make annotation
        if need_annotation_flag:
            need_annotation_flag = False
            gff[k]['gbkey'] = gene_entry_type
            m = re.search('product=([^;]+);?', a[8])
            if m:
                gff[k]['product'] = m.group(1)


def read_circreads (r1file):
    if re.search('.gz$', r1file):
        f1 = gzip.open(r1file)
    else:
        f1 = open(r1file)
    for k in SeqIO.parse(f1, "fasta"):
        if args.flagr2 == True:
            crtseqs[k.id]['seq1'] = k.seq.reverse_complement()
        else:
            crtseqs[k.id]['seq1'] = k.seq
        crtseqs[k.id]['r1hitcount'] = 0
        crtseqs[k.id]['rnatype'] = 'unknown'
    f1.close()


def read_circblasttable (r1file):
    """
    read blast result.
    """
    if re.search('.gz$', r1file):
        f1 = gzip.open(r1file)
    else:
        f1 = open(r1file)
    for line in f1:
        line = line.rstrip('\n')
        w = line.split('\t')
        if len(w) > 13:
            continue
        seq_id = w[0]
        chro = w[1]	# the chromosome of reads mapped to genome
        readlength = int(w[12])     # length of read
        alignlength = int(w[3])     # length of aligned sequence including gap
        ### reverse seq if necessary
        qstart_value = int(w[6])
        qend_value = int(w[7])
        sstart_value = int(w[8])
        send_value = int(w[9])
        if args.flagr2 == True:
            qstart_value = readlength - int(w[7]) + 1
            qend_value = readlength - int(w[6]) + 1
            sstart_value = int(w[9])
            send_value = int(w[8])
        if (readlength-alignlength) > minarm:    # split alignment
            if crtseqs[seq_id]['r1hitcount'] == 0:
                crtseqs[seq_id]['q1start'] = qstart_value
                crtseqs[seq_id]['q1end'] = qend_value
                crtseqs[seq_id]['s1start'] = sstart_value
                crtseqs[seq_id]['s1end'] = send_value
                crtseqs[seq_id]['s1acc'] = chro
                crtseqs[seq_id]['r1hitcount'] += 1
            #The second aligned part should be 15 nt away from the 1st aligned part, avoid multiple alignment
            elif abs(crtseqs[seq_id]['q1start'] - qstart_value) > minarm \
             and abs(crtseqs[seq_id]['q1end'] - qend_value) > minarm\
             and abs(crtseqs[seq_id]['s1start'] - sstart_value) > minarm\
             and abs(crtseqs[seq_id]['s1end'] - send_value) > minarm\
             and (crtseqs[seq_id]['q1start'] - qstart_value) * (crtseqs[seq_id]['q1end'] - qend_value) > 0 \
             and (crtseqs[seq_id]['s1start'] - sstart_value) * (crtseqs[seq_id]['s1end'] - send_value) > 0 \
             and crtseqs[seq_id]['s1acc'] == chro:	# one hit is not allowed to fall in the other hit, the genome should be the same
                crtseqs[seq_id]['r1hitcount'] += 1
                if crtseqs[seq_id]['q1start'] < qstart_value:     # hit1 is on the left, hit2 on the right of read
                    crtseqs[seq_id]['q2start'] = qstart_value
                    crtseqs[seq_id]['q2end'] = qend_value
                    crtseqs[seq_id]['s2start'] = sstart_value
                    crtseqs[seq_id]['s2end'] = send_value
                else:    # exchange
                    crtseqs[seq_id]['q2start'] = crtseqs[seq_id]['q1start']
                    crtseqs[seq_id]['q2end'] = crtseqs[seq_id]['q1end']
                    crtseqs[seq_id]['s2start'] = crtseqs[seq_id]['s1start']
                    crtseqs[seq_id]['s2end'] = crtseqs[seq_id]['s1end']
                    crtseqs[seq_id]['q1start'] = qstart_value
                    crtseqs[seq_id]['q1end'] = qend_value
                    crtseqs[seq_id]['s1start'] = sstart_value
                    crtseqs[seq_id]['s1end'] = send_value
                crtseqs[seq_id]['s2acc'] = chro
        else: # labeled if at least one alignment is normal
            crtseqs[seq_id]['rnatype'] = 'normal'
            # normal_reads.
    f1.close()


def read_circout(f):
    ''' read circ.out file from action = circ
    '''
    circinfo = collections.OrderedDict()
    for l in open(f):
        if l[0] == '#':
            continue
        l = l.rstrip()
        a = l.split('\t')
        k = a[0]	# the name of read
        if k not in circinfo.keys():
            circinfo[k] = {}
        circinfo[k]['chro'] = a[1]
        circinfo[k]['type'] = a[3]
        circinfo[k]['strandness'] = a[4]
        circinfo[k]['s1start'] = int(a[8]) # start > end in cis strand
        circinfo[k]['s1end'] = int(a[9])
        circinfo[k]['s2start'] = int(a[12])
        circinfo[k]['s2end'] = int(a[13])
        circinfo[k]['overlap'] = int(a[14])
        circinfo[k]['line'] = l # original infomation of line
    return circinfo


def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "N": "N", "Y":"Y"}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement[base] for base in reverse_sequence)
    return reverse_complement_sequence


def output_circ ():
    # for circ.out and circ.bed, seq1 is full coor, seq2 is defect coor. All overlap pair is give to seq1
    out_fh = open('%s_circ.out' % (basename), 'w')	# the mapping res of each reads, reads have 2 or more hits, but including hybrid.
    bed_fh = open('%s_circ.bed' % (basename), 'w')	# Info of reads mapped to genome much like splicing in bed formate.
    r1_fh = open('%s_circ_R1.fasta'%(basename), 'w')
    r2_fh = open('%s_circ_R2.fasta'%(basename), 'w')
    out_fh.write("#readname\tChro\tRead length\tType\tStrandness\tIntron length\tq1 start\tq1 end\ts1 start\ts1end\t"
                             "q2 start\tq2 end\ts2 start\ts2 end\tOverlap\tJunction seq\n")
    for k, v in crtseqs.items():
        if (v['rnatype'] != 'normal' and v['r1hitcount'] == 2) :
            dist_err = v['q1end'] - v['q2start']+1
            overlap = 0
            junc = v['q1end']
            if dist_err >= 1:
                overlap = dist_err
                v['q2start'] = v['q2start'] + overlap # q2 correct
                inputseq = v['seq1'][junc-10:junc] + "|" + v['seq1'][junc:junc+10]
            else:
                inputseq = v['seq1'][junc-10:junc] + "|" + v['seq1'][junc:junc - dist_err].lower() +"|" + v['seq1'][junc - dist_err:junc - dist_err+10]
            if v['s1start'] < v['s1end'] and v['s2start'] < v['s2end']:
                v['strandness'] = "++"
                v['s2start'] = v['s2start'] + overlap    # hit2 corrected
                rnalength = abs(v['s2start'] - v['s1end'] - 1)    # length of intron or circRNA
                if abs(rnalength) < maxintron:
                    if v['s1end'] < v['s2start']:
                        v['rnatype'] = 'exon'    # joined exons
                        # 5' arm is full, 3' arm is reduced
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s1start'] - 1, v['s1end'], k + '_E1', rnalength, '+'))
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s2start'] - 1, v['s2end'], k + '_E2', rnalength, '+'))
                        r1_fh.write('>%s\n%s\n' % (k + '_E', v['seq1'][v['q1start'] - 1:v['q1end']]))
                        r2_fh.write('>%s\n%s\n' % (k + '_E', v['seq1'][v['q2start'] - 1:v['q2end']]))
                    elif v['s1start'] > v['s2end']:
                        v['rnatype'] = 'circ'
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s1start'] - 1, v['s1end'], k + '_C1', rnalength, '+'))
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s2start'] - 1, v['s2end'], k + '_C2', rnalength, '+'))
                        r1_fh.write('>%s\n%s\n' % (k + '_C', v['seq1'][v['q1start'] - 1:v['q1end']]))
                        r2_fh.write('>%s\n%s\n' % (k + '_C', v['seq1'][v['q2start'] - 1:v['q2end']]))
                    else:
                        # This part is right, no need for checking again
                        if v['s1end'] <= v['s2end']:
                            v['rnatype'] = 'overlap2'
                        elif v['s1start'] <= v['s2end']:
                            v['rnatype'] = 'overlap3'
                else:
                    v['rnatype'] = 'hybrid'
            elif v['s1start'] > v['s1end'] and v['s2start'] > v['s2end']:
                v['strandness'] = "--"
                v['s2start'] = v['s2start'] - overlap    # hit2 corrected
                rnalength = abs(v['s1end'] - v['s2start'] - 1)    # length of intron (+) or circRNA (-)
                if rnalength < maxintron:
                    if v['s1end'] > v['s2start']:
                        v['rnatype'] = 'exon'
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s1end'] - 1, v['s1start'], k + '_E1', rnalength, '-'))
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s2end'] - 1, v['s2start'], k + '_E2', rnalength, '-'))
                        r1_fh.write('>%s\n%s\n' % (k + '_E', v['seq1'][v['q1start'] - 1:v['q1end']]))
                        r2_fh.write('>%s\n%s\n' % (k + '_E', v['seq1'][v['q2start'] - 1:v['q2end']]))
                    elif v['s1start'] < v['s2end']:
                        v['rnatype'] = 'circ'
                        # offset -1 right, the range is greedy
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s1end'] - 1, v['s1start'], k + '_C1', rnalength, '-'))
                        bed_fh.write(
                            '%s\t%d\t%d\t%s\t%s\t%s\n' % (v['s1acc'], v['s2end'] - 1, v['s2start'], k + '_C2', rnalength, '-'))
                        r1_fh.write('>%s\n%s\n' % (k + '_C', v['seq1'][v['q1start'] - 1:v['q1end']]))
                        r2_fh.write('>%s\n%s\n' % (k + '_C', v['seq1'][v['q2start'] - 1:v['q2end']]))
                        if v['s1end'] > v['s2end']:
                            v['rnatype'] = 'overlap2'
                        elif v['s1start'] > v['s2end']:
                            v['rnatype'] = 'overlap3'
                else:
                    v['rnatype'] = 'hybrid'
            elif v['s1start'] < v['s1end'] and v['s2start'] > v['s2end']:
                rnalength = 0
                v['strandness'] = "+-"
                v['rnatype'] = 'hybrid'
            elif v['s1start'] > v['s1end'] and v['s2start'] < v['s2end']:
                rnalength = 0
                v['strandness'] = "-+"
                v['rnatype'] = 'hybrid'
            out_fh.write('%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n' % (k, v['s1acc'], len(v['seq1']),v['rnatype'],\
             v['strandness'], rnalength, v['q1start'],v['q1end'],v['s1start'], v['s1end'],v['q2start'],v['q2end'], v['s2start'],\
             v['s2end'], overlap, inputseq))
    r1_fh.close()
    r2_fh.close()
    out_fh.close()
    bed_fh.close()

def analyze_circout():
    '''
    Analysis BHB junction
    Re output read1 and read2 with BHB split if score > 7, if score <7, reads were also output in fasta file.
    '''
    juncinfo = collections.OrderedDict()
    # If BHB cutoff is satisfied, the region of junction is modified
    r1_fh = open(basename + '_BHB_R1.fasta', 'w')
    r2_fh = open(basename + '_BHB_R2.fasta', 'w')
    for read, v in circinfo.items():
        if v['type'] not in ('circ', 'exon'):
            continue
        s1end_full = v['s1end']
        s2start_defect = v['s2start']
        strand = ''
        # creat one juction with no split
        if v['strandness'] == "++":
            strand = '+'
        elif v['strandness'] == "--":
            strand = '-'
        # save junction infomation
        upcoor = min(s1end_full, s2start_defect)	# the boundary coor of intron in upstream genome coor
        downcoor = max(s1end_full, s2start_defect)
        # No matter exon or intron, we note all splicing by intron start and end.
        if v['type'] == 'exon':
            upcoor += 1
            downcoor -= 1

        # shift the split and compare the scores
        # judging standard
        # 1. 5(1'/bp) - 4(2/bp) - 5(1'/bp)
        # 2. at least 3 bp in middle , s2start_defect)    # down stream coor
        #         if v['type'] == 'exon':helix
        scores = []
        flags = []
        junc_structures = []
        for shift in range(v['overlap']+1):
            if v['strandness'] == "++":
                junc_structure = analyze_junc(s1end_full-shift, s2start_defect-shift, v['strandness'], v['type'], v['chro'])
            elif v['strandness'] == "--":
                junc_structure = analyze_junc(s1end_full+shift, s2start_defect+shift, v['strandness'], v['type'], v['chro'])
            if junc_structure == None:
                break
            junc_structures.append(junc_structure)
        if not junc_structures:
            # Reads at the tail of genome seq
            continue
        for j in junc_structures:
            scores.append(j[-1])
            flag1 = j[1][:5].count('|') + j[1][:5].count('o')
            flag2 = j[1][8:12].count('|') + j[1][8:12].count('o')
            flag3 = j[1][15:].count('|') + j[1][15:].count('o')
            flags.append([flag1,flag2,flag3])
        maxscore = scores[0]
        maxindex = 0 # start from 0
        for i in range(len(scores)):
            if scores[i] > maxscore:
                maxindex = i
                maxscore = scores[i]
                #use tuple as key, (chromosome, upcoor, downcoor, strandness, type)
        if strand == "+":
            upcoor -= maxindex
            downcoor -= maxindex
        else:
            upcoor += maxindex
            downcoor += maxindex
        k = (v['chro'], upcoor, downcoor, strand, v['type'])
        if upcoor >= downcoor:
            continue # Avoid aligning reads continuously to the genome, but there are some parts in the middle of the reads that cannot be aligned.
        suffix = '_C' if v['type'] == 'circ' else '_E'

        r1_length = len(r1[read + suffix])
        r1_fh.write('>' + read + suffix + '\n')
        r1_fh.write(str(r1[read+suffix][: r1_length-maxindex]) + '\n')
        r2_fh.write('>' + read + suffix + '\n')
        r2_fh.write(str(r1[read+suffix][r1_length-maxindex :]) + str(r2[read+suffix][:]) + '\n')
        if k not in juncinfo.keys():
            juncinfo[k] = {}
            juncinfo[k]['chro'] = v['chro']
            juncinfo[k]['upcoor'] = upcoor
            juncinfo[k]['downcoor'] = downcoor
            juncinfo[k]['readcount'] = 1
            juncinfo[k]['nodup_readcount'] = 1
            juncinfo[k]['type'] = v['type']
            juncinfo[k]['strand'] = strand
            juncinfo[k]['overlap'] = v['overlap'] # overlap can be multi number, we only choose the first num
            juncinfo[k]['s1start_s2end'] = [(v['s1start'], v['s2end'])]
            juncinfo[k]['splitsize'] = maxindex
            juncinfo[k]['score'] = maxscore
            juncinfo[k]['seq1'] = junc_structures[maxindex][0]
            juncinfo[k]['pairing'] = junc_structures[maxindex][1]
            juncinfo[k]['seq2'] = junc_structures[maxindex][2]
            juncinfo[k]['flag1'] = flags[maxindex][0]
            juncinfo[k]['flag2'] = flags[maxindex][1]
            juncinfo[k]['flag3'] = flags[maxindex][2]
            juncinfo[k]['intron seq'] = chr[juncinfo[k]['chro']][juncinfo[k]['upcoor']-1: juncinfo[k]['downcoor']]
            juncinfo[k]['signal_rate'] = 0
            if juncinfo[k]['strand'] == '-':
                juncinfo[k]['intron seq'] = reverse_complement(juncinfo[k]['intron seq'])
            juncinfo[k]['gbkey'] = ""
            juncinfo[k]['product'] = ""
            juncinfo[k]['gene_ID'] = ""
            if if_output_splicing_reads:
                juncinfo[k]['reads'] = read
        else:
            juncinfo[k]['readcount'] += 1
            if if_output_splicing_reads:
                juncinfo[k]['reads'] += ',' + read
            if (v['s1start'], v['s2end']) not in juncinfo[k]['s1start_s2end']:
                juncinfo[k]['nodup_readcount'] += 1
                juncinfo[k]['s1start_s2end'].append((v['s1start'], v['s2end']))
    return juncinfo


def cal_spliced_rate_new():
    '''The missing value of normal cov = 0
    The ligation rate is calculated as 2*h/(2*h + n5 + n3) for individual datasets, where h is number of hybrid reads
    and n5 and n3 are minimal number of normal reads in a ±3 nt range around 5' and 3' ligation sites, respectively.
    Such counting of normal reads aims to pick bona fide normal reads and also address uncertainty in  ligation site
    determination due to the presence of overlap. Both hybrid and normal reads were counted prior to deduplication.
    Average of non-zero ligation rates from all datasets were reported.
    '''
    for k in juncinfo.keys():
        juncinfo[k]['spliced_rate_new'] = 0
        juncinfo[k]['normal_cov'] = 0
        if juncinfo[k]['nodup_readcount'] > 1:
            n5_covs = [genomecov.get((k[0], i), 0) for i in range(k[1]-3, k[1]+3)]
            n3_covs = [genomecov.get((k[0], i), 0) for i in range(k[2]-2, k[2]+4)]
            n5, n3 = min(n5_covs), min(n3_covs)
            juncinfo[k]['normal_cov'] = int((n5 + n3)/2)
            juncinfo[k]['spliced_rate_new'] = 2 * juncinfo[k]['readcount'] / (2 * juncinfo[k]['readcount'] + n5 + n3)


def output_junc(gff_file, bedgraph):
    bed_fh = open("%s_junction.bed" %(basename), 'w')
    l = '\t'.join(['#chro','start','end','name','readcount(nodup)','strand','type','length','readcount','maxscore','overlap',\
        'splitsize', 'helix1_count', 'helix2_count', 'helix3_count','seq1', 'pairing', 'seq2', 'intron seq'])
    if bedgraph != None:
        l += '\t' + 'normal_cov' + '\t' + 'ligation_rate'
    if gff_file != None:
        l += '\t' + 'gene_ID' + '\t' + 'gene_type'+ '\t' + 'gene_name'
    l += '\tSignal_rate\tAbundance'
    if if_output_splicing_reads:
        l += '\tReads'

    bed_fh.write(l+'\n')
    count = 0
    juncinfo_keys_coor_sorted = sorted(juncinfo.keys(), key = lambda x: juncinfo[x]['upcoor'])
    for k in sorted(juncinfo_keys_coor_sorted, key = lambda x: juncinfo[x]['chro']):
        v = juncinfo[k]
        count += 1
        juncinfo[k]['name'] = 'junc' + str(count)
        if juncinfo[k]['type'] == 'circ':
            juncinfo[k]['name'] += '_C'
        elif juncinfo[k]['type'] == 'exon':
            juncinfo[k]['name'] += '_E'
        bed_fh.write('%s\t%d\t%d\t%s\t%d\t%s\t' % (v['chro'], v['upcoor'] - 1, v['downcoor'], v['name'], v['nodup_readcount'], v['strand']))
        bed_fh.write('%s\t%d\t%d\t%.1f\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s'%
                     (v['type'],v['downcoor'] -v['upcoor']+1,v['readcount'], v['score'],v['overlap'],v['splitsize'],
                      v['flag1'],v['flag2'],v['flag3'], v['seq1'], v['pairing'],v['seq2'], v['intron seq']))
        if bedgraph != None:
            bed_fh.write('\t%d\t%.4f' % (v['normal_cov'], v['spliced_rate_new']))
        if gff_file != None:
            bed_fh.write('\t' + v['gene_ID'] + '\t' + v['gbkey'] + '\t' + v['product'])
        bed_fh.write('\t' + str(round(v['signal_rate'], 4)) + '\t0' )
        if if_output_splicing_reads:
            bed_fh.write('\t' + v['reads'])
        bed_fh.write('\n')
    bed_fh.close()


def analyze_junc(s1end, s2start, strandness, rnatype,chromosome):
    seq5 = '' # from 5' to 3' in RNA
    seq3 = '' # from 5' to 3' in RNA
    ### fixing coordiante
    # trans 1-based to 0-based, so -1 for s1end and s2start
    # the sequence range can not be out of the choromosome range!
    if rnatype == 'exon' and strandness == '++':
        if s1end-1-10 < 0 or s2start-1+6 > len(chr[chromosome]) :
            return None
        seq5 = chr[chromosome][s1end-1-10:s1end-1+7]
        seq3 = chr[chromosome][s2start-1-11:s2start-1+6]
    elif rnatype == 'exon' and strandness == '--':
        if s1end-1+11 > len(chr[chromosome]) or s2start-1-5 < 0:
            return None
        seq5 = chr[chromosome][s1end-1-6: s1end-1+11].reverse_complement()
        seq3 = chr[chromosome][s2start-1-5:s2start-1+12].reverse_complement()
    elif rnatype == 'circ' and strandness == '++':
        if s2start-1-11 < 0 or s1end-1+7 > len(chr[chromosome]):
            return None
        seq5 = chr[chromosome][s2start-1-11:s2start-1+6]
        seq3 = chr[chromosome][s1end-1-10:s1end-1+7]
    elif rnatype == 'circ' and strandness == '--':
        if s2start-1-5 < 0 or s2start-1+12 > len(chr[chromosome]) or s1end-1-6 < 0 or s1end-1+11 > len(chr[chromosome]) :
            return None
        seq5 = chr[chromosome][s2start-1-5: s2start-1+12].reverse_complement()
        seq3 = chr[chromosome][s1end-1-6: s1end-1+11].reverse_complement()
    #print seq5, seq3
    align5 = seq5[:5] + "---" +seq5[5:]
    align3 = seq3[:5] + "---" + seq3[5:]
    align3r = align3[::-1]
    basepair = list(" "*20)
    for i in (0,1,2,3,4,8,9,10,11,15,16,17,18,19):
        basepair[i] = checkpair(align5[i], align3r[i])
    bp_central = "".join(basepair[8:12])
    bp_central_score = 2 * bp_central.count("|") + 0.5 * bp_central.count("o")
    bp_up = "".join(basepair[0:5])
    bp_up_score = bp_up.count("|") + 0.5 * bp_up.count("o")
    bp_down = "".join(basepair[15:20])
    bp_down_score = bp_down.count("|") + 0.5 * bp_down.count("o")
    bp_all = "".join(basepair)
    bp_total_score = bp_up_score + bp_central_score + bp_down_score
    return align5, bp_all, align3r, bp_total_score


def add_note(if_anti = False):
    for kjunc, vjunc in juncinfo.items():
        for kgff, vgff in gff.items():
            if kjunc[0] == kgff[0] and is_overlap(kjunc[1], kjunc[2], kgff[1], kgff[2]):
                for key in ('gbkey', 'product', 'gene_ID'):
                    if key in vgff.keys():
                        if kjunc[3] == kgff[3]:
                            vjunc[key] += vgff[key] + ','
                        elif kjunc[3] != kgff[3] and if_anti:
                            if key == 'gbkey':
                                vjunc[key] += 'anti ' + vgff[key] + ','
                            else:
                                vjunc[key] += vgff[key] + ','
        vjunc['gbkey'] = vjunc['gbkey'][:-1]
        vjunc['product'] = vjunc['product'][:-1]
        vjunc['gene_ID'] = vjunc['gene_ID'][:-1]


def is_overlap(s1, e1, s2, e2):
    if max(s1, s2) <= min(e1, e2):
        return True
    else:
        return False


def checkpair(b1, b2):
    comb = b1 + b2
    if comb in ("GC", "CG", "AT", "TA"):
        return "|"
    elif comb in ("GT", "TG"):
        return "o"
    else:
        return " "


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Circular RNA ")
    parser.add_argument('-a', '--action', required=True, choices=['circ', 'junction'])
    parser.add_argument('-flagr2', '--flagr2', help="flag for forward or reverse reads r1 or r2", action='store_true', default=False)
    parser.add_argument('-r1', '--r1', help="read 1, fasta or fastq")
    parser.add_argument('-j1', help="Junction read1")
    parser.add_argument('-j2', help="Junction read2")
    parser.add_argument('-g', '--genome', help="fasta file for genome")
    parser.add_argument('-c', '--circ', help="read circ result from action = circ")
    parser.add_argument('--output_reads', action='store_true',  help="Output the spliced reads if need default true")
    parser.add_argument('-basename', '--basename', help="basename of files")
    parser.add_argument('-tab1', '--table_blast1', help="blast table fmt 6")
    parser.add_argument('-s', '--score', default=9, help=">= score,default=5", )
    parser.add_argument('-gff', '--gff', default=None, help="read genome GFF or GFF3 file, default = None")
    parser.add_argument('-bam', '--bam',default = None, help="normal alignment sorted bam file, default = None")
    parser.add_argument('-j', '--junction', help = 'junction.bed file from action = junction')
    parser.add_argument('-bedgraph', '--bedgraph', help='Bedgraph file for normal alignment')
    parser.add_argument('--bed', help='.bed file')
    parser.add_argument('--anti',  action='store_true', help='if anti notification is need')
    parser.add_argument('-o', '--output', default='result.txt', help='Output file name')

    args = parser.parse_args()
    basename = args.basename
    if args.action == "circ":
        # Output junction reads.
        read_circreads(args.r1)
        read_circblasttable(args.table_blast1)
        output_circ()

    if args.action == "junction":
        if_output_splicing_reads = args.output_reads
        chr = read_genome(args.genome)
        circinfo = read_circout(args.circ)
        genomecov = read_bedgraph(args.bedgraph)
        r1 = read_junc_reads(args.j1)
        r2 = read_junc_reads(args.j2)
        juncinfo = analyze_circout() # analyze BHB motif, and output read fasta file
        if args.gff != None:
            read_gff(args.gff)
            add_note(args.anti)
        cal_spliced_rate_new()
        output_junc(args.gff, args.bedgraph)
