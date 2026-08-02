import collections
import re
import argparse

SEARCH_RANGE = 10
SEARCH_CUTOFF = 0


def simpleSpecies(species):
    return species[0] + species[2:5]


def str2bool(value):
    if isinstance(value, bool):
        return value
    return str(value).lower() in ['true', 't', 'yes', 'y', '1']


def colname(name):
    return name.replace('#', '').strip().lower()


def get_col(header, names):
    header_dict = {colname(v): i for i, v in enumerate(header)}
    for name in names:
        if colname(name) in header_dict:
            return header_dict[colname(name)]
    return None


def get_value(w, header, names, default=''):
    idx = get_col(header, names)
    if idx is None or idx >= len(w):
        return default
    return w[idx]


def newscore (s1, s2, p):
    """score for BHB motifs"""
    def bp(b1, b2):
        # score for a base pair
        score = 0
        if b1 == "A" and b2 == "T":
            score = 1
        if b1 == "T" and b2 == "A":
            score = 1
        if b1 == "G" and b2 == "C":
            score = 1.5
        if b1 == "C" and b2 == "G":
            score = 1.5
        if b1 == "G" and b2 == "T":
            score = 1
        if b1 == "T" and b2 == "G":
            score = 1
        return score

    p=p.replace('o','|')
    score_1 = searchall('||', p[0:5]) + bp(s1[0],s2[0]) + bp(s1[1],s2[1]) + bp(s1[2],s2[2]) + 1.5*bp(s1[3],s2[3]) +2*bp(s1[4],s2[4])
    score_2 = 2*searchall('||', p[8:12]) + 2*bp(s1[8],s2[8]) + 2*bp(s1[9],s2[9]) + 2*bp(s1[10],s2[10]) + 2*bp(s1[11],s2[11])
    score_3 = searchall('||', p[15:20]) + 2*bp(s1[15],s2[15]) + 1.5*bp(s1[16],s2[16]) + bp(s1[17],s2[17]) + bp(s1[18],s2[18]) + bp(s1[19],s2[19])
    return max(score_1, score_3) + score_2


def searchall (query, template):
    len_query = len(query)
    len_template = len(template)
    count = 0
    for i in range(len_template - len_query + 1):
      if template[i:i+len_query] == query:
          count +=1
    return count


def read_input (inputfile, species):
    # read processing data from Wu SL.
    name_count = 0
    with open(inputfile) as f:
        header = next(f).rstrip('\n').split('\t')

        for line in f:
            line = line.replace('"', '')
            w = line.rstrip('\n').split('\t')
            name_count += 1
            seqid = get_value(w, header, ['Name'])
            crtseqs[seqid]['name'] = species + '_' +  str(name_count)
            crtseqs[seqid]['chr'] = get_value(w, header, ['Chr'])
            crtseqs[seqid]['start'] = int(get_value(w, header, ['Start'])) # zero based
            crtseqs[seqid]['end'] = int(get_value(w, header, ['End']))
            crtseqs[seqid]['junc'] = seqid[5:]
            crtseqs[seqid]['count'] = int(get_value(w, header, ['Total_nodup_count']))
            crtseqs[seqid]['strand'] = get_value(w, header, ['Strandness'])
            crtseqs[seqid]['rna_type'] = get_value(w, header, ['Type']).replace('linear', 'exon')
            crtseqs[seqid]['length'] = int(get_value(w, header, ['Length']))
            crtseqs[seqid]['maxscore'] = float(get_value(w, header, ['maxscore']))
            crtseqs[seqid]['overlap'] = int(get_value(w, header, ['Overlap']))
            crtseqs[seqid]['split_size'] = int(get_value(w, header, ['Splitsize']))
            crtseqs[seqid]['gene_id'] = get_value(w, header, ['Gene_ID'])
            crtseqs[seqid]['gene_type'] = get_value(w, header, ['Gene_type'])
            crtseqs[seqid]['gene_product'] = get_value(w, header, ['Gene_product'])
            crtseqs[seqid]['h1'] = int(get_value(w, header, ['Helix1_count']))
            crtseqs[seqid]['h2'] = int(get_value(w, header, ['Helix2_count']))
            crtseqs[seqid]['h3'] = int(get_value(w, header, ['Helix3_count']))
            crtseqs[seqid]['seq1'] = get_value(w, header, ['Seq1'])
            pairing = get_value(w, header, ['Pairing'])
            crtseqs[seqid]['pairing'] = pairing.replace(' ', 'x')
            crtseqs[seqid]['seq2'] = get_value(w, header, ['Seq2'])
            crtseqs[seqid]['ligation_rate_average'] = float(get_value(w, header, ['ligation_rate_average', 'Ave_spliced_rate'], 0))
            crtseqs[seqid]['spliced_rate'] = crtseqs[seqid]['ligation_rate_average']
            crtseqs[seqid]['ligation_rate_sd'] = float(get_value(w, header, ['ligation_rate_SD', 'SD_spliced_rate'], 0))
            #crtseqs[seqid]['overlap_seq'] = get_value(w, header, ['Overlap_seq'])
            crtseqs[seqid]['start_raw'] = 0 # make start bigger, not the star before split
            crtseqs[seqid]['end_raw'] = 0
            if  crtseqs[seqid]['strand'] == '+':
                crtseqs[seqid]['start_raw'] = crtseqs[seqid]['start'] + crtseqs[seqid]['split_size']
                crtseqs[seqid]['end_raw'] = crtseqs[seqid]['end'] + crtseqs[seqid]['split_size']
            elif crtseqs[seqid]['strand'] == '-':
                crtseqs[seqid]['start_raw'] = crtseqs[seqid]['start'] - crtseqs[seqid]['split_size'] + crtseqs[seqid]['overlap']
                crtseqs[seqid]['end_raw'] = crtseqs[seqid]['end'] - crtseqs[seqid]['split_size'] + crtseqs[seqid]['overlap']

            crtseqs[seqid]['length_target'] = -1
            crtseqs[seqid]['length_minus'] = -1
            crtseqs[seqid]['length_plus'] = -1
            crtseqs[seqid]['length_exon'] = 0
            crtseqs[seqid]['newscore'] = newscore(crtseqs[seqid]['seq1'], crtseqs[seqid]['seq2'], pairing)

            # define intergenetic RNA
            if crtseqs[seqid]['gene_type'] == '':
                crtseqs[seqid]['gene_type'] = 'intergene'

            # Reclassification of gene_type into 4 rna_class
            if crtseqs[seqid]['gene_type'] in ['CDS', 'CDS,CDS']:
                crtseqs[seqid]['rna_class'] = 'CDS'
            elif crtseqs[seqid]['gene_type'] in ['rRNA', 'rRNA,rRNA', 'rRNA,tRNA,rRNA', 'tRNA,rRNA', 'rRNA,tRNA']:
                crtseqs[seqid]['rna_class'] = 'rRNA'
            elif crtseqs[seqid]['gene_type'] in ['tRNA']:
                crtseqs[seqid]['rna_class'] = 'tRNA'
            else:
                crtseqs[seqid]['rna_class'] = 'other'
            # reassign 7S, 5S
            for q in ['7S', '5S']:
                if re.search(q, crtseqs[seqid]['gene_product']):
                    crtseqs[seqid]['rna_class'] = 'other'


def giveBhbFlag_250827():
    # Modify rule 2 in BHB
    # all antisense as LC

    for k, v in crtseqs.items():
        v['BHB'] = 'N'
        if v['rna_class'] == 'anti':
            v['BHB'] = 'LC'
            continue

        if v['maxscore'] >= 10 and v['newscore'] >= 20 and v['h2'] >= 3:
            if v['maxscore'] >= 15  or \
                (v['newscore'] >= 23 and v['maxscore'] >= 11 and v['overlap'] <=3 and v['count'] >= 4*v['length_minus'] and v['count'] >= 4*v['length_plus']) or \
                (v['count'] >= 10 and v['length_minus'] >= 2 and v['count'] > 2*v['length_minus'] and v['length_minus'] >= 4*v['length_plus']) or \
                (v['length_exon'] >= 2 and v['overlap'] <= 3):
                v['BHB'] = "Y"
                continue

        if  (v['overlap'] >=3)  or (v['spliced_rate'] >= 0.8 and v['rna_type'] in ['linear', 'exon']) :
            v['BHB'] = 'LC'
            


def giveAntijuncFlag(is_specific):
    if is_specific:
        junc_list = sorted(crtseqs.keys(), key=lambda x: (crtseqs[x]['chr'], crtseqs[x]['start_raw'], crtseqs[x]['end_raw'], crtseqs[x]['rna_type'] , crtseqs[x]['strand']))
        junc_len = len(junc_list)
        count_pair = 0
        for k in range(junc_len - 1):
            # antisense circ or exon
            if (crtseqs[junc_list[k + 1]]['strand'] != crtseqs[junc_list[k]]['strand'] and
                    crtseqs[junc_list[k + 1]]['count'] > 1 and
                    crtseqs[junc_list[k]]['count'] > 1 and
                    crtseqs[junc_list[k + 1]]['start_raw'] == crtseqs[junc_list[k]]['start_raw'] and
                    crtseqs[junc_list[k + 1]]['end_raw'] == crtseqs[junc_list[k]]['end_raw'] and
                    crtseqs[junc_list[k + 1]]['rna_type'] == crtseqs[junc_list[k]]['rna_type']):
                count_pair += 1
                if crtseqs[junc_list[k + 1]]['count'] > crtseqs[junc_list[k]]['count']:
                    crtseqs[junc_list[k]]['rna_class'] = 'anti'
                else:
                    crtseqs[junc_list[k+1]]['rna_class'] = 'anti'
    ### unspeicific data
    else:
        for k, v in crtseqs.items():
            gene_types = v['gene_type'].split(',')
            if gene_types:
                anti_flag = True
                for temp_type in gene_types:
                    if 'anti ' not in temp_type:
                        anti_flag = False
                        break
                if anti_flag:
                    v['rna_class'] = 'anti'


def exon_pair_nooutput ():
    # label circ/exon pair
    junc_list = sorted(crtseqs.keys(), key=lambda x: (crtseqs[x]['chr'], crtseqs[x]['start'], crtseqs[x]['end'], crtseqs[x]['strand']))
    junc_len = len(junc_list)
    count_pair = 0
    #for i, v in enumerate(junc_list):
    for k in range(junc_len-1):
        if (crtseqs[junc_list[k+1]]['rna_type'] != crtseqs[junc_list[k]]['rna_type'] and
                crtseqs[junc_list[k + 1]]['chr'] == crtseqs[junc_list[k]]['chr'] and
                crtseqs[junc_list[k+1]]['strand'] == crtseqs[junc_list[k]]['strand'] and
                crtseqs[junc_list[k+1]]['start'] == crtseqs[junc_list[k]]['start'] and
                crtseqs[junc_list[k+1]]['end'] == crtseqs[junc_list[k]]['end']):
            count_pair += 1
            crtseqs[junc_list[k]]['length_exon'] = crtseqs[junc_list[k+1]]['count']
            crtseqs[junc_list[k+1]]['length_exon'] = crtseqs[junc_list[k]]['count']
            #print ('%d\t%s\t%d\t%d\t%s\t%d' % (k, junc_list[k], crtseqs[junc_list[k]]['length_exon'],k+1, junc_list[k+1], crtseqs[junc_list[k+1]]['length_exon']))


def search_near ():
    # find a hit, search neighbor circRNAs with shorter intron formed after degradation and ligation.
    count_circ = 0
    junc_list = sorted(crtseqs.keys(), key=lambda x: (crtseqs[x]['start'],crtseqs[x]['end']))
    for i, seqid in enumerate(junc_list):
        if crtseqs[seqid]['count'] > SEARCH_CUTOFF :
            target_start = crtseqs[seqid]['start']
            target_end = crtseqs[seqid]['end']
            target_length = crtseqs[seqid]['length']
            length_minus = 0
            length_plus = 0
            length_target = 0
            count_circ += 1
            #search position +/- 100
            for j in range(max(i-100,0), min(i+100,len(junc_list))):
                if (crtseqs[junc_list[j]]['start'] in range(target_start-SEARCH_RANGE,target_start+SEARCH_RANGE) and
                    crtseqs[junc_list[j]]['end'] in range(target_end-SEARCH_RANGE,target_end+SEARCH_RANGE) and
                    crtseqs[junc_list[j]]['strand'] == crtseqs[seqid]['strand']
                    ):
                    if crtseqs[seqid]['rna_type'] == 'circ':
                        if crtseqs[junc_list[j]]['rna_type'] == 'circ' and crtseqs[junc_list[j]]['length'] < target_length:
                            length_minus += crtseqs[junc_list[j]]['count']
                        elif crtseqs[junc_list[j]]['rna_type'] == 'circ' and crtseqs[junc_list[j]]['length'] > target_length:
                            length_plus += crtseqs[junc_list[j]]['count']
                        elif crtseqs[junc_list[j]]['rna_type'] == 'circ' and crtseqs[junc_list[j]]['length'] == target_length:
                            length_target += crtseqs[junc_list[j]]['count']
                    else:
                        if crtseqs[junc_list[j]]['rna_type'] in ['exon', 'linear'] and crtseqs[junc_list[j]]['length'] > target_length:
                            length_minus += crtseqs[junc_list[j]]['count']
                        elif crtseqs[junc_list[j]]['rna_type'] in ['exon', 'linear'] and crtseqs[junc_list[j]]['length'] < target_length:
                            length_plus += crtseqs[junc_list[j]]['count']
                        elif crtseqs[junc_list[j]]['rna_type'] in ['exon', 'linear'] and crtseqs[junc_list[j]]['length'] == target_length:
                            length_target += crtseqs[junc_list[j]]['count']
            crtseqs[seqid]['length_target'] = length_target - crtseqs[seqid]['count']
            crtseqs[seqid]['length_minus'] = length_minus
            crtseqs[seqid]['length_plus'] = length_plus


def outputRc2(output_file_name):
    # high confidence and interesting hits
    out_all = open(output_file_name, 'w')
    head_line = '\t'.join([
        'Species', 'chr', 'start', 'end', 'junc_ID', 'j_main', 'strand', 'RNA_type', 'length',
        'BHB_score', 'overlap', 'H1', 'H2', 'H3', 'BHB', 'j_shift', 'j_minus', 'j_plus',
        'j_pair', 'ligation_rate_average', 'ligation_rate_SD', 'seq1', 'pairing', 'seq2',
        'RNA_class', 'gene_class', 'gene_ID', 'gene_product', 'New_score'
    ])
    out_all.write(head_line + '\n')
    for seqid in sorted(crtseqs.keys(), key=lambda x: (crtseqs[x]['count']), reverse=True):
        v = crtseqs[seqid]
        print_line = '\t'.join(map(str, (
            species[0] + '. ' + species[2:], v['chr'], v['start'], v['end'], v['name'], v['count'],
            v['strand'], v['rna_type'].replace('exon', 'linear'), v['length'], v['maxscore'],
            v['overlap'], v['h1'], v['h2'], v['h3'], v['BHB'], v['length_target'], v['length_minus'],
            v['length_plus'], v['length_exon'], round(v['ligation_rate_average'], 4),
            round(v['ligation_rate_sd'], 4), v['seq1'], v['pairing'], v['seq2'], v['rna_class'],
            v['gene_type'], v['gene_id'], v['gene_product'], v['newscore']
        )))
        if (v['count'] > 1 and v['overlap'] < 10):
            out_all.write (print_line + '\n')
    out_all.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get Set2 (read count > 1, overlap < 10 nt)")
    parser.add_argument('--species', type=str,  required=True, help="Bed files for merge")
    parser.add_argument('--set1', type=str, required=True, help="Bed files for merge")
    parser.add_argument('--specific', default='True', help='If the data is strand specific data, specific data by default' )
    args = parser.parse_args()

    output_file_name = args.species + '_Set2.bed'
    species = args.species
    simple_species = simpleSpecies(species)
    set1_file = args.set1
    is_specific = str2bool(args.specific)

    crtseqs = collections.defaultdict(dict)
    read_input(set1_file, simple_species)
    giveAntijuncFlag(is_specific)
    exon_pair_nooutput()
    search_near()
    giveBhbFlag_250827()  # After search_near() and exon_pair_nooutput()
    outputRc2(output_file_name)


