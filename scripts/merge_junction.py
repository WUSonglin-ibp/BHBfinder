import collections
import argparse
import re
import numpy as np
import os

n_info = 4
def readCircs(files):
    # read info and save uniq circ info
    sample_num = len(files)
    sample_names = []
    circ_info = {}

    for i, file in enumerate(files):
        sample_names.append(re.search('(.+)_junction', file).group(1))
        name_count = 0
        for l in open(file):
            l = l.rstrip('\n')
            if not l:
                continue
            if l[0] == '#':
                continue
            a = l.split('\t')
            # key: chr, start, end, strand, type
            k = (a[0], int(a[1]), int(a[2]), a[5], a[6])
            if k not in circ_info.keys():
                name_count += 1
                circ_info[k] = {}
                circ_info[k]['total_nodup'] = 0
                circ_info[k]['nodup_flag'] = False
                circ_info[k]['sample_info'] = [0] * sample_num * n_info  # 'Nodup_count', 'Read_count', 'Normal', 'Spliced_rate', 'Signal_rate', 'Abundance'
                circ_info[k]['output_line'] = '\t'.join( a[5:8] + a[9:19] + a[21:24] )
            circ_info[k]['total_nodup'] += int(a[4])
            circ_info[k]['sample_info'][i*n_info : (i+1)*n_info] = a[4], a[8], a[19], a[20]
    ## sort circ info according to chr and start
    sorted_circ_info = collections.OrderedDict()
    def custom_sort(item):
        return item[0][:2]
    for k, v in sorted(circ_info.items(), key=custom_sort):
        sorted_circ_info[k] = v
    return sorted_circ_info, sample_num, sample_names

def outputMerge(output_file, samples):
#add standard deviation of splicing rate
    fh = open(output_file, 'w')
    header = '\t'.join(
        ['#Chr', 'Start', "End", "Name", 'Total_nodup_count', "Strandness", "Type", "Length", "maxscore", "Overlap",
         "Splitsize", "Helix1_count", "Helix2_count", "Helix3_count", "Seq1", "Pairing", "Seq2", 'Intron_seq',
         "Gene_ID", "Gene_type", "Gene_product", "Ave_spliced_rate", "SD_spliced_rate"])
    header_basenames = ['Nodup_count', 'Read_count', 'Normal', 'Spliced_rate']
    for sample in samples:
        for base in header_basenames:
            header += '\t' + base + '/' + os.path.basename(sample)
    fh.write(f'{header}\n')

    name_count = 1
    for k, v in circ_info.items():
        # cal ave splicedrate
        spliced_rate_new = 0
        spliced_rate_sd = 0
        count = 0
        spliced_rate_list = []
        for i in range(sample_num):
            #Bug: spliced_rate = 0.0000 for deduplicated hybrid read =1 and normal read = 0 are included for averaging. 
            #if int(v['sample_info'][i*n_info]) > 0: 
            if int(v['sample_info'][i*n_info]) > 1: # deduplicated hybrid read > 1, same as selection in junction.bed
                count += 1
                #spliced_rate_new += float(v['sample_info'][n_info * i + 3])
                spliced_rate_list.append(float(v['sample_info'][n_info * i + 3]))
        if count> 0:
            #spliced_rate_new = spliced_rate_new/count
            spliced_rate_new = np.average(spliced_rate_list)
            spliced_rate_sd = np.std(spliced_rate_list)

        # output
        fh.write('\t'.join(map(str, [k[0], k[1], k[2], 'junc_'+str(name_count),  v['total_nodup']] )))
        fh.write('\t' + v['output_line'] +  '\t' + str(round(spliced_rate_new, 4)) + '\t' +  str(round(spliced_rate_sd, 4)) )
        fh.write('\t' + '\t'.join(map(str, v['sample_info'])) + '\n')
        name_count += 1
    fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Circular RNA ")
    parser.add_argument('--beds', type=str, nargs='+', required=True, help="Bed files for merge")
    parser.add_argument('-o', help="Output file", default='merged_junction.bed')

    args = parser.parse_args()

    circ_info, sample_num, sample_names = readCircs(args.beds)
    outputMerge(args.o, sample_names)

