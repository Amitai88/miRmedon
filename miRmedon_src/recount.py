import pickle

def recount(e_miRbase_path):
    e_miRbase = pickle.load(open(e_miRbase_path, 'rb'))
    sig_sites = {}
    with open('editing_info.txt', 'r') as file:
        for line in file:
            if not line.startswith('miRNA'):
                miR, pos, _, _, _, _, _, p_value = line.rstrip('\n').split('\t')
                p_value = float(p_value)
                pos = int(pos)-1
                if p_value < 0.05:
                    if miR not in sig_sites.keys():
                        sig_sites[miR] = [pos]
                    else:
                        sig_sites[miR] += [pos]
    counts = {}
    with open('counts.txt', 'r') as count_file:
        for line in count_file:
            haplotype, count = line.rstrip('\n').split('\t')
            counts[haplotype] = float(count)
    crt_counts = {}
    for haplotype, count in counts.items():
        if haplotype.find('_e') != -1:
            miR = haplotype[:haplotype.find('_')]
            put_editing_sites = e_miRbase[miR][haplotype]['editingSites']
            try:
                valid_sites = sig_sites[miR]
            except KeyError:
                valid_sites = []
            sites_intersect = set(put_editing_sites).intersection(set(valid_sites))
            if len(put_editing_sites) == len(sites_intersect):
                crt_counts[haplotype] = count
            else:
                haplotypes = list(e_miRbase[miR].keys())
                editing_sites = [tuple(e_miRbase[miR][haplotype]['editingSites']) for haplotype in haplotypes]
                alt_haplotype = {e_sites: haplotype for e_sites, haplotype in zip(editing_sites, haplotypes)}
                sites_intersect = tuple(sites_intersect)
                try:
                    alt_haplotype = alt_haplotype[sites_intersect]
                    if alt_haplotype not in crt_counts.keys():
                        crt_counts[alt_haplotype] = count
                    else:
                        crt_counts[alt_haplotype] += count
                except KeyError:
                    continue
        else:
            if haplotype not in crt_counts.keys():
                crt_counts[haplotype] = count
            else:
                crt_counts[haplotype] += count
    with open('crt_counts.txt', 'w') as crt_counts_file:
        crt_counts_file.write('miRNA_form' + '\t' + 'count' + '\t' + 'sequence' + '\t' + 'editing_sites' + '\n')
        for haplotype, count in crt_counts.items():
            miR = haplotype[:haplotype.find('_')]
            seq = e_miRbase[miR][haplotype]['sequence']
            put_editing_sites = ', '.join(str(x+1) for x in e_miRbase[miR][haplotype]['editingSites'])
            crt_counts_file.write(haplotype + '\t' + str(count) + '\t' + seq + '\t' + put_editing_sites + '\n')

