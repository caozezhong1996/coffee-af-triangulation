# -*- coding: utf-8 -*-
"""Search BBJ raw GWAS by genomic coordinate for the 6 missing variants."""
import zipfile, gzip

targets = {
    ('15', '75027880'): 'rs2472297',
    ('19', '41353107'): 'rs56113850',
    ('22', '24870527'): 'rs17842490',
    ('8',  '33790200'): 'rs78267637',
    ('10', '135315795'): 'rs117810762',
    ('15', '75174251'): 'rs117968677',
}
# also count nearby variants (same chr, +/-5kb) to see if region is covered at all
nearby = {v: 0 for v in targets.values()}
pos_of = {v: (int(k[0]), int(k[1])) for k, v in targets.items()}

z = zipfile.ZipFile(r'coloc\bbj_AF.zip')
name = [n for n in z.namelist() if n.endswith('.txt.gz')][0]
found = {}
n_total = 0
with z.open(name) as f:
    with gzip.open(f, 'rt') as g:
        hdr = g.readline().strip().split('\t')
        # columns: v CHR POS SNID Allele1 Allele2 AC_Allele2 AF_Allele2 imputationInfo N BETA SE Tstat p.value
        for line in g:
            n_total += 1
            p = line.split('\t')
            try:
                pos_val = int(float(p[2]))
            except ValueError:
                continue
            key = (p[1], str(pos_val))
            if key in targets:
                found[targets[key]] = dict(zip(['CHR','POS','SNID','A1','A2','AF_A2','imp','N','BETA','SE','p'],
                                               [p[1],p[2],p[3],p[4],p[5],p[7],p[8],p[9],p[10],p[11],p[13].strip()]))
            else:
                chr_ = p[1]
                for v,(c,pp) in pos_of.items():
                    if chr_ == str(c) and abs(pos_val-pp) <= 5000:
                        nearby[v] += 1

print('total variant lines scanned:', n_total)
print()
for v in targets.values():
    if v in found:
        print(v, 'FOUND:', found[v])
    else:
        print(v, 'MISSING at exact coordinate; variants within +/-5kb:', nearby[v])
