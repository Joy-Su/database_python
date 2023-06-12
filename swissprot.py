import os
import subprocess
import click
import requests
import gzip
import pandas as pd

#@click.command()
#@click.option('--path', help = '数据库的整体路径')

def Swissprot(path):
    path = os.path.abspath(path)
    if not os.path.exists(os.path.join(path,'Swissprot')):
        os.mkdir(os.path.join(path,'Swissprot'))
    path = os.path.join(path,'Swissprot')
    '''download fasta'''
    r = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz')
    with open(os.path.join(path,'uniprot_sprot.fasta.gz'),'wb') as f:
        f.write(r.content)
    g_file = gzip.GzipFile(os.path.join(path,'uniprot_sprot.fasta.gz'))
    open(os.path.join(path,'uniprot_sprot.fasta'),'wb+').write(g_file.read())
    g_file.close()
    os.remove(os.path.join(path,'uniprot_sprot.fasta.gz'))
    # '''download uniprot_dat'''
    # r = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz')
    # with open(os.path.join(path,'uniprot_sprot.dat.gz'),'wb') as f:
    #     f.write(r.content)
    # g_file = gzip.GzipFile(os.path.join(path,'uniprot_sprot.dat.gz'))
    # open(os.path.join(path,'uniprot_sprot.dat'),'wb+').write(g_file.read())
    # g_file.close()
    # os.remove(os.path.join(path,'uniprot_sprot.dat.gz'))

    '''下载各物种'''
    #'&& echo -e "archaea\nbacteria\nviruses\nfungi\nplants\nhuman\ninvertebrates\nmammals\nrodents\nvertebrates" > uni_spe.txt'
    uni_spe = ['archaea','bacteria','viruses','fungi','plants','human','invertebrates','mammals','rodents','vertebrates']
    '''&& cat uni_spe.txt|while read a; \
            do /home/fanyucai/software/axel/axel-v2.15/bin/axel -n 20 \
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_"$a".dat.gz"; \
            gunzip "uniprot_sprot_$a.dat.gz"; \
            perl script/deal_uniprot_sprot.pl "uniprot_sprot_"$a".dat" "uniprot_sprot_"$a".dat-format.txt"; done'''

    def download_gz(file):
        r = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_'+file+'.dat.gz')
        with open(os.path.join(path,'uniprot_sprot_'+file+'.dat.gz'),'wb') as f:
            f.write(r.content)
        g_file = gzip.GzipFile(os.path.join(path,'uniprot_sprot_'+file+'.dat.gz'))
        open(os.path.join(path,'uniprot_sprot_'+file+'.dat'),'wb+').write(g_file.read())
        g_file.close()
        os.remove(os.path.join(path,'uniprot_sprot_'+file+'.dat.gz'))
        hash_dict = {}
        with open(os.path.join(path, 'uniprot_sprot_' + file + '.dat'), 'r') as infile:
            for line in infile:
                line = line.strip()
                line = line.strip('\n')
                line = line.split(None, 1)
                if line[0] == 'ID' and len(line)>1:
                    ID = line[1].split(None)[0]
                    hash_dict[ID] = []
                    i = 0
                    j = 0
                    k = 0
                if line[0] == 'AC' and len(line) > 1:
                    AC = line[1].split(';')[0]
                    hash_dict[ID].append(AC)
                if line[0] == 'DE' and len(line)>1:
                    if 'RecName' in line[1]:
                        if k == 0:
                            hash_dict[ID].append(line[1])
                            k += 1
                if line[0] == 'GN' and len(line) > 1:
                    # GN=line[1]
                    if j == 0:
                        hash_dict[ID].append(line[1])
                        j += 1
                    else:
                        hash_dict[ID][-1] += ' '
                        hash_dict[ID][-1] += line[1]
                if line[0] == 'OS' and len(line) > 1:
                    if i == 0:
                        hash_dict[ID].append(line[1])
                        i += 1
                    else:
                        hash_dict[ID][-1] += ' '
                        hash_dict[ID][-1] += line[1]
        format = pd.DataFrame.from_dict(hash_dict, orient='index').reset_index()
        format.to_csv(os.path.join(path,'uniprot_sprot_' + file + '.dat-format.txt'),sep='\t',index=False,header=None)
    list(map(download_gz,uni_spe))
    #'&& cat uniprot_sprot_archaea.dat-format.txt uniprot_sprot_bacteria.dat-format.txt uniprot_sprot_viruses.dat-format.txt|awk -F"\t" \'{print "sp|"$2"|"$1}\'> MICRO.txt'
    with open(os.path.join(path,'uniprot_sprot_archaea.dat-format.txt'),'r') as f1, open(os.path.join(path,'uniprot_sprot_bacteria.dat-format.txt'),'r') as f2, \
            open(os.path.join(path, 'uniprot_sprot_viruses.dat-format.txt'),'r') as f3, open(os.path.join(path, 'MICRO.txt'),'w') as outfile:
        for line in f1:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')
        for line in f2:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')
        for line in f3:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')

    '''&& cat uniprot_sprot_human.dat-format.txt uniprot_sprot_invertebrates.dat-format.txt uniprot_sprot_mammals.dat-format.txt uniprot_sprot_rodents.dat-format.txt 
    uniprot_sprot_vertebrates.dat-format.txt|awk -F"\t" \'{print "sp|"$2"|"$1}\'>ANIMALS.txt'''
    with open(os.path.join(path, 'uniprot_sprot_human.dat-format.txt'),'r') as f1, open(os.path.join(path, 'uniprot_sprot_invertebrates.dat-format.txt'),'r') as f2, \
        open(os.path.join(path, 'uniprot_sprot_mammals.dat-format.txt'),'r') as f3, open(os.path.join(path, 'uniprot_sprot_rodents.dat-format.txt'),'r') as f4, \
        open(os.path.join(path, 'uniprot_sprot_vertebrates.dat-format.txt'),'r') as f5, open(os.path.join(path, 'ANIMALS.txt'),'w') as outfile:
        for line in f1:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')
        for line in f2:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')
        for line in f3:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')
        for line in f4:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')
        for line in f5:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')

    #'&& awk -F"\t" \'{print "sp|"$2"|"$1}\' uniprot_sprot_fungi.dat-format.txt>FUNGI.txt'
    with open(os.path.join(path, 'uniprot_sprot_fungi.dat-format.txt'),'r') as f, open(os.path.join(path, 'FUNGI.txt'),'w') as outfile:
        for line in f:
            line = line.strip().strip('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')

    #'&& awk -F"\t" \'{print "sp|"$2"|"$1}\' uniprot_sprot_plants.dat-format.txt>PLANTS.txt'
    with open(os.path.join(path, 'uniprot_sprot_plants.dat-format.txt'),'r') as f, open(os.path.join(path, 'PLANTS.txt'),'w') as outfile:
        for line in f:
            line = line.strip().split('\t')
            outfile.write('sp|'+line[1]+'|'+line[0]+'\n')
    cmd = 'cd %s' %(path)
    cmd += '&& /home/dna/software/Bin/seqtk subseq uniprot_sprot.fasta MICRO.txt>MICRO.fa && /home/dna/software/Bin/seqtk subseq uniprot_sprot.fasta PLANTS.txt>PLANTS.fa && /home/dna/software/Bin/seqtk subseq uniprot_sprot.fasta ANIMALS.txt>ANIMALS.fa && /home/dna/software/Bin/seqtk subseq uniprot_sprot.fasta FUNGI.txt>FUNGI.fa'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in ANIMALS.fa -d ANIMALS && /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in MICRO.fa -d MICRO && /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in PLANTS.fa -d PLANTS && /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in FUNGI.fa -d FUNGI'
    subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    Swissprot()
