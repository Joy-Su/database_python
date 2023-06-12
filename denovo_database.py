#!/home/sujieyi/anaconda3/envs/python3.9/bin/python3.9

import click
import os
import subprocess
import nr
import swissprot
import go
import pfam
import kegg
import nt
import requests
from bs4 import BeautifulSoup

@click.command()
@click.option('--path', help = '数据库的整体路径')
@click.option('--d', help = '需要整理的库nr,swissprot,pfam,go,kegg,nt')
@click.option('--check', type =click.Choice(['T','F']), help = 'T or F ; 确认是否进行更新检索;如果T：检查是否为最新版本；如果F：不检查；   default:F', default = 'F')


def denovo(path, d, check):
    path = os.path.abspath(path)
    d = d.lower().split(',')
    if str(check) == 'F':
        if 'nr' in d:
            print ('??')
            nr.NR(path)
        if 'swissprot' in d:
            swissprot.Swissprot(path)
        if 'pfam' in d:
            pfam.pfam(path)
        if 'go' in d:
            go.go(path)
        if 'kegg' in d:
            kegg.kegg(path)
        if 'nt' in d:
            nt.NT(path)
    else:
        if 'nr' in d:
            url = 'https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/'
            label = 'nr.gz'
            txt = '/public/mid/rna/denovo_database_2022/New/nr.txt'
            check_update(url, label, txt, nr.NR, 'nr', path)
        if 'swissprot' in d:
            url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/'
            label = 'uniprot_sprot.dat.gz'
            txt = '/public/mid/rna/denovo_database_2022/New/swissprot.txt'
            check_update(url, label, txt, swissprot.Swissprot, 'swissprot', path)
        if 'pfam' in d:
            url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/'
            label = 'Pfam-A.hmm.gz'
            txt = '/public/mid/rna/denovo_database_2022/New/pfam.txt'
            check_update(url, label, txt, pfam.pfam, 'pfam', path)
        if 'go' in d:
            url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/'
            label = 'uniprot_sprot.dat.gz'
            txt = '/public/mid/rna/denovo_database_2022/New/swissprot.txt'
            check_update(url, label, txt, go.go, 'swissprot', path)
        if 'kegg' in d:
            kegg.kegg(path)
        if 'nt' in d:
            url = 'https://ftp.ncbi.nlm.nih.gov/blast/db/'
            label = 'nt.00.tar.gz'
            txt = '/public/mid/rna/denovo_database_2022/New/nt.txt'  #这个记录nr.00.tar.gz更新的时间即可
            check_update(url, label, txt, nt.NT, 'nt', path)


def check_update(url, label, txt, database, DATABASE, path_):
    r = requests.get(url)
    html_content = r.text
    soup = BeautifulSoup(html_content, 'html.parser')
    if DATABASE in ['nr','swissprot','nt']:
        url_time = '-'.join(soup.find('a', href=label).next_sibling.strip().split(None)[0:2])
    if DATABASE in ['pfam']:
        url_time = '-'.join(soup.find('a',href=label).find_next('td', align='right').text.strip().split(None)[0:2])
    ls_time = []
    with open(txt,'r') as f:
        for line in f:
            line = line.strip('\n')
            ls_time.append(line)
    if url_time in ls_time:
        print ('nothing to be updated!')
        pass
    else:
        print('检测到新版本！开始更新-----')
        database(path_)
if __name__ == '__main__':
    denovo()