import re
import click
import os
import pandas as pd
import numpy as np
import requests
import json

#@click.command()
#@click.option('--path',help='数据库的整体路径')

def kegg(path):
    path = os.path.abspath(path)
    if not os.path.exists(os.path.join(path, 'KEGG')):
        os.mkdir(os.path.join(path,'KEGG'))
    path = os.path.join(path,'KEGG')
    if not os.path.exists(os.path.join(path, 'genes')):
        os.mkdir(os.path.join(path, 'genes'))
    if not os.path.exists(os.path.join(path, 'genes','Animals')):
        os.mkdir(os.path.join(path, 'genes','Animals'))
    if not os.path.exists(os.path.join(path, 'genes','Plants')):
        os.mkdir(os.path.join(path, 'genes','Plants'))
    if not os.path.exists(os.path.join(path, 'genes','Micro')):
        os.mkdir(os.path.join(path, 'genes','Micro'))
    if not os.path.exists(os.path.join(path, 'genes','Fungi')):
        os.mkdir(os.path.join(path, 'genes','Fungi'))
    Animals_fasta = [i for i in os.listdir(os.path.join(path, 'genes','Animals')) if i.endswith('.fasta')]
    Plants_fasta = [i for i in os.listdir(os.path.join(path, 'genes','Plants')) if i.endswith('.fasta')]
    Micro_fasta = [i for i in os.listdir(os.path.join(path, 'genes','Micro')) if i.endswith('.fasta')]
    Fungi_fasta = [i for i in os.listdir(os.path.join(path, 'genes','Fungi')) if i.endswith('.fasta')]
    def cat_fa(ls, name): ##传入列表

        for i in ls:
                with open(os.path.join(path,'genes',name,i),'r') as f, open(os.path.join(path,name,'.fa'),'a') as out:
                    for line in f:
                        out.write(line+'\n')
    cat_fa(Animals_fasta, 'Animals')
    cat_fa(Plants_fasta, 'Plants')
    cat_fa(Micro_fasta, 'Micro')
    cat_fa(Fungi_fasta, 'Fungi')

    ls = ['Animals.fa', 'Plants.fa', 'Micro.fa', 'Fungi.fa']
    with open(os.path.join(path,'kegg,fa'),'w') as out:
        out.write('id KO\n')
        for i in ls:
            with open(os.path.join(path,i),'r') as f:
                for line in f:
                    if line.startswith('>'):
                        line = line.replace('>','').split(None)
                        out.write(line[0]+' '+line[1]+'\n')
    # url = "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir="
    # response = requests.get(url)
    # data = response.json()
    # with open(os.path.join("ko00001.json"), "w") as f:
    #     json.dump(data, f)

    cmd = 'cd %s' %path
    cmd += '&& curl "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=" -o ./ko00001.json'
    cmd += '&& /home/sujieyi/anaconda3/envs/python3.9/bin/python3.9 /home/sujieyi/script/database/kegg/KEGG_anno.py'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in Fungi.fa -d FUNGI'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in Animals.fa -d ANIMAL'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in Plants.fa -d PLANT'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in Micro.fa -d MICRO'
    subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    kegg()
