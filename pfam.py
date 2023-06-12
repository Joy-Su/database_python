import os
import subprocess
import click
import requests
import gzip
import pandas as pd

#@click.command()
#@click.option('--path',help='数据库的整体路径')

def pfam():
    path = os.path.join(path)
    if not os.path.exists(os.path.join(path, 'pfam')):
        os.makedirs(os.path.join(path, 'pfam'))
    path = os.path.join(path,'pfam')
    r = requests.get('https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz')
    with open(os.path.join(path, 'Pfam-A.hmm.gz'), 'wb') as file:
        file.write(r.content)
    g_file = gzip.GzipFile(os.path.join(path,'Pfam-A.hmm.gz'))
    open(os.path.join(path, 'Pfam-A.hmm'), 'wb+').write(g_file.read())
    g_file.close()
    os.remove(os.path.join(path, 'Pfam-A.hmm.gz'))
    cmd = '/home/fanyucai/software/hmmer/hmmer-v3.1b2/bin/hmmpress %s' % (os.path.join(path, 'Pfam-A.hmm'))
    subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    pfam()