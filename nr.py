import subprocess
import click
import requests
import zipfile
import os
import shutil
import gzip
import pandas as pd
import collections

#@click.command()
#@click.option('--path', help = '数据库的整体路径')
def NR(path):
    path = os.path.abspath(path)
    if not os.path.exists(os.path.join(path,'NR')):
        os.makedirs(os.path.join(path,'NR'))
    path = os.path.join(path, 'NR')
    """ download taxdmp """
    r = requests.get('http://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip')
    with open(os.path.join(path,'taxdmp.zip'), 'wb') as f:
        f.write(r.content)
    zf = zipfile.ZipFile(os.path.join(path,'taxdmp.zip'))
    zf.extractall()
    os.remove(os.path.join(path, 'taxdmp.zip'))
    """ download nr """
    r = requests.get('https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz')
    with open(os.path.join(path, 'nr.gz'), 'wb') as f:
        f.write(r.content)
    g_file = gzip.GzipFile(os.path.join(path, 'nr.gz'))
    open (os.path.join(path, 'nr'), 'wb+').write(g_file.read())
    g_file.close()
    os.remove(os.path.join(path, 'nr.gz'))

    shutil.copyfile('/home/sujieyi/script/database/nr/fungi_taxid.txt', os.path.join(path, 'fungi_taxid.txt'))

    with open(os.path.join(path, 'nodes.dmp'), 'r') as infile, open(os.path.join(path, 'nodes.txt'), 'w') as outfile:
        for line in infile:
            line = line.replace('\t|','')
            columns = line.split('\t')
            outfile.write(columns[0] + '\t' + columns[4] + '\n')

    fungi_id = pd.read_csv(path+'/fungi_taxid.txt',sep='\t',header=None,names=['id'])
    nodes_txt = pd.read_csv(path+'/nodes.txt',sep='\t',header=None,names=['id','type'])
    match_fungi = {}
    nodes_txt['id'].apply(lambda x: match_fungi.update({x: 12}) if x in fungi_id['id'].tolist() else match_fungi.update({x: nodes_txt[nodes_txt['id'] == x]['type'].values[0]}))
    match_fungi = pd.DataFrame.from_dict(match_fungi, orient='index',columns=['type']).reset_index().rename(columns={'index':'id'})
    match_fungi.to_csv(os.path.join(path, 'nodes_new.txt'),sep='\t',index=False,header=False)
    with open(os.path.join(path, 'division.dmp'),'r') as f, open(os.path.join(path, 'division.txt'),'w') as out:
        for line in f:
            line = line.replace('\t|','')
            column = line.split('\t')
            out.write(column[0]+'\t'+column[1]+'\n')
        out.write('12'+'\t'+'FUN')
    division = pd.read_csv(os.path.join(path, 'division.txt'),sep='\t',header=None,names=['type','name'])
    division2nodes = pd.merge(match_fungi,division,on='type')[['id','type','name']]
    division2nodes.to_csv(os.path.join(path, 'division2nodes.txt'),sep='\t',index=False,header=False)
    #'&& grep \">\" nr|sed \"s/ /\t/\"|sed \"s/ \[/\t/g\" |sed \"s/]//g\"|sed \"s/>//\" > all-anno.txt'
    with open(os.path.join(path, 'nr')) as f, open(os.path.join(path, 'all-anno.txt'),'w') as out:
        for line in f:
            if line.startswith('>'):
                line = line.replace(' ','\t',1).replace(' [','\t').replace(']','').replace('>','').replace('\n','')
                out.write(line+'\n')
    #'&& sed \"1,2d\" names.dmp|cut -f1,3,7 | awk -F\"\t\" -v OFS=\"\t\" "\$3~/scientific name/{print \$1,\$2}" > taxonomyID_species.txt'
    # '&& sed \'s/"//g\' taxonomyID_species.txt > taxonomyID_species-new.txt'
    with open(os.path.join(path,'names.dmp'),'r') as f, open(os.path.join(path,'taxonomyID_species-new.txt'),'w') as out:
        lines = f.readlines()[2:]
        for line in lines:
            line = line.replace('"','')
            line = line.split('\t')
            if line[6] == 'scientific name':
                out.write(line[0]+'\t'+line[2]+'\n')
    #'&& awk -F\"\t\" -v OFS=\"\t\" \"NR==FNR{a[\$1]=\$1;b[\$1]=\$3;next}{if(\$1==a[\$1])print  \$2,b[\$1]}\"  division2nodes.txt taxonomyID_species-new.txt|sort -u>name2id.txt'
    taxonomy = pd.read_csv(os.path.join(path, 'taxonomyID_species-new.txt'),sep='\t',header=None,names=['id','tax'])
    name2id = pd.merge(taxonomy,division2nodes[['id','name']],how='inner',on='id')[['tax','name']]
    name2id.to_csv(os.path.join(path,'name2id.txt'),sep='\t',index=False,header=False)

    #'&& awk -F\"\t\" -v OFS=\"\t\" \"NR==FNR{a[\$1]=\$1;b[\$1]=\$2;next}{if(\$3==a[\$3])print \$1,b[\$3]}\" name2id.txt all-anno.txt > final_anno.txt'
    all_anno = pd.read_csv(os.path.join(path,'all-anno.txt'),sep='\t',header=None,error_bad_lines=False,names=['id','description','tax'])
    final_anno = pd.merge(all_anno,name2id,on='tax')[['id','name']]
    final_anno.to_csv(os.path.join(path,'final_anno.txt'),sep='\t',index=False,header=False)

    #'&& echo -e \"BCT\nPHG\nVRL\nENV\nUNA\nSYN\" >MICRO && echo -e \"PLN\nUNA\nSYN\" >PLANT && echo -e \"MAM\nINV\nPRI\nROD\nVRT\nSYN\nUNA\" >ANIMAL && echo -e \"FUN\nSYN\nUNA\" > FUNGI'
    MICRO = ['BCT','PHG','VRL','ENV','UNA','SYN']
    PLANT = ['PLN','UNA','SYN']
    ANIMAL = ['MAM','INV','PRI','ROD','VRT','SYN','UNA']
    FUNGI = ['FUN','SYN','UNA']

    with open(os.path.join(path,'MICRO'),'w') as f:
        for i in MICRO:
            f.write(i+'\n')
    with open(os.path.join(path,'PLANT'),'w') as f:
        for i in PLANT:
            f.write(i+'\n')
    with open(os.path.join(path,'ANIMAL'),'w') as f:
        for i in ANIMAL:
            f.write(i+'\n')
    with open(os.path.join(path,'FUNGI'),'w') as f:
        for i in FUNGI:
            f.write(i+'\n')

    #'&& awk -F"\t" -v OFS="\t" "NR==FNR{a[\$1]=\$1}NR>FNR{if(a[\$2]==\$2)print \$1}" FUNGI final_anno.txt > FUNGI.gilist'
    MICRO_gilist = final_anno[final_anno['name'].isin(MICRO)]['id']
    MICRO_gilist.to_csv(os.path.join(path,'MICRO.gilist'),sep='\t',index=False,header=False)
    PLANT_gilist = final_anno[final_anno['name'].isin(PLANT)]['id']
    PLANT_gilist.to_csv(os.path.join(path,'PLANT.gilist'),sep='\t',index=False,header=False)
    ANIMAL_gilist = final_anno[final_anno['name'].isin(ANIMAL)]['id']
    ANIMAL_gilist.to_csv(os.path.join(path,'ANIMAL.gilist'),sep='\t',index=False,header=False)
    FUNGI_gilist = final_anno[final_anno['name'].isin(FUNGI)]['id']
    FUNGI_gilist.to_csv(os.path.join(path,'FUNGI.gilist'),sep='\t',index=False,header=False)
    cmd = 'cd %s' %(path)
    cmd += '&& /home/dna/software/Bin/seqtk subseq -l 60 nr FUNGI.gilist > FUNGI.fa'
    cmd += '&& /home/dna/software/Bin/seqtk subseq -l 60 nr ANIMAL.gilist > ANIMAL.fa'
    cmd += '&& /home/dna/software/Bin/seqtk subseq -l 60 nr PLANT.gilist > PLANT.fa'
    cmd += '&& /home/dna/software/Bin/seqtk subseq -l 60 nr MICRO.gilist > MICRO.fa'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in FUNGI.fa -d FUNGI'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in ANIMAL.fa -d ANIMAL'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in PLANT.fa -d PLANT'
    cmd += '&& /home/rna/softwares/diamond/diamond-2.0.13/bin/diamond makedb --in MICRO.fa -d MICRO'
    subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    NR()












