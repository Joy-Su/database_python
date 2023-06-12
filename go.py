import re
import click
import os
import pandas as pd
import numpy as np
import requests
import gzip
#@click.command()
#@click.option('--path',help='数据库的整体路径')

def go(path):
    path = os.path.abspath(path)
    if not os.path.exists(os.path.join(path,'GO')):
        os.makedirs(os.path.join(path,'GO'))
    swiss_path = os.path.join(path, 'Swissprot')
    path = os.path.join(path,'GO')
    '''download uniprot_dat'''
    r = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz')
    with open(os.path.join(path, 'uniprot_sprot.dat.gz'), 'wb') as f:
        f.write(r.content)
    g_file = gzip.GzipFile(os.path.join(path, 'uniprot_sprot.dat.gz'))
    open(os.path.join(path, 'uniprot_sprot.dat'), 'wb+').write(g_file.read())
    g_file.close()
    os.remove(os.path.join(path, 'uniprot_sprot.dat.gz'))
    ID = {}
    AC = {}
    GO_ID = {}
    GO_Cata = {}
    GO_term = {}
    GO_Evid = {}
    num = 1
    with open(os.path.join(path, 'uniprot_sprot.dat'),'r') as f, open(os.path.join(path, 'go2term.pre'),'w') as out1:
        for line in f:
            line = line.strip().replace('\n','')
            if line.startswith('ID   '):
                line = line.split(None)[1]
                if len(line)>1:
                    ID[num] = line
            elif line.startswith('AC   ') and line.endswith(';'):
                line = line.split('   ')[1]
                if AC.get(num,'') == '':
                    AC[num] = line
                else:
                    AC[num] += '#'
                    AC[num] += line
            elif re.match('^DR   GO; (GO:\d{7}); (C|F|P):(.+); ([A-Z]{1,}:.+)\.', line):
                a = re.match('^DR   GO; (GO:\d{7}); (C|F|P):(.+); ([A-Z]{1,}:.+)\.', line)
                out1.write(a.group(1)+'\t'+a.group(2)+'\t'+a.group(3)+'\t'+a.group(4)+'\n')
                if GO_ID.get(num,'') == '':
                    GO_ID[num] = a.group(1)
                else:
                    GO_ID[num] += ','
                    GO_ID[num] += a.group(1)
                if GO_Cata.get(num,'') == '':
                    GO_Cata[num] = a.group(2)
                else:
                    GO_Cata[num] += '|'
                    GO_Cata[num] += a.group(2)
                if GO_term.get(num, '') == '':
                    GO_term[num] = a.group(3)
                else:
                    GO_term[num] += '|'
                    GO_term[num] += a.group(3)
                if GO_Evid.get(num, '') == '':
                    GO_Evid[num] = a.group(4)
                else:
                    GO_Evid[num] += '|'
                    GO_Evid[num] += a.group(4)
            elif line.startswith('//'):
                num += 1
    with open(os.path.join(path, 'all.txt'), 'w') as out2, open(os.path.join(path, 'Uni2GO.pre'), 'w') as out3:
        out2.write('NO.'+'\t'+'ID'+'\t'+'UniPro_ID'+'\t'+'GO_ID'+'\t'+'GO_Category'+'\t'+'GO_Term'+'\t'+'GO_Evidence'+'\n')
        for k in sorted(AC.keys()):
            AC[k] = AC[k].replace('; ',';').replace('#',';')
            out2.write(str(k)+'\t'+ID[k]+'\t'+AC[k].strip(';')+'\t'+GO_ID.get(k,'')+'\t'+GO_Cata.get(k,'')+'\t'+GO_term.get(k,'')+'\t'+GO_Evid.get(k,'')+'\n')
            b = [x for x in AC[k].split(';') if x!='']
            for c in b:
                out3.write(ID[k]+'\t'+c+'\t'+GO_ID.get(k,'')+'\t'+GO_Cata.get(k,'')+'\t'+GO_term.get(k,'')+'\t'+GO_Evid.get(k,'')+'\n')

    go = {}
    ID = {}
    GO_ID = {}
    GO_Cata = {}
    GO_term = {}
    GO_Evid = {}
    with open(os.path.join(path,'Uni2GO.pre'),'r') as f, open(os.path.join(path, 'Uni2GO.txt'),'w')as out:
        out.write('UniPro_ID'+'\t'+'ID'+'\t'+'GO_ID'+'\t'+'GO_Category'+'\t'+'GO_Term'+'\t'+'GO_Evidence'+'\n')
        for a in f:
            a = a.strip().split('\t')
            if len(a) > 2:
                go[a[0]] = a[2]
                if ID.get(a[1], '') == '':
                    ID[a[1]] = a[0]
                    GO_ID[a[1]] = a[2]
                    GO_Cata[a[1]] = a[3]
                    GO_term[a[1]] = a[4]
                    GO_Evid[a[1]] = a[5]
                else:
                    ID[a[1]] += ';' + a[0]
                    GO_ID[a[1]] += ',' + a[2]
                    GO_Cata[a[1]] += '|' + a[3]
                    GO_term[a[1]] += '|' + a[4]
                    GO_Evid[a[1]] += '|' + a[5]
            else:
                go[a[0]] = ''
                ID[a[1]] = a[0]
                GO_ID[a[1]] = ''
                GO_Cata[a[1]] = ''
                GO_term[a[1]] = ''
                GO_Evid[a[1]] = ''
        for i in sorted(ID.keys()):
            out.write(str(i)+'\t'+ID[i]+'\t'+GO_ID[i]+'\t'+GO_Cata[i]+'\t'+GO_term[i]+'\t'+GO_Evid[i]+'\n')

    #&& perl script/deal_uniprot_sprot.pl uniprot_sprot.dat format.txt
    hash_dict = {}
    with open(os.path.join(path, 'uniprot_sprot.dat'),'r') as f:
        index = 0
        for line in f:
            line = line.strip()
            line = line.strip('\n')
            line = line.split(None, 1)
            if line[0] == 'ID' and len(line) > 1:
                index += 1
                if index>1:
                    ###目前发现不存在GN的,对这种其中插入空值
                    print(list(hash_dict)[-1])
                    if len(hash_dict[list(hash_dict)[-1]])==3:
                        hash_dict[list(hash_dict)[-1]].insert(2, ' ')
                ID = line[1].split(None)[0]
                hash_dict[ID] = []
                i = 0
                j = 0
                k = 0
                w = 0
            if line[0] == 'AC' and len(line) > 1:
                AC = line[1].split(';')[0]
                if w == 0:
                    hash_dict[ID].append(AC)
                    w += 1
            if line[0] == 'DE' and len(line) > 1:
                if 'RecName' in line[1]:
                    if k == 0:
                        hash_dict[ID].append(line[1])
                        k += 1
            if line[0] == 'GN' and len(line) > 1:
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
    format.to_csv(os.path.join(path,'format.txt'),sep='\t',index=False,header=None)
    '''&& sed "s/RecName: Full=//;s/Name=//;s/ORFNames=//;s/OrderedLocusNames=//;s/ \t/\t/g;s/ \t/\t/;s/;\t/\t/g;s/. $//" format.txt | 
        sed "1d" |awk -F"\t" -v OFS="\t" "{print \$1,\$2,\$3,\$4,\$5,\$6}" >anno_swi.txt'''
    format[1] = format[1].apply(lambda x: re.sub('RecName: Full=', '', x)).apply(lambda x: re.sub(';', '', x))
    format[2] = format[2].apply(lambda x: re.sub('Name=', '', x)).apply(lambda x: re.sub('ORFNames=', '', x)).apply(lambda x: re.sub('OrderedLocusNames=', '', x)).apply(lambda x: re.sub(';', '', x))
    format[3] = format[3].apply(lambda x:re.sub('.$','',x))
    format.to_csv(os.path.join(path, 'anno_swi.txt'),sep='\t',index=False, header=None)
    #&& grep ">" uniprot_sprot.fasta|awk "{print \$1}"|sed "s/>//" |awk -F"|" -v OFS="\t" "{print \$3,\$0}" >list
    with open(os.path.join(swiss_path, 'uniprot_sprot.fasta'),'r') as f, open(os.path.join(path, 'list'),'w') as out:
        for line in f:
            if line.startswith('>'):
                line = line.split(' ')[0].replace('>','')
                out.write(line.split('|')[2]+'\t'+line+'\n')
    #awk -F"\t" -v OFS="\t" "NR==FNR{a[\$1]=\$1;b[\$1]=\$0;next}{if(\$1==a[\$1]) print b[\$1],\$0}" list anno_swi.txt|awk -F"\t" -v OFS="\t" "{print \$2,\$5,\$6,\$7}" |cut -f1,2 >anno_swissprot.txt
    mylist = pd.read_csv(os.path.join(path,'list'),sep='\t',names=['ID1','ID2'])
    anno_swissprot = pd.merge(mylist, format, left_on='ID1',right_on='index')[['ID2',1]]
    anno_swissprot.to_csv(os.path.join(path, 'anno_swissprot.txt'),sep='\t',index=False,header=None)
    #'&& awk -F"\t" -v OFS="\t" \'NR==FNR{a[$2]=$2;b[$2]=$4}NR>FNR{if(a[$1]==$1)print b[$1],$0}\'  %s/anno_swi.txt Uni2GO.txt >tmp.txt'
    Uni2Go = pd.read_csv(os.path.join(path, 'Uni2GO.txt'),sep='\t')
    tmp = pd.merge(format,Uni2Go,right_on='UniPro_ID',left_on=0,how='inner')[[2,0,'index','GO_ID','GO_Category','GO_Term','GO_Evidence']]
    tmp.to_csv(os.path.join(path, 'tmp.txt'),sep='\t',index=False,header=None)
    #&& awk -F"[|\t]" -v OFS="\t" \'NR==FNR{a[$2]=$2;b[$1+"|"+$2+"|"+$3]=$2}NR>FNR{if(a[$1]==$1)print b[$1],$0}\' %s/anno_swissprot.txt tmp.txt >tmp2.txt' % (Swiss_path)
    #awk -F'[|\t]' -v OFS='' 'NR==FNR{a[$2]=$2;c[$2]=$3;b[$2]=$4}NR>FNR{if(a[$2]==$2)print "sp|",a[$2],"|",c[$2],"\t",b[$2],"\t",$0}' Swissprot/anno_swissprot.txt Swissprot/tmp.txt 改成这个才对
    anno_swissprot_tmp = pd.concat([anno_swissprot, anno_swissprot['ID2'].str.split('|', expand=True)], axis=1)
    anno_swissprot_tmp.columns = ['ID2',0,1,2,3]
    GO_anno = pd.merge(anno_swissprot_tmp,tmp,left_on=3,right_on='index')[['ID2','GO_ID','GO_Term','GO_Category']]
    GO_anno = GO_anno[GO_anno['GO_ID'].apply(lambda x: x is not np.NaN)]
    GO_anno.to_csv(os.path.join(path, 'GO_anno.txt'), sep='\t', index=False,header=['id','GO_ID','GO_Term','GO_Category'])

if __name__ == '__main__':
    go()
