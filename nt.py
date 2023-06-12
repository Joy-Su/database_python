import subprocess
import requests
import tarfile
import os


def NT(path):
    path = os.path.abspath(path)
    if not os.path.exists(os.path.join(path, 'NT')):
        os.makedirs(os.path.join(path, 'NT'))
    path = os.path.join(path, 'NT')
    ls = [str(num).zfill(2) for num in range(0, 100)]
    with open(os.path.join(path,'a.list'),'w') as f:
        for i in ls:
            f.write(str(i)+'\n')
    # ls = [str(num).zfill(2) for num in range(1, 100)]
    # for i in ls:
    #     url = 'https://ftp.ncbi.nlm.nih.gov/blast/db/nt.' + str(i) +'.tar.gz'
    #     r = requests.get(url)
    #     with open(os.path.join(path, 'nr.'+str(i)+'.tar.gz'), 'wb') as f:
    #         f.write(r.content)
    #     t_file = tarfile.open(os.path.join(path, 'nr.'+str(i)+'.tar.gz'))
    #     t_file.extractall(os.path.join(path, 'nr.'+str(i)+'.tar.gz'))
    #     t_file.close()
    #     os.remove('nr.'+str(i)+'.tar.gz')
    #cmd = 'cd %s && cat a.list|while read a; do /home/rna/.aspera/connect/bin/ascp -i /home/rna/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:blast/db/nt.${a}.tar.gz ./;done' % (path)
    cmd = "cd %s && cat a.list|while read a; do /home/fanyucai/software/axel/axel-v2.15/bin/axel -n 20 https://ftp.ncbi.nlm.nih.gov/blast/db/nt.${a}.tar.gz;tar -zxvf nt.${a}.tar.gz;done" % (path)
    subprocess.call(cmd, shell=True)
if __name__ == '__main__':
    NT()