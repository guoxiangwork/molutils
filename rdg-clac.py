from openbabel import pybel
import subprocess
from pathlib import Path


def cc(wd):
    mol = next(pybel.readfile(format='xyz',filename=wd+'/cluster.xyz'))
    mol.write(format='xyz',filename=wd+'/molden.xyz',overwrite=True)

    subprocess.run(['xtb','molden.xyz','--molden','--chrg','1','--uhf','3'],cwd=wd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    input='\n'.join(['2','1','2','2,3','0','7','1','-10','q'])
    command=['multiwfn','molden.input']
    process=subprocess.Popen(command,cwd=wd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output,error=process.communicate(input=input.encode())
    output=output.decode().strip().split('\n')
    for line in output:
        if line.startswith(' Reduced density gradient with promolecular approximation:'):
            rdg_xtb = float(line.split()[-1])
            rdg_rev = (rdg_xtb-0.1497)/0.7321
            print(rdg_rev)
            return rdg_rev

for wd in Path('.').iterdir():
    print(str(wd))
    cc(str(wd))