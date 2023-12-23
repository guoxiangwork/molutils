from rdkit import Chem
from wanhuamolml.descriptors2D.molecular_chain_length import no_br_length


with open('1497.txt','r') as f:
    line = f.read()
    smi_str = line.split(', ')
    smi = [x.strip("'") for x in smi_str]

cansmi = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in smi]


# 去重
newcansmi = []
for smi in cansmi:
    if smi not in newcansmi:
        newcansmi.append(smi)

print('正则化smiles后剩', len(newcansmi))

cansmi_no = []
for smi in newcansmi:
    if no_br_length(smi):
        cansmi_no.append(smi)
oldsmi = open('remove_duplication-562.smi','r').read().splitlines()

print('去除N或O后剩', len(cansmi_no))

fl = []
for cs in cansmi_no:
    if cs not in oldsmi:
        #print(cs)
        fl.append(cs)
print('去除上次的后剩', len(cansmi_no))
fl1 = []
for smi in fl:
    if 'c1' not in smi:
        fl1.append(smi)
print('去重含有苯环后剩',len(fl1))

print('\n'.join(fl1),file=open('24.smi','w'))