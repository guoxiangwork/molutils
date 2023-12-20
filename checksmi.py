smiles_list = open('remove_duplication-562.smi','r').read().split()
i=1
for idx,smi in enumerate(smiles_list):
    if '[C]' in smi:
        print(smi,idx)
        i+=1
print(i)