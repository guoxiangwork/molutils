with open('bvc6.log','r')  as f:
    for line in f:
        if float(line.split()[0])>0.373:
            print(line.split()[1])