import csv

# f = open("../LocusCoverageperIndividual_nSSR_FullLength.txt","r")
# print(f.read())
with open('../LocusCoverageperIndividual_nSSR_FullLength.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ')
    # for row in spamreader:
    #     print(', '.join(row))

