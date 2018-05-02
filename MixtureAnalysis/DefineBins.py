chr_19_size = 61431566

bin = 0

bins = []

output = open("bins.txt", 'w')

while bin <= chr_19_size:
    bin = bin + 100
    bins.append(bin)
    # output.write(str(bin))
    # output.write("\n")
    #


print("Done Generating... writing...")

with open("chr19_bins.txt", 'w') as f:
    for bin in bins:
        f.write(str(bin))
        f.write("\n")

print("Done")
