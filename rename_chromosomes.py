import sys

# Сопоставление RefSeq -> Chromosome
refseq_to_chr = {
    "NC_037328.1": "1", "NC_037329.1": "2", "NC_037330.1": "3", "NC_037331.1": "4", "NC_037332.1": "5", 
    "NC_037333.1": "6", "NC_037334.1": "7", "NC_037335.1": "8", "NC_037336.1": "9", "NC_037337.1": "10", 
    "NC_037338.1": "11", "NC_037339.1": "12", "NC_037340.1": "13", "NC_037341.1": "14", "NC_037342.1": "15", 
    "NC_037343.1": "16", "NC_037344.1": "17", "NC_037345.1": "18", "NC_037346.1": "19", "NC_037347.1": "20", 
    "NC_037348.1": "21", "NC_037349.1": "22", "NC_037350.1": "23", "NC_037351.1": "24", "NC_037352.1": "25", 
    "NC_037353.1": "26", "NC_037354.1": "27", "NC_037355.1": "28", "NC_037356.1": "29", "NC_037357.1": "X", 
    "NC_082638.1": "Y", "NC_006853.1": "MT"
}

# Читаем входные и выходные файлы из аргументов командной строки
input_bed = sys.argv[1]
output_bed = sys.argv[2]

with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue
        
        fields = line.strip().split("\t")
        if fields[0] in refseq_to_chr:
            fields[0] = refseq_to_chr[fields[0]]
        
        outfile.write("\t".join(fields) + "\n")

print(f"Файл {output_bed} создан с переименованными хромосомами.")
