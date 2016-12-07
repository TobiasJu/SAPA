import csv
import sys

print "converting file to SAPA format"

counter = 0
clinvar = open('converted_hg38_pathogenic.csv', 'w')
clinvar.write("Chr,Pos,Ref,Alt,Clinical significance,Gene\n")
with open("clinvar_pathogenic_snp.txt") as csvfile:
    variant_lines = csv.reader(csvfile, delimiter='\t', quotechar='"')
    # skip header if there
    has_header = csv.Sniffer().has_header(csvfile.read(100))
    csvfile.seek(0)  # rewind
    incsv = csv.reader(csvfile)
    if has_header:
        header = next(variant_lines)
        print header
    for variant_line in variant_lines:
        export_string = ""
        ref_alt = variant_line[0].split(">")
        if len(ref_alt) == 2:
            ref = ref_alt[0][-1]
            alt = ref_alt[1][0]
            chr = variant_line[6]
            pos = variant_line[7]
            clin_sig = variant_line[4].split("(")
            buildver = variant_line[8]
            gene = variant_line[1]

            if buildver == "GRCh38" and chr and pos:  # and clin_sig == "Pathogenic"
                counter += 1
                print chr
                print pos
                print ref
                print alt
                print clin_sig[0]
                print gene
                # print buildver
                export_string += chr+","+pos+","+ref+","+alt+","+clin_sig[0]+","+gene+"\n"
                clinvar.write(export_string)

        if counter == 100:
            break
clinvar.close()
print "END"