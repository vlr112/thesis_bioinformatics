import sys

external_file = sys.argv[1] #"external_af.txt"
vcf_file = sys.argv[2]  #"input.vcf"
output_file = sys.argv[3] #"updated.vcf"

# Read the external file and store the AF values in a dictionary
af_dict = {}
with open(external_file, "r") as file:
    for line in file:
        variant_id, af = line.strip().split("\t")
        af_dict[variant_id] = af

# Update the AF values in the VCF file
with open(vcf_file, "r") as input_vcf, open(output_file, "w") as output_vcf:
    for line in input_vcf:
        if line.startswith("#"):
            # Write header lines as they are
            output_vcf.write(line)
        else:
            fields = line.strip().split("\t")
            position = fields[1]  # Assuming the position is in the second column
            # if position == '5005':
            #     break
            if position in af_dict:
                af = float(af_dict[position])
                if af > 0.5:
                    af = 1 - af
                info_field = fields[7]
                # print('info_field', info_field)
                updated_info = info_field.replace(info_field, "AF1=" + str(af))
                # print(updated_info)
                fields[7] = updated_info
            output_vcf.write("\t".join(fields) + "\n")
