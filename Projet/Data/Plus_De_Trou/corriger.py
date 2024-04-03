
with open('Plus_De_Trou.txt', 'r') as input_file:
    with open('Plus_De_Trou_m.txt', 'w') as output_file:
        for line in input_file:

            parts = line.split(";")

            parts.insert(9, "[")
            parts.insert(11, "[")

            modified_string = ";".join(parts)

            output_file.write(modified_string)

