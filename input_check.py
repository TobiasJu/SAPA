# class for checking the input of the annotated SNPs
# convert the abbreviations to full names

def check_input(input_string):
    if input_string == "D":
        return "Deleterious"
    elif input_string == "T":
        return "Tolerated"
    else:
        return input_string