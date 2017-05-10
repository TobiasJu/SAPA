# class for checking the input of the annotated SNPs
# convert the abbreviations to full names

def check_input(input_string):
    if input_string == "D":
        return "Deleterious"
    elif input_string == "T":
        return "Tolerated"
    elif input_string == "P":
        return "Possibly damaging"
    elif input_string == "B":
        return "Benign"
    elif input_string == "D":
        return "damaging"
    elif input_string == "N":
        return "Neutral"
    elif input_string == "U":
        return "Unknown"
    elif input_string == "A":
        return "disease_causing_automatic"
    elif input_string == "N":
        return "polymorphism"
    elif input_string == "P":
        return "polymorphism_automatic"
    elif input_string == "H":
        return "High"
    elif input_string == "M":
        return "Medium"
    elif input_string == "L":
        return "Low"
    else:
        return input_string