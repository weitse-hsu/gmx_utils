def convert_res_code(res_code):
    """
    Converts an amino acid code between three-letter and one-letter formats.
    Non-standard amino acids or invalid codes are denoted by 'X'.

    Parameters
    ----------
    res_code : str
        The input amino acid code (either three-letter like 'ALA' or one-letter like 'A').

    Returns
    -------
    converted_code : str
        The converted amino acid code (one-letter if input is three-letter, or three-letter if input is one-letter).
        Returns 'X' for invalid or non-standard codes.
    """
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    }

    one_to_three = {v: k for k, v in three_to_one.items()}

    res_code = res_code.upper()
    if len(res_code) == 3:
        converted_code = three_to_one.get(res_code, 'X')
    elif len(res_code) == 1:
        converted_code = one_to_three.get(res_code, 'X')
    else:
        raise ValueError(f"Invalid amino acid code: {res_code}. It must be either a three-letter or one-letter code.")

    return converted_code
