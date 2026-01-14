from collections import OrderedDict


def parse_ndx(ndx_file):
    """
    Parse a GROMACS index (.ndx) file and return a dictionary of groups.

    Parameters
    ----------
    ndx_file : str
        Path to the GROMACS index file.

    Returns
    -------
    groups : dict
        An ordered dictionary where keys are group names and values are lists of atom indices.
    group_str : str
        A string representation of the groups for easy viewing.
    """
    groups = OrderedDict()
    current_group = None
    with open(ndx_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[') and line.endswith(']'):
                current_group = line[1:-1].strip()
                groups[current_group] = []
            elif current_group is not None:
                indices = line.split()
                groups[current_group].extend(int(idx) for idx in indices)

    name_w = max((len(name) for name in groups), default=0)
    label_w = name_w + 1
    count_w = len(str(max((len(atoms) for atoms in groups.values()), default=0)))
    idx_w = len(str(len(groups) - 1))

    group_str = ''
    for i, (name, atoms) in enumerate(groups.items()):
        label = f"{name}:"
        group_str += f"  ({i}) {'':<{idx_w - len(str(i))}}{label:<{label_w}} {len(atoms):>{count_w}} atoms\n"

    return groups, group_str
