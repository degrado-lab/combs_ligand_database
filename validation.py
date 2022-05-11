"""Unfortunately the validation xml files from the RCSB
appear to contain errors in occupancy. Either that or I don't
understand it.  Either way, omitting occupancy here."""
import gzip
import pandas as pd
import xml.etree.ElementTree as et


#dtypes of dataframe columns in :func:`parse_validation_rep`
val_dtype_dict = {'chain': 'category',
                  'resname': 'category',
                  'resnum': int,
                  'altcode': 'category',
                  'icode': 'category',
                  'rscc': float,
                  'rsr': float,
                  'rsrz': float,
                  'clash': bool,
                  'phi': float,
                  'psi': float,
                  'rama': 'category',
                  'rota': 'category'}


# dataframe columns in :func:`parse_validation_rep`
val_col_names = list(val_dtype_dict.keys())


def get_data(tree):
    e = tree.findall('Entry')[0].attrib
    R = e['PDB-R']
    RFREE = e['PDB-Rfree']
    COMPLETENESS = e['DataCompleteness']
    RESOLUTION = e['PDB-resolution']
    FOFC_CORR = e['Fo_Fc_correlation']
    IOVERSIGI = e['IoverSigma'].split('(')[0]
    ANISO = e['DataAnisotropy']
    ABS_PER_RFREE = e['absolute-percentile-DCC_Rfree']
    ABS_PER_RSRZ = e['absolute-percentile-percent-RSRZ-outliers']
    ABS_PER_CLASH = e['absolute-percentile-clashscore']
    ABS_PER_RAMA = e['absolute-percentile-percent-rama-outliers']
    ABS_PER_ROTA = e['absolute-percentile-percent-rota-outliers']
    CLASHSC = e['clashscore']
    PDBACC = e['pdbid']
    PER_RSRZ = e['percent-RSRZ-outliers']
    PER_RAMA = e['percent-rama-outliers']
    PER_ROTA = e['percent-rota-outliers']
    REL_PER_RFREE = e['relative-percentile-DCC_Rfree']
    REL_PER_CLASH = e['relative-percentile-clashscore']
    REL_PER_RSRZ = e['relative-percentile-percent-RSRZ-outliers']
    REL_PER_RAMA = e['relative-percentile-percent-rama-outliers']
    REL_PER_ROTA = e['relative-percentile-percent-rota-outliers']
    return (PDBACC, RFREE, R, RESOLUTION, COMPLETENESS,
            CLASHSC, ABS_PER_RFREE, ABS_PER_RSRZ, ABS_PER_CLASH,
            ABS_PER_RAMA, ABS_PER_ROTA, REL_PER_RFREE, REL_PER_CLASH,
            REL_PER_RSRZ, REL_PER_RAMA, REL_PER_ROTA, PER_RSRZ,
            PER_RAMA, PER_ROTA, IOVERSIGI, FOFC_CORR, ANISO)


def parse_validation_rep(path_to_rep):
    """

    Parameters
    ----------
    path_to_rep

    Returns
    -------

    """
    with gzip.open(path_to_rep, 'rb') as infile:
        tree = et.parse(infile)
    attribs = []
    for ms in tree.findall('ModelledSubgroup'):
        children = list(ms)
        clash = False
        if children:
            for child in children:
                if child.tag == 'clash':
                    clash = True
        attrib = ms.attrib
        attrib['clash'] = clash
        attribs.append(attrib)
    validation_df = pd.DataFrame(attribs, columns=val_col_names)
    data = get_data(tree)
    return validation_df.astype(dtype=val_dtype_dict), data


def get_data(tree):
    e = tree.findall('Entry')[0].attrib
    return_vals = []
    for key in ['pdbid', 'PDB-R', 'PDB-Rfree', 'PDB-resolution', 
                'DataCompleteness', 'clashscore', 
                'absolute-percentile-DCC_Rfree', 
                'absolute-percentile-percent-RSRZ-outliers', 
                'absolute-percentile-clashscore', 
                'absolute-percentile-percent-rama-outliers', 
                'absolute-percentile-percent-rota-outliers', 
                'relative-percentile-DCC_Rfree', 
                'relative-percentile-clashscore', 
                'relative-percentile-percent-RSRZ-outliers', 
                'relative-percentile-percent-rama-outliers', 
                'relative-percentile-percent-rota-outliers',
                'percent-RSRZ-outliers', 
                'percent-rama-outliers', 'percent-rota-outliers', 
                'IoverSigma', 'Fo_Fc_correlation', 'DataAnisotropy']:
        try:
            return_vals.append(e[key])
        except:
            return_vals.append(None)
    return return_vals
