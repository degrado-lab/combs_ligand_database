import os
import sys
import pickle
import argparse
import traceback
import numpy as np
import prody as pr
import pandas as pd

from copy import deepcopy

from validation import *
from rotalyze import *
from probe import *

from probe import timeout

"""
Updated pdb files and validation reports should be downloaded via the 
pdb ftp server:

> rsync -rlpt -v -z --delete --port=33444 
  rsync.rcsb.org::ftp_data/structures/divided/pdb/ $LOCAL_PDB_MIRROR_PATH

> rsync -rlpt -v -z --delete --include="*/" --include="*.xml.gz" --exclude="*"  
  --port=33444 rsync.rcsb.org::ftp/validation_reports/ $LOCAL_VALIDATION_PATH

These two paths should be provided using the -e and -v arguments to this 
script, respectively.
"""

resnames_aa_20 = ['CYS', 'ASP', 'SER', 'GLN', 'LYS',
                  'ILE', 'PRO', 'THR', 'PHE', 'ASN',
                  'GLY', 'HIS', 'LEU', 'ARG', 'TRP',
                  'ALA', 'VAL', 'GLU', 'TYR', 'MET',
                  'MSE']
non_prot_sel = 'not resname ' + ' and not resname '.join(resnames_aa_20)
    
valid_cols = ['biounit', 'pdb_acc', 'R_free', 'R_obs', 'resolution', 
              'completeness', 'clashscore', 'abs_perc_rfree', 'abs_perc_rsrz', 
              'abs_perc_clashscore', 'abs_perc_rama_outl', 
              'abs_perc_rotamer_outl', 'rel_perc_rfree', 'rel_perc_clashscore', 
              'rel_perc_rsrz_outl', 'rel_perc_rama_outl', 
              'rel_perc_rotamer_outl', 'perc_rsrz_outl', 'perc_rama_outl', 
              'perc_rotamer_outl', 'ioversigi', 'fofc_corr', 'anisotropy', 
              'avg_bb_beta', 'sigma_bb_beta', 'ligs']


def get_segs_chains_resnums(atomgroup, selection, resnames=False):
    """Get a set containing tuples of segments, chains, and resnums.

    Parameters
    ----------
    atomgroup : prody.atomic.atomgroup.AtomGroup
        ProDy AtomGroup for a protein structure.
    selection : str
        String specifying subset of atoms in the ProDy selection algebra.
    resnames : bool
        If True, return residue names in each tuple as well.

    Returns
    -------
    segs_chains_resnums : set
        Set containing tuples of segments, chains, and resnums for 
        each residue that matches the selection string.
    """
    sel = atomgroup.select(selection)
    if sel is None:
        return set()
    if resnames:
        return set(zip(sel.getSegnames(), sel.getChids(), sel.getResnums(), 
                       sel.getResnames()))
    else:
        return set(zip(sel.getSegnames(), sel.getChids(), sel.getResnums()))


def get_author_assigned_biounits(pdb_file):
    """Given a gzipped PDB file, obtain the author-assigned biounits.

    Parameters
    ----------
    pdb_file : str
        Path to gzipped PDB file for which to obtain author-assigned biounits.

    Returns
    -------
    biomols : list
        List of integer indices of biological assemblies.
    """
    _str = 'AUTHOR DETERMINED BIOLOGICAL UNIT'
    biomols = []
    with gzip.open(pdb_file, 'rt') as infile:
        for line in infile:
            if line[:10] == 'REMARK 350':
                break
        for line in infile:
            if line[:10] != 'REMARK 350':
                break
            if 'BIOMOLECULE:' in line:
                biomol = int(line.strip().split('BIOMOLECULE:')[-1])
                continue
            if _str in line:
                biomols.append(biomol)
    return biomols


@timeout(5)
def get_bio(path):
    """Given the path to a gzipped PDB file, return its biological assemblies.

    Parameters
    ----------
    path : str
        Path to gzipped PDB file for which to return biological assemblies.

    Returns
    -------
    bio : prody.AtomGroup or list
        ProDy AtomGroup or list of ProDy AtomGroups for the biological 
        assemblies of the structure.
    """
    return pr.parsePDB(path, biomol=True)


def write_biounits(ent_gz_paths, pdb_tmp_dir, water_csv_path=None, 
                   max_ligands=None, write=True):
    """For a list of ent.gz files, write the author-assigned biounits to PDB.

    Parameters
    ----------
    ent_gz_paths : list
        List of paths to ent.gz files for PDB structures.
    pdb_tmp_dir : str
        Temporary directory at which to output unzipped ent files.
    water_csv_path : str
        Path to CSV containing water molecules to include, with columns 
        'pdb_code' (the entries of which are formatted pdbXXXX, with 
        XXXX being the four-letter accession code) and 'resnum'.
    max_ligands : int, optional
        Maximum number of heteroatom (i.e. non-protein, non-nucleic, and 
        non-water) residues to permit in a biological assembly.
    write : bool, optional
        If False, do not write the biological assemblies to PDB files.

    Returns
    -------
    bio_paths : list
        List of paths (within pdb_tmp_dir) to biounit PDB files.
    chain_pair_dicts : list
        List of dicts, one for each PDB file, that assign the original 
        chain ID from the asymmetric unit to each chain in the biounit.
    """
    if pdb_tmp_dir[-1] != '/':
        pdb_tmp_dir += '/'
    if water_csv_path:
        water_df = pd.read_csv(water_csv_path)
    else:
        water_df = None
    bio_paths = []
    chain_pair_dicts = []
    for i, path in enumerate(ent_gz_paths):
        try:
            bio = get_bio(path)
            pdb_code = path.split('/')[-1][:7]
            if type(bio) != list:
                bio = [bio]
            bio_list = [k + 1 for k in range(len(bio))]
            author_assigned = get_author_assigned_biounits(path)
            if len(author_assigned) > 0:
                bio_list = [int(b.getTitle().split()[-1]) for b in bio]
                bio = [bio[bio_list.index(i)] for i in author_assigned]
                bio_list = author_assigned
            for i, b in enumerate(bio):
                water_sel = 'not water'
                if water_df is not None:
                    subdf = water_df[water_df.pdb_code == pdb_code]
                    water_resnums = subdf['resnum'].values
                    if len(water_resnums):
                        water_sel += ' or resname HOH and resnum ' + \
                                     ' '.join([str(n) for n in water_resnums])
                bio[i] = b.select(water_sel).toAtomGroup()
            n_near = [len(get_segs_chains_resnums(b, 
                      non_prot_sel + ' within 4 of protein')) 
                      for b in bio]
            if type(max_ligands) is int:
                n_ligands = \
                    [len(get_segs_chains_resnums(b, 'not water hetero')) 
                     for b in bio]
                bio = [b for b, nl, nn in zip(bio, n_ligands, n_near) 
                       if nl < max_ligands and nc > 0]
            else:
                bio = [b for b, nn in zip(bio, n_near) if nn > 0] 
            for i, b in enumerate(bio):
                chids = b.getChids()
                segs = b.getSegnames() 
                new_chids = np.ones(len(chids), dtype='object')
                new_chids[:] = np.nan
                unique_seg_chids = sorted(set(list(zip(segs, chids))))
                for j, (seg, chid) in enumerate(unique_seg_chids):
                    if j < len(unique_seg_chids):
                        new_chids[(segs == seg) & (chids == chid)] = \
                            unique_seg_chids[j][1]
                    else:
                        break
                orig_chids = deepcopy(chids)
                chids[new_chids != np.nan] = new_chids[new_chids != np.nan]
                mask = (new_chids == np.nan)
                if mask.any():
                    print('***************************')
                    print(pdb, 'is more than 90 chains!')
                    chids[mask] = '?'
                    b.setChids(chids)
                    bio[i] = b.select('not chain ?').copy()
                else:
                    b.setChids(chids)
                chain_pair_dict = dict(zip(chids[~mask], orig_chids[~mask]))
                bio_path = pdb_tmp_dir + pdb_code[3:].upper() + '_biounit_' + \
                           str(bio_list[i]) + '.pdb'
                if write:
                    pr.writePDB(bio_path, bio[i])
                bio_paths.append(bio_path)
                chain_pair_dicts.append(chain_pair_dict)
        except Exception:
            print('**************************************************')
            traceback.print_exc(file=sys.stdout)
    return bio_paths, chain_pair_dicts


def reduce_pdbs(pdb_list, reduce_path, hetdict_path=None):
    """Add hydrogens to a list of PDB files using the Reduce program.

    Parameters
    ----------
    pdb_list : list
        List of paths to PDB files to be reduced.
    reduce_path : str
        Path to reduce binary.
    hetdict_path : str
        Path to het_dict specifying for the Reduce program how ligands 
        should be protonated.
    """
    for p in pdb_list:
        cmd = [reduce_path, '-TRIM', p, '>', p + '_trimreduce', 
               ';', reduce_path]
        if hetdict_path is not None:
            cmd += ['-DB', hetdict_path]
        cmd += ['-BUILD', p + '_trimreduce', '>', p + '_reduce', ';', 
                'rm', p + '_trimreduce', ';', 'mv', p + '_reduce', p]
        os.system(' '.join(cmd))


def transfer_O3(struct):
    """Transfer O3' atoms within nucleic acids to the next nucleotide.

    Parameters
    ----------
    struct : prody.atomic.atomgroup.AtomGroup
        The ProDy AtomGroup for which to transfer O3' atoms.

    Returns
    -------
    struct_trans : prody.atomic.atomgroup.AtomGroup
        The AtomGroup with O3' atoms transferred.
    """
    p_atoms = struct.select("element P")
    o3_atoms = struct.select("name O3'")
    if p_atoms is None or o3_atoms is None:
        return struct
    p_coords = p_atoms.getCoords()
    o3_coords = o3_atoms.getCoords()
    for at, at_coords in zip(o3_atoms, o3_coords):
        p_dists = np.linalg.norm(p_coords - at_coords, axis=1)
        if np.min(p_dists) < 2.:
            near_idx = int(np.argmin(p_dists))
            seg = p_atoms[near_idx].getSegname()
            chid = p_atoms[near_idx].getChid()
            resnum = p_atoms[near_idx].getResnum()
            resname = p_atoms[near_idx].getResname()
            if seg:
                at.setSegname(seg)
            at.setChid(chid)
            at.setResnum(resnum)
            at.setResname(resname)
            at.setName("OP3")
    return struct
 

def pdbs_to_pkl(pdb_list, chain_pair_dicts, validation_dir, 
                probe_outdir, rotalyze_outdir, prody_pkl_outdir, 
                validation_pkl_outdir, probe_path, rotalyze_path, 
                retry=False):
    """Use the Probe software to evaluate ligand contacts in a list of PDBs.
       Then, use the phenix.rotalyze software to evaluate rotamers for each.
       Last, execute DSSP on each and pickle the resulting ProDy objects.

    Parameters
    ----------
    pdb_list : list
        List of paths to PDB files to be probed.
    chain_pair_dicts : list
        List of dicts, one for each PDB file, that assign the original 
        chain ID from the asymmetric unit to each chain in the biounit.
    validation_dir : str
        Path to directory containing xml.gz files with validation report 
        data for all PDB structures.
    probe_outdir : str
        Path to directory at which to output pickled probe dataframes.
    rotalyze_outdir : str
        Path to directory at which to output pickled rotalyze dataframes.
    prody_pkl_outdir : str
        Path to directory at which to output pickled ProDy objects.
    validation_pkl_outdir : str
        Path to directory at which to output pickled dataframes containing 
        validation information, backbone B-factor statistics, and a list of 
        all non-protein residues for the biological assembly.
    probe_path : str
        Path to probe binary.
    rotalyze_path : str
        Path to rotalyze binary.
    retry : bool
        If True, run as if the code has already run but did not complete.
    """
    def resi_pm_n(resi, n=1): # get selection string for resi +/- n
        return '> {} and resnum < {}'.format(str(resi - n), str(resi + n))
    selstr_prot = 'protein within 4 of (segment {} and chain {} and resnum {})'
    selstr_prot_noseg = 'protein within 4 of (chain {} and resnum {})'
    selstr_general = '(segment {} and chain {} and resnum {})'
    selstr_general_noseg = '(chain {} and resnum {})'
    if probe_outdir[-1] != '/':
        probe_outdir += '/'
    if rotalyze_outdir[-1] != '/':
        rotalyze_outdir += '/'
    if prody_pkl_outdir[-1] != '/':
        prody_pkl_outdir += '/'
    if validation_pkl_outdir[-1] != '/':
        validation_pkl_outdir += '/'
    if validation_dir[-1] != '/':
        validation_dir += '/'
    n_pdbs = len(pdb_list)
    for n, p, chain_pair_dict in zip(range(n_pdbs), pdb_list, 
                                     chain_pair_dicts):
        print('pdb', str(n + 1), 'of', str(n_pdbs))
        val_data = []
        basename = '/'.join(p.split('/')[:-1]) + '/'
        name = p.split('/')[-1][:-4]
        pdb_name = name[:4].lower()
        if retry and name + '.pkl' in \
                os.listdir(prody_pkl_outdir + '/' + pdb_name[1:3] + '/'):
            continue
        ### create hierarchical subdirectories to ease filesystem load
        for outdir in [probe_outdir, rotalyze_outdir, prody_pkl_outdir, 
                       validation_pkl_outdir]:
            if not os.path.exists(outdir + '/' + pdb_name[1:3]):
                os.mkdir(outdir + '/' + pdb_name[1:3])
        ###
        val_path = validation_dir + '/' + pdb_name + '/' + pdb_name + \
                   '_validation.xml.gz'
        try:
            val_df, data = parse_validation_rep(val_path)
            b = pr.parsePDB(p)
            b = transfer_O3(b) # transfer O3' atoms in nucleotides
        except Exception:
            traceback.print_exc(file=sys.stdout)
            cmd = ['rm', basename + name + '.pdb;']
            os.system(' '.join(cmd))
            continue
        segs_chains_resnums_resnames = \
            get_segs_chains_resnums(b, non_prot_sel, True)
        if not len(segs_chains_resnums_resnames):
            cmd = ['rm', basename + name + '.pdb;']
            os.system(' '.join(cmd))
            continue
        contact_segs = [] # segment names of non-protein residues that contact 
                          # the protein according to probe
        contact_chids = [] # chain IDs of non-protein residues that contact 
                           # the protein according to probe
        prev_sites = [] # alphabetized list of resnames for a ligand and all 
                        # protein residues that contact it, used to remove 
                        # redundant protein-ligand contacts within a biounit
        selstrs = [] # selection strings which will collectively define the 
                     # subset of atoms that make it into the pickled AtomGroup
        lig_set = set()
        for seg_lig, chid_lig, resi_lig, resn_lig in segs_chains_resnums_resnames:
            # identify residues that contact the ligands 
            resi_lig_for_sel = resi_pm_n(resi_lig, 1)
            if seg_lig:
                selstrs.append(selstr_general.format(seg_lig, chid_lig, 
                                                     resi_lig_for_sel))
                scrr_prot = get_segs_chains_resnums(b, 
                    selstr_prot.format(seg_lig, chid_lig, resi_lig_for_sel), 
                    True)
            else:
                selstrs.append(selstr_general_noseg.format(chid_lig, 
                                                           resi_lig_for_sel))
                scrr_prot = get_segs_chains_resnums(b, 
                    selstr_prot_noseg.format(chid_lig, resi_lig_for_sel), 
                    True)
            if not len(scrr_prot):
                continue
            _, _, aa_resi, aa_resn = [list(tup) for tup in zip(*scrr_prot)]
            site = sorted([resn_lig] + aa_resn + 
                          [str(resi_lig)] + [str(resi) for resi in aa_resi])
            if site in prev_sites:
                continue
            else:
                prev_sites.append(site)
            # run probe software to analyze binding site
            try:
                probe_df = parse_probe(p, segname1=seg_lig, 
                                       chain1=chid_lig, 
                                       resnum1=resi_lig, 
                                       path_to_probe=probe_path)
            except TimeoutError:
                break
            if not len(probe_df):
                continue
            probe_df.loc[probe_df.chain1 == chid_lig, 'orig_chain1'] = \
                chain_pair_dict[chid_lig]
            probe_df = \
                pd.merge(probe_df, val_df, 
                         left_on=['orig_chain1', 'resname1', 'resnum1'], 
                         right_on=['chain', 'resname', 'resnum'], 
                         suffixes=['', '1'])
            if not len(probe_df):
                continue
            for chain2 in probe_df.chain2.unique():
                probe_df.loc[probe_df.chain2 == chain2, 'orig_chain2'] = \
                    chain_pair_dict[chain2]
            probe_df = \
                pd.merge(probe_df, val_df, 
                         left_on=['orig_chain2', 'resname2', 'resnum2'], 
                         right_on=['chain', 'resname', 'resnum'], 
                         suffixes=['1', '2'])
            if not len(probe_df):
                continue
            probe_df = probe_df.loc[:, ~probe_df.columns.duplicated()]
            probe_outpath = probe_outdir + '/' + pdb_name[1:3] + '/' + name + \
                            '_' + seg_lig + '_' + chid_lig + '.pkl'
            if os.path.exists(probe_outpath):
                with open(probe_outpath, 'rb') as inpkl:
                    prev_df = pickle.load(inpkl)
                probe_df = pd.concat([prev_df, probe_df], ignore_index=True)
            probe_df.to_pickle(probe_outpath)
            lig_set.add(probe_df.resname1.iloc[0])
            # check which residues contact the ligand according to probe
            for seg_prot, chid_prot, resi_prot, resn_prot in scrr_prot:
                if not ((probe_df['chain2'] == chid_prot) & 
                        (probe_df['resnum2'] == resi_prot)).any():
                    continue
                contact_chids.append(chid_prot)
                contact_segs.append(seg_prot)
                resi_prot_for_sel = resi_pm_n(resi_prot, 2)
                if seg_prot:
                    _str = selstr_general.format(seg_prot, chid_prot,
                                                 resi_prot_for_sel)
                else:
                    _str = selstr_general_noseg.format(chid_prot, 
                                                       resi_prot_for_sel)
                if _str not in selstrs:
                    selstrs.append(_str)
        if len(contact_chids):
            avg_bb_beta = b.select('name N CA C O').getBetas().mean()
            sigma_bb_beta = b.select('name N CA C O').getBetas().std()
            ligs = ' '.join(lig_set)
            val_data.append(tuple([name] + data + 
                                  [avg_bb_beta, sigma_bb_beta, ligs]))
            parse_rotalyze(p, rotalyze_path, 
                           rotalyze_outdir + '/' + pdb_name[1:3])
            pr.execDSSP(p, outputdir=basename)
            pr.parseDSSP(basename + name + '.dssp', b)
            if len(contact_segs) and '' not in contact_segs:
                final_selstr = ' or '.join(selstrs) + ' or (name CA CB and ('
                addends = []
                for seg_prot, chid_prot in zip(contact_segs, contact_chids):
                    addends.append(
                        ('(segment {} and chain {})').format(seg_prot, 
                                                             chid_prot))
                final_selstr += ' or '.join(addends) + '))'
            else:
                final_selstr = ' or '.join(selstrs) + \
                               ' or (chain {} and name CA CB)'
                contact_chids = ' '.join(set(contact_chids))
                final_selstr = final_selstr.format(contact_chids)
            compressed_b = b.select(final_selstr).toAtomGroup()
            with open(prody_pkl_outdir + '/' + pdb_name[1:3] + '/' + \
                      name + '.pkl', 'wb') as outfile:
                pickle.dump(file=outfile, obj=compressed_b)
            validation_df = pd.DataFrame(val_data, columns=valid_cols)
            with open(validation_pkl_outdir + '/' + pdb_name[1:3] + '/' + \
                      name + '.pkl', 'wb') as outfile:
                pickle.dump(file=outfile, obj=validation_df)
        if os.path.exists(basename + name + '.dssp'):
            cmd = ['rm', basename + name + '.pdb;',
                   'rm', basename + name + '.dssp']
        else:
            cmd = ['rm', basename + name + '.pdb']
        os.system(' '.join(cmd))


def ent_gz_dir_to_combs_db_files(ent_gz_dir, validation_dir, 
                                 prody_pkl_outdir, rotalyze_outdir, 
                                 probe_outdir, validation_outdir, pdb_tmp_dir, 
                                 reduce_path, probe_path, rotalyze_path, 
                                 max_ligands=25, water_csv_path=None, 
                                 pdb_het_dict=None, retry=False):
    """Generate input files for COMBS database generation from ent.gz files.

    Parameters
    ----------
    ent_gz_dir : str
        Path to directory containing ent.gz files from which to generate input 
        files for COMBS database generation.
    validation_dir : str
        Path to directory containing xml.gz files with validation report 
        data for all PDB structures.
    prody_pkl_outdir : str
        Path to directory at which to output pickled ProDy files.
    rotalyze_outdir : str
        Path to directory at which to output pickled rotalyze dataframes.
    probe_outdir : str
        Path to directory at which to output pickled probe dataframes.
    validation_outdir : str
        Path to directory at which to output pickled validation dataframes.
    pdb_tmp_dir : str
        Temporary directory at which to output unzipped ent files.
    reduce_path : str
        Path to reduce binary.
    probe_path : str
        Path to probe binary.
    rotalyze_path : str
        Path to rotalyze binary.
    max_ligands : int
        Maximum number of heteroatom (i.e. non-protein, non-nucleic, and 
        non-water) residues to permit in a biological assembly.
    water_csv_path : str
        Path to CSV containing water molecules to include, with columns 
        'pdb_code' (the entries of which are formatted pdbXXXX, with 
        XXXX being the four-letter accession code) and 'resnum'.
    pdb_het_dict : str
        Path to het_dict specifying for the Reduce program how ligands 
        should be protonated.
    retry : bool
        Run as if the code has already been run but did not complete.
    """
    for _dir in [ent_gz_dir, prody_pkl_outdir, rotalyze_outdir, 
                 probe_outdir, validation_outdir, pdb_tmp_dir]:
        if not os.path.exists(_dir):
            os.mkdir(_dir)
    if ent_gz_dir[:-1] != '/':
        ent_gz_dir += '/'
    ent_gz_paths = [ent_gz_dir + path for path in os.listdir(ent_gz_dir)]
    bio_paths, chain_pair_dicts = write_biounits(ent_gz_paths, pdb_tmp_dir, 
                                                 water_csv_path, 
                                                 write=(not retry))
    if not retry:
        reduce_pdbs(bio_paths, reduce_path, pdb_het_dict)
    pdbs_to_pkl(bio_paths, chain_pair_dicts, validation_dir, 
                probe_outdir, rotalyze_outdir, prody_pkl_outdir, 
                validation_outdir, probe_path, rotalyze_path, retry)


def parse_args():
    argp = argparse.ArgumentParser()
    argp.add_argument('-e', '--ent-gz-dir', help="Path to directory "
                      "containing ent.gz files from which to generate input "
                      "files for COMBS database generation.")
    argp.add_argument('-v', '--validation-dir', help="Path to directory "
                      "containing xml.gz files with validation report "
                      "data for all PDB structures.")
    argp.add_argument('-o', '--prody-pkl-outdir', help="Path to directory at "
                      "which to output pickled ProDy files.")
    argp.add_argument('-r', '--rotalyze-outdir', help="Path to directory at "
                      "which to output pickled rotalyze dataframes.")
    argp.add_argument('-p', '--probe-outdir', help="Path to directory at "
                      "which to output pickled probe dataframes.")
    argp.add_argument('-a', '--validation-outdir', help="Path to directory at "
                      "which to output pickled validation dataframes.")
    argp.add_argument('-t', '--pdb-tmp-dir', 
                      default='/wynton/scratch/rian.kormos/tmp_pdbs/', 
                      help="Temporary directory at which to output unzipped "
                      "ent files.")
    argp.add_argument("--reduce-path", help="Path to reduce binary.")
    argp.add_argument("--probe-path", help="Path to probe binary.")
    argp.add_argument("--rotalyze-path", help="Path to rotalyze binary.")
    argp.add_argument('-m', '--max-ligands', type=int, default=25, 
                      help="Maximum number of heteroatom (i.e. non-protein, "
                      "non-nucleic, and non-water) residues to permit in a "
                      "biological assembly.")
    argp.add_argument('-w', '--water-csv-path', 
                      help="Path to CSV file containing water molecules to "
                      "be added to the database (optional, see docstrings for "
                      "more information).")
    argp.add_argument('-d', '--pdb-het-dict', help="Path to het_dict "
                      "specifying for the Reduce program how ligands "
                      "should be protonated (optional).")
    argp.add_argument('--retry', action='store_true', 
                      help="Run as if the code has already been run but "
                      "did not complete (i.e. finish generating the files "
                      "that did not generate in an earlier run.")
    return argp.parse_args()


if __name__ == "__main__":
    args = parse_args()
    _reduce = '/wynton/home/degradolab/rkormos/reduce/reduce_src/reduce'
    _probe = '/wynton/home/degradolab/rkormos/probe/probe'
    _rotalyze = ('/salilab/diva1/programs/x86_64linux/phenix-1.19.1.4122/'
                 'phenix-1.19.1-4122/build/bin/phenix.rotalyze')
    if args.reduce is None:
        args.reduce = _reduce
    if args.probe is None:
        args.probe = _probe
    if args.rotalyze is None:
        args.rotalyze = _rotalyze
    ent_gz_dir_to_combs_db_files(args.ent_gz_dir, args.validation_dir, 
                                 args.prody_pkl_outdir, args.rotalyze_outdir, 
                                 args.probe_outdir, args.validation_outdir, 
                                 args.pdb_tmp_dir, args.reduce_path, 
                                 args.probe_path, args.rotalyze_path, 
                                 args.max_ligands, args.water_csv_path, 
                                 args.pdb_het_dict, args.retry)
