'''Functions to define various groups of sites of interest in NP. Written by Mike Doud, updated March 2015.'''

from Bio.PDB import *

def RNABindingGrooveSites():
    '''These are the sites inferred to be the RNA-binding groove of NP in Ye et al, Nature 2006.'''
    return [65,150,152,156,174,175,195,199,213,214,221,236,355,357,361,391,148,184,216]

def VariableSites():
    """There are the sites that have different amino-acid identities between Aichi/1968 and PR/1934.
    The amino-acid sequences in *aichi68_identity* and *pr8_identity* are from translated nucleotide 
    sequences from the plasmid maps for the PR8 and Aichi68 NP plasmids used in these experiments.
    """
    different_sites = []
    for i in range(2,499):
        if aichi68_identity(i) != pr8_identity(i):
            different_sites.append(i)
    return different_sites

def aichi68_identity(site):
    '''Returns the amino-acid identity of Aichi68 NP at position `site` (indexed as 1,2,3...)'''
    site = int(site)
    aichi68_seq = "MASQGTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYKRVDRKWMRELVLYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDTTYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRKTRSAYERMCNILKGKFQTAAQRAMMDQVRESRNPGNAEIEDLIFLARSALILRGSVAHKSCLPACVYGPAVASGYDFEKEGYSLVGIDPFKLLQNSQVYSLIRPNENPAHKSQLVWMACNSAAFEDLRVLSFIRGTKVSPRGKLSTRGVQIASNENMDAMESSTLELRSRYWAIRTRSGGNTNQQRASAGQISVQPAFSVQRNLPFDKPTIMAAFTGNTEGRTSDMRAEIIRMMEGAKPEEMSFQGRGVFELSDERAANPIVPSFDMSNEGSYFFGDNAEEYDN"
    return aichi68_seq[site-1]

def pr8_identity(site):
    '''Returns the amino-acid identity of PR8 NP at position `site` (indexed as 1,2,3...)'''
    site = int(site)
    pr8_seq =     "MASQGTKRSYEQMETDGERQNATEIRASVGKMIGGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRVNGKWMRELILYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELVRMIKRGINDRNFWRGENGRKTRIAYERMCNILKGKFQTAAQKAMMDQVRESRNPGNAEFEDLTFLARSALILRGSVAHKSCLPACVYGPAVASGYDFEREGYSLVGIDPFRLLQNSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVLSFIKGTKVLPRGKLSTRGVQIASNENMETMESSTLELRSRYWAIRTRSGGNTNQQRASAGQISIQPTFSVQRNLPFDRTTIMAAFNGNTEGRTSDMRTEIIRMMESARPEDVSFQGRGVFELSDEKAASPIVPSFDMSNEGSYFFGDNAEEYDN"
    return pr8_seq[site-1]

def SitesContactingVariableSites():
    '''Returns a list of all sites that are in contact with variable sites (can include other variable sites)
    and a list of conserved sites that are in contact with variable sites (ie the first list but excluding variable sites)'''

    # to replicate this, specify the full path to the 2IQH_monomerC.pdb file and run the enclosed code:
    '''
    pdbfile = '2IQH_monomerC.pdb'
    dist_cutoff = 4.5
    contacts = CalculateContacts(pdbfile, dist_cutoff)
    variable_sites = VariableSites()

    all_contact_sites_list = [] # all sites (including variable and conserved) that are in contact with any site that is variable between Aichi68 and PR8
    conserved_contact_sites_list = [] # conserved sites that are in structural contact with any site that is variable between Aichi68 and PR8

    for diffsite in variable_sites:
        if diffsite in contacts:
            unadded_contacts = [c for c in contacts[diffsite] if c not in all_contact_sites_list]
            for c in unadded_contacts:
                if c not in variable_sites:
                    conserved_contact_sites_list.append( c )
                all_contact_sites_list.append( c )

    return (conserved_contact_sites_list, all_contact_sites_list)
    '''
    return ([30, 35, 38, 290, 291, 292, 47, 49, 53, 105, 107, 99, 66, 94, 96, 111, 328, 331, 332, 337, 157, 190, 197, 256, 221, 227, 228, 230, 231, 254, 159, 261, 301, 302, 306, 469, 274, 304, 307, 330, 333, 326, 327, 349, 350, 352, 380, 383, 324, 354, 357, 358, 103, 370, 318, 374, 405, 418, 419, 412, 420, 426, 446, 449, 421, 439, 453, 454, 451, 459, 468, 472, 474], [30, 35, 38, 290, 291, 292, 47, 49, 53, 105, 107, 99, 66, 94, 96, 111, 328, 331, 332, 337, 157, 190, 197, 253, 256, 221, 227, 228, 230, 231, 194, 254, 257, 159, 261, 301, 302, 306, 469, 274, 304, 307, 330, 333, 326, 327, 349, 350, 352, 380, 383, 324, 354, 357, 358, 103, 370, 318, 374, 405, 418, 419, 412, 420, 423, 426, 446, 449, 421, 422, 450, 439, 453, 454, 455, 451, 452, 459, 468, 472, 474])

def CalculateContacts(pdbfile, dist_cutoff):
    """Calculates all residues in contact with one another given a PDB structure.
    Residues are in contact if their distance apart is less than dist_cutoff.

    pdbfile -> PDB structure file for which to find residue pairs forming contacts
    dist_cutoff -> Maximum distance a pair of residues can be separated by and still
       consider that residue pair in contact. Typically set cutoff to 7 A for comparisons
       between beta carbons, but smaller thresholds (4.5 A) can be used here for all possible 
       distances between sidechain atoms distal to the alpha carbon.
    contacts -> This is a dictionary of lists storing residues making contacts.
       contacts[res1] = [res2, res3...], where res# are integer
       residue numers from the PDB structure. res1 is in contact with res2 and res3.

    Written by Orr Ashenberg.
    """
    contacts = {}
    distances = CalculateDistances(pdbfile)
    for res1 in distances:
        contacts[res1] = []
        for res2 in distances[res1]: 
            if distances[res1][res2] < dist_cutoff:
                contacts[res1].append(res2)
    return contacts

def CalculateDistances(pdbfile):
    """Calculate distances between residues given a PDB structure.

    The distance separating the residues within a pair of residues is calculated using the (x,y,z)
    coordinates for all sidechain atoms in each residue. 

    For each residue, the distances are calculated using all of the atoms distal to the CA atom.
    If a residue is glycine, the CA atom is used instead. 

    The distances between the same residues are not calculated (distance between residue x and residue x). 
    This function uses the PDB module from BioPython to parse the PDB structure and to measure the distance between residues.

    The distance returned for each pair of residues is the minimum distance between any pair of sidechain atoms considered for
    those residues.

    pdbfile -> PDB structure file for which to measure pairwise residue distances
    distances -> Dictionary containing residue-residue pairwise distances.
        distances[resnum1][resnum2] = distance, where resnum1 and resnum2 are integer residue
        numbers from the PDB file. The dictionary is symmetric, ie distances[resnum1][resnum2] = 
        distances[resnum2][resnum1]. The distances dictionary is the return value of this function.

    Written by Orr Ashenberg, September 2014. Modified by Mike Doud, October 2014, to calculate the minimum distance between
    all possible sidechain atoms for each pair of residues.
    """
    parser = PDBParser()
    structure = parser.get_structure('temp', pdbfile)

    distances = {}
    for residue1 in structure.get_residues():
        resname1 = residue1.get_resname()
        resnum1 = residue1.get_id()[1]
        distances[resnum1] = {}

        # For Glycine we will measure distance from the alpha carbon; for all other residues we will consider all distances
        # from all atoms distal to the alpha carbon on the side chain:
        if resname1 == 'GLY':
            atomids_1_to_check = ['CA']
        else:
            atomids_1_to_check = [ a.get_name() for a in residue1.get_unpacked_list() if a.get_name() not in ['N','CA','O','C'] ]

        #print "resname1 = %s and we are checking atoms %s" % (resname1, atomids_1_to_check)
        for residue2 in structure.get_residues():
            if residue1 is not residue2:
                resname2 = residue2.get_resname()
                resnum2 = residue2.get_id()[1]
                if resname2 == 'GLY':
                    atomids_2_to_check = ['CA']
                else:
                    atomids_2_to_check = [ a.get_name() for a in residue2.get_unpacked_list() if a.get_name() not in ['N','CA','O','C'] ]

                #print "resname2 = %s and we are checking atoms %s" % (resname2, atomids_2_to_check)
                try:
                    this_residue_pair_distances = [] # hold all pairwise distances for the atoms under consideration in these two residues
                    for atom1 in atomids_1_to_check:
                        for atom2 in atomids_2_to_check:
                            this_residue_pair_distances.append( residue1[atom1] - residue2[atom2] )
                    distances[resnum1][resnum2] = min(this_residue_pair_distances)
                    #print resname1, residue1, resnum1, atomid1, resname2, residue2, resnum2, atomid2, distances[resnum1][resnum2]

                except KeyError:
                    warnings.warn("Missing CA or CB atom in CalculateDistances %s %d %s %s %d %s %s" \
                                      % (resname1, resnum1, atomid1, resname2, resnum2, atomid2, pdbfile))
    return distances
