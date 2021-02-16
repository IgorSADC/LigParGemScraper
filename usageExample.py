from LigParGemScraper import LigParGemScraper
import requests

def get_smiles_from_name(name):
    '''
    A small function to get the smiles from the molecule name using the PUGREST API
    '''
    return requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/TXT").text.split('\n')[0]
    
l = LigParGemScraper(file_path = 'Triclosan')
l.update_smiles(get_smiles_from_name("Triclosan"))
r = l.make_connection()
l.file_name = 'Triclosan'
l.get_values(['PDB'])