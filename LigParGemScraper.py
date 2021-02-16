from SecurePathProxy import SecurePathProxy
import re
import requests
import os
class LigParGemScraper:
    """
    This class acts like a small state machine. It connects to the LigParGem Website. 
    Update whathever you want and call the functions make_connection and get_values.
    You can pass any list of values to the get_values and the code will get it for your molecule.
    Remember you need to choose between the PDB file or the SMILES string.
    """
    def __init__(self, file_path = '.'):
        #The input comes with the default benzene:
        self.input_dict = {'smiData': 'c1ccccc1', 
                             'chargetype': 'cm1a',
                             'dropcharge': ' 0 ',
                             'checkopt': ' 0 ',
                             'submit':'Submit'}
        #And without any file:
        self.__use_smiles()
        self.charge_types = ['cm1abcc', 'cm1a']
        self.file_path = file_path
        
    def change_charge_type(self, charge_type_index):
        assert charge_type_index < len(self.charge_types) - 1
        self.update_input_dict({'chargetype' : self.charge_types[charge_type_index] })
        
        
    def update_input_dict(self, new_dict):
        for k in new_dict:
            self.input_dict[k] = new_dict[k]
    
    def change_molecule_optimization_interations(self, new_interations):
        assert( isinstance(new_interations, int) and new_interations <=3 and new_interations >= 0)
        self.update_input_dict({'checkopt': f' {new_interations} '})
        
    def change_molecule_charge(self, new_charge):
        assert(isinstance(new_charge, int) and new_charge >= -2 and new_charge <= 2)
        self.update_input_dict({'dropcharge': f' {new_charge} '})
        
    def update_smiles(self, new_smiles):
        assert(isinstance(new_smiles, str))
        self.__use_smiles()
        self.update_input_dict({'smiData': new_smiles })
    
    def update_pdb(self, pdb_name):
        with SecurePathProxy(os, self.file_path) as p:
            assert (os.path.exists( pdb_name))
            self.__use_pdb(pdb_name)
        
    def make_connection(self):
        """
        This function make the first connection to the LibParGen website.
        """
        self.s = requests.Session()
        self.s.get('http://zarbi.chem.yale.edu/ligpargen/')
        self.s.headers.update({'Host' : 'zarbi.chem.yale.edu',
                 'Origin' : 'http://zarbi.chem.yale.edu',
                 'Referer' : 'http://zarbi.chem.yale.edu/ligpargen/',
                 'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
                 'User-Agent' : 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/88.0.4324.150 Safari/537.36 Edg/88.0.705.68'})
        
        r = self.s.post('http://zarbi.chem.yale.edu/cgi-bin/results_lpg.py', data = self.input_dict, files= self.files )
        self.__extract_file_out(r)
        
        return r
    
    def get_values(self, value_list = ['PDB', 'XML', 'GRO']):
        """
        This function should be called only after make_connection. The default list of values are ['PDB', ' XML', 'GRO'].
        You can override it to whatever you want.
        """
        assert(hasattr(self, 'file_out')) #This function should only be called after make_connection this assert garantees it
        
        self.s.headers.update({'Referer' : 'http://zarbi.chem.yale.edu/cgi-bin/results_lpg.py'})
        for v in value_list:
            input_dict = {'fileout':self.file_out + f'.{v.lower()}',
                         'go': value_list}
            r = self.s.post('http://zarbi.chem.yale.edu/cgi-bin/download_lpg.py', data = input_dict)
            self.__save_file(r.text, v)
            
    def set_file_path(self, file_path):
        self.file_path = file_path
    
    def __use_smiles(self):
        self.files = {'molpdbfile': ('', '', 'application/octet-stream')}
    
    def __use_pdb(self, pdb_name):
        self.update_input_dict({'smiData': '' })
        self.files = {'molpdbfile': (pdb_name, open(pdb_name, 'rb'), 'application/octet-stream')}
        
    def __extract_file_out(self, r):
        self.file_out = re.search('name="fileout" value="(.+)"', r.text)[1].split('.')[0]
        
        
    def __save_file(self, file_content, extension):
            with SecurePathProxy(os, self.file_path) as p :
                file_name = self.file_out.split('/tmp/')[-1] + f'.{extension.lower()}' if not hasattr(self,'file_name') else self.file_name + f'.{extension.lower()}'
                with open(file_name, 'w') as f:
                    f.write(file_content)



        