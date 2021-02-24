import pandas as pd
class LigParGenTxtGenerator:
    def __init__(self, original_pdb, ligpargen_pdb, ligpargen_itp):
        self.original_pdb, self.ligpargen_pdb, self.ligpargen_itp = original_pdb, ligpargen_pdb, ligpargen_itp
        self.create_hash_map()
        self.convert_dict = {'H': 6, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C':6, 'N':7, 'O':8, 'F':9, 'S':16, 'Cl':17}
    
    def create_hash_map(self):
        self.mol_db = pd.DataFrame(columns=['REMARK', 'ID', 'ATOM', 'ID_MOL', 'X', 'Y', 'Z', '_'])
        self.mol_ligpargen_db =  pd.DataFrame(columns=['REMARK', 'ID', 'ATOM', '_','ID_MOL', 'X', 'Y', 'Z'])
        self.__read_original_pdb()
        self.__read_ligpargen_pdb()
        self.mol_ligpargen_db['ATOM_IN_GAUSSIAN'] =  self.mol_ligpargen_db['ATOM'].str.extract('([A-Z a-z])')
        self.mapped_id = {int(self.mol_db.loc[c].ID) : self.__extract_mapping(int(c)) for c in self.mol_db.index}
        self.mol_db['OPLS_ID'] = self.mol_db.apply(lambda x : self.mapped_id[int(x.ID)], axis=1)
        self.get_values_from_itp()
        
    def get_values_from_itp(self):
        opls_db, atoms_db = self.__read_itp()
        temporary_df = pd.DataFrame(columns=['EPSILON', 'SIGMA'])
        for i in range(self.mol_db.shape[0]):
            temporary_df.loc[temporary_df.shape[0]] = self.__get_eplison_sigma(opls_db, atoms_db, self.mol_db.iloc[i])
        self.mol_db = self.mol_db.join(temporary_df)
        
    def save_txt_file(self, title = 'SomeDiceIput', filename='molecule.txt', charges = None):
        self.mol_db['Charges'] = float(0) if charges == None else charges
        with open(filename, 'w') as new_file:
            new_file.write(f'1 {title} \n')
            new_file.write(f'{self.mol_db.shape[0]} na x y z q epsilon sigma \n')
            for i in range(self.mol_db.shape[0]):
                current_row = self.mol_db.iloc[i]
                new_file.write(f'{current_row.ID} {self.__get_atomic_number_from_symbol(current_row.ATOM)} {current_row.X} {current_row.Y} {current_row.Z} {current_row.Charges} {current_row.EPSILON:.3f} {current_row.SIGMA} \n')
                
            new_file.write('$end')
                
    
    def __get_atomic_number_from_symbol(self, s):
        '''TEMPORARY FUNCTION. TO BE MOVED TO MOLEKING'''
        return self.convert_dict[s]
        

    def __get_eplison_sigma(self, opls_db, atoms_db, row):
        ID = int(row.ID)
        c_value = opls_db[opls_db['OPLS_ID'] == atoms_db[atoms_db['ATOM'] == self.mapped_id[ID]]['TYPE'].iloc[0]]
        return float(c_value['EPSILON'].values[0]) * 0.23901, float(c_value['SIGMA'].values[0]) * 10
        
        
    
    def __read_itp(self):
        opls_db = pd.DataFrame(columns = ['OPLS_ID', 'ATOM', 'MASS', '_', '_', 'SIGMA', 'EPSILON'])
        atoms_db = pd.DataFrame(columns = ['NR', 'TYPE', 'RESNR', 'RESIDUE', 'ATOM', 'CGNR', 'CHARGE', 'MASS'])
        with open(self.ligpargen_itp, 'r') as file:
            capture_opls = False
            capture_atoms = False
            for i, c in enumerate(file):
                if('[ moleculetype ]' in c): capture_opls = False
                if(capture_opls):
                    opls_db.loc[opls_db.shape[0]] =  c.split()

                elif(capture_atoms):
                    if(len(c.split())!= atoms_db.shape[1]): continue
                    atoms_db.loc[atoms_db.shape[0]] =  c.split()

                if('[ atomtypes ]' in c): capture_opls = True

                if('[ atoms ]' in c): capture_atoms = True
                if('[ bonds ]' in c): 
                    capture_atoms = False
                    break
        return (opls_db, atoms_db)
        
        
        
    def __read_original_pdb(self):
        with open(self.original_pdb, 'r') as file:
            for i, c in enumerate(file):
                if('END' in c): break
                if not 'HETATM' in c: continue
                self.mol_db.loc[self.mol_db.shape[0]] =  c.split()
                
    def __read_ligpargen_pdb(self):
        with open(self.ligpargen_pdb, 'r') as file:
            for i, c in enumerate(file):
                if('END' in c): break
                if not 'ATOM' in c: continue
                self.mol_ligpargen_db.loc[self.mol_ligpargen_db.shape[0]] =  c.split()

    def __extract_mapping(self, loc_value, criterium = 0.1):
        atom = self.mol_db.iloc[loc_value].ATOM
        x =  float(self.mol_db.loc[loc_value].X)
        y =  float(self.mol_db.loc[loc_value].Y)
        z =  float(self.mol_db.loc[loc_value].Z)
        targets = self.mol_ligpargen_db[self.mol_ligpargen_db['ATOM_IN_GAUSSIAN'] == atom]
        targets = targets[abs(targets['X'].astype(float) - x) <= criterium ]
        if(targets.shape[0] == 1): return targets.iloc[0]['ATOM']
        targets = targets[abs(targets['Y'].astype(float) - y) <= criterium ]
        if(targets.shape[0] == 1): return targets.iloc[0]['ATOM']
        targets = targets[abs(targets['Z'].astype(float) - z) <= criterium ]
        if(targets.shape[0] == 1): return targets.iloc[0]['ATOM']
        raise Exception('Please check your data or change the criterium.') 