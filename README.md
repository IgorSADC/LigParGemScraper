# LigParGemScraper
A small piece of code for programatically access ligParGem

The LigParGemScraper works as a State Machine. The reason for that is to allow users to override an already declared scraper to make a new one without creating a new object.
Everytime you set de smiles or the pdb the code is going to reset the other option. You need to choose between a smiles or pdb as it works on the site.

#Make Connection
After finishing your configuration you should call the function make_connection. The get_values function only work after this.

#Get_values
The get_values function expects a list with every file that you want to retrieve from the WebSite. An example is ['PDB', 'GRO']!

# File Path
The default file path is used to find the pdb and to save the files you pass to the 'get_values()' function. The files are always saved on the same folder of the PDB in case you use the PDB.
If you are using the smiles the code can actually create the folder for you if it doen't exist!
