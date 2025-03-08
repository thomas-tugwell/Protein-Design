# Protein Design
A collection of scripts to direct and help the flow of my design pipeline and to evaluate and filter the results.
The purpose of this repository is to enhance compatibility between popular design tools such as RFDiffusion, LigandMPNN, Alphafold, Chai, ESMfold, etc... as well as introduce new tools to aid others in their designs.

<img width="519" alt="image" src="https://github.com/user-attachments/assets/39f1a01a-fce5-478d-bed2-08798007799f" />

# Usage
--Under Construction--

## Docking with space_dock.py

<img width="1069" alt="image" src="https://github.com/user-attachments/assets/ea286793-bc28-40f5-884f-7a23d226c645" />

Need to install Bio module:
```
pip install Bio
```

You need to edit these lines to supply your structures as separate pds.
```
pdb_mobile  = "3_9_chainsfrom1_noend.pdb"
pdb_station = "pent_pore_chains_noend.pdb"
```


Modify these lines to match your input chains (both inputs have same chains) and select one residue from each pore to define the planes that will be aligned.
In a future version, the code will be expanded to allow mismatched length oligomers and chain identities that are different between inputs. For now, both inputs should contain the same chains. 

```
chains = ["A", "B", "C", "D", "E"]
res1, res2 = 75, 109  # Mobile (res 75), Stationary (res 109)
```

The input pdbs should begin renumbering starting at each new chain:

```
ATOM      1  N   MET A   1     -12.015 -12.546  -7.489  1.00 59.29           N  
ATOM      2  CA  MET A   1     -13.205 -11.866  -6.985  1.00 60.33           C  
ATOM      3  C   MET A   1     -13.804 -12.621  -5.804  1.00 62.80           C  
ATOM      4  CB  MET A   1     -14.247 -11.714  -8.094  1.00 49.42           C  
ATOM      5  O   MET A   1     -13.944 -13.844  -5.851  1.00 56.70           O  
ATOM      6  CG  MET A   1     -15.396 -10.787  -7.731  1.00 44.03           C  
ATOM      7  SD  MET A   1     -16.614 -10.610  -9.093  1.00 54.10           S  
ATOM      8  CE  MET A   1     -17.911 -11.738  -8.511  1.00 37.02           C  
ATOM      9  N   ARG A   2     -14.298 -11.856  -4.636  1.00 77.00           N  
ATOM     10  CA  ARG A   2     -14.804 -12.546  -3.454  1.00 78.72           C  
ATOM     11  C   ARG A   2     -16.328 -12.500  -3.403  1.00 80.20           C  
ATOM     12  CB  ARG A   2     -14.217 -11.933  -2.181  1.00 69.49           C
.
.
.
ATOM    699  CB  ALA A  96     -39.147 -18.831   6.950  1.00 78.12           C  
ATOM    700  O   ALA A  96     -38.086 -21.853   6.636  1.00 75.66           O  
TER     700      ALA A  96                                                      
ATOM    701  N   MET B   1     -25.006 -19.820   6.224  1.00 43.99           N  
ATOM    702  CA  MET B   1     -25.785 -18.862   7.001  1.00 49.89           C  
ATOM    703  C   MET B   1     -26.451 -17.835   6.091  1.00 46.78           C  
ATOM    704  CB  MET B   1     -24.898 -18.154   8.027  1.00 37.77           C  
ATOM    705  O   MET B   1     -25.769 -17.088   5.388  1.00 43.70           O  
ATOM    706  CG  MET B   1     -25.312 -18.399   9.469  1.00 35.57           C  
ATOM    707  SD  MET B   1     -24.333 -17.406  10.661  1.00 25.05           S  
ATOM    708  CE  MET B   1     -25.161 -17.859  12.211  1.00 27.55           C  
```

Once the planes are calculated and aligned, these lines define the distance along the normal vector between the midpoints of each pore.
```
n_steps = 20
offset_values = np.linspace(20.0, -20.0, n_steps)
```
