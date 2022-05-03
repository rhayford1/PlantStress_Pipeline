# PlantStress_Pipeline
Plant Stress Gene and Annotation
This repo holds the source code for the pipeline to retrieve plant stress genes and annotation from the scientific literature

Requirements

Programming language – Python

Instructions

1.	Install MongoDB

2.	Create a virtual environment for the scripts 


3.	Run the first script, “load_pmids” which is the main script. This will prompt user for the json input file and the “pmid_extractor’ functions

4.	The “pmid_extractor” script imports the other script, “generate database” which dumps the data in the database.


5.	As the script goes through each pmid entry in the json input file, it saves it to the database if the parsed response has a gene mention and an annotation

