"""Converts a COBRA .mat model to SBML. Used in place of MATLAB's writeCbModel
SBML export, whose OutputSBML MEX binary has no Apple Silicon build.

Usage: python3 mat_to_sbml.py <model.mat> <model.xml>
"""
import sys, warnings
warnings.filterwarnings('ignore')
import cobra.io

mat_in, xml_out = sys.argv[1], sys.argv[2]
model = cobra.io.load_matlab_model(mat_in)

# MATLAB's SBML writer defaults every gene product's sboTerm to SBO:0000243;
# write_sbml_model does not.
for gene in model.genes:
    gene.annotation.setdefault('sbo', 'SBO:0000243')

cobra.io.write_sbml_model(model, xml_out)
