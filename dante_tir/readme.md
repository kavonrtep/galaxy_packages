# Wrapper for the DANTE-TIR tools
- Documentation - https://github.com/kavonrtep/dante_tir
- Installation using conda https://anaconda.org/petrnovak/dante_tir

Galaxy DANTE_TIR tool specification:
 - input DANTE gff23 file
 - input correcondig FASTA file
 - cpu used - based on GALAXY_SLOTS

- output in galaxy - 
  -   
- DANTE_TIR gff3
- DANTE_TIR_final_fasta

# Galaxy xml file specification documentation
https://docs.galaxyproject.org/en/latest/dev/schema.html


# versions - current DANTE_TIR version is 0.2.0

specify version in macros.xml file, if DANTE_TIR version is 0.2.0 then Galaxy tool version is 0.2.0.1
create also .shed.yml for the tools


# dependencies
- installed from conda-forge, bioconda, r, and petrnovak channels
