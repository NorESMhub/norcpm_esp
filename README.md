# NorCPM_esp

# NorCPM Structure

### Scientific description:
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1169902.svg)](https://doi.org/10.5281/zenodo.1169902)

### Overview
  NorCPM do ensemble data assimilation with NorESM1 and 2. NorCPM_ESP apply the process as an ESP component in NorESM2. Which modify model data in memory to reduce overhead time from restarting model.

  For now NorCPM_ESP is only DA for BLOM component.
  
  It is still under development.

### Usage
  The content in cime_config/config_compsets.xml is needed to insert to NorESM/cime_config/config_compsets.xml

  Check the test script in 'testscript/'.

### Portability
  It only tested on Betzy.

  The SST DA needs monthly file named as YYYY_MM_DD.nc and placed at specific structure of directory.

