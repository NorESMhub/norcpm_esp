# NorCPM_esp
### Scientific description:
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1169902.svg)](https://doi.org/10.5281/zenodo.1169902)

### Overview
  NorCPM_esp runs as a ESP (External System Processing) component, which is introduced since CESM2. Also in NorESM2

  This component runs after all other componets at every time step, and takes no effect by default.

  In NorCPM_esp, this is used to run data assimilation. Including read observation data, calculation EnKF and feedback to model.

  For the EnKF (Ensemble Kalman Filter), which is necessary to have multiple ensemble members. 

  The multiple instances function is also introduced in CESM2. Which runs multiple copies of model.


### Required settings
  0. There are some setting need to be done in NorESM, please see doc/noresm_settings.md
  1. Compset need be set to NHISTNCPM
  2. NINST must larger than 1. 10 is suggested for scientific meaning.
  3. NTASKS can be default but need to be multiple to NINST. 
     However, it is suggested to use one node per instance (Betzy)
     The NTASKS_ESP need to cover all CPUs.

  Detail setup from NorESM 2.1.3 can be found in doc/INSTALL.sh
  Example run script can be found in doc/norcpm_esp_test.sh

### Usage
  Please start from [NorESM repo](https://github.com/NorESMhub/NorESM). 
  In file Externals.cfg, [norcpm] section, set the required from False to True.
  Then use manage_externals/checkout_externals to download components.
  The test script can be found in NorESM/components/norcpm/doc/create_norcpm_esp.sh.
  Most of the settings are in env_run.xml, start with NORCPM_.

### Portability
  It only tested on Betzy.

  The SST DA needs monthly file named as YYYY_MM_DD.nc and placed at specific structure of directory.

### Issues
  1. Only OISST(daily), EN4 profiles(monthly, temperature and salnity) supported.
  2. Data assimilation(DA) only when middle day of model (43200) instead of day begin (00000) due to the restart time.
  3. EN4 DA takes more than 20 secs. Which is due to parallel method.
