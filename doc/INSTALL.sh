#!/bin/bash -e
if [ ! -f "release-noresm2.1.3.tar.gz" ];then
    wget -c https://github.com/NorESMhub/NorESM/archive/refs/tags/release-noresm2.1.3.tar.gz
fi
if [ ! -d "NorESM-release-noresm2.1.3" ];then
    tar xzf release-noresm2.1.3.tar.gz
fi
cd NorESM-release-noresm2.1.3/

## necessary patch
### Externals.cfg, include norcpm as a component
if [ -z "$(grep norcpm Externals.cfg)" ]; then
    cat >> Externals.cfg << EOF
[norcpm]
tag = n2.0.8.4
protocol = git
repo_url = https://github.com/NorESMhub/NorCPM_ESP
local_path = components/norcpm
required = true
EOF
fi

## download components 
manage_externals/checkout_externals

### cime/config/cesm/config_files.xml
if [ -z "$(grep norcpm cime/config/cesm/config_files.xml)" ]; then
    sed -i -e'219i<value component="norcpm"    >$SRCROOT/components/norcpm</value>' cime/config/cesm/config_files.xml
fi
### cime_config/config_compsets.xml
if [ -z "$(grep norcpm cime_config/config_compsets.xml)" ]; then
    sed -i -e'559i<!-- NorCPM as esp HIST compset. Add this compset -->\
  <compset>\
    <alias>NHISTNCPM</alias>\
    <lname>HIST_CAM60%NORESM_CLM50%BGC-CROP_CICE%NORESM-CMIP6_BLOM%ECO_MOSART_SGLC_SWAV_NORCPM_BGC%BDRDDMS</lname>\
  </compset>\
  ' \
    cime_config/config_compsets.xml
fi

### patch for BLOM multi instances
if [ -z "$(grep '        pg_blom = OcnInParamGen.from_namelist_xml(xml_fil)' components/blom//cime_config/buildnml)" ]; then
    sed -i -e'116i\        pg_blom = OcnInParamGen.from_namelist_xml(xml_fil)' components/blom//cime_config/buildnml
fi
if [ -z "$(grep NTASKS_PER_INST_OCN components/blom//cime_config/buildcpp)" ];then
    sed -i -e's/NTASKS_OCN/NTASKS_PER_INST_OCN/' components/blom//cime_config/buildcpp
fi

