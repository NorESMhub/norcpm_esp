## Before checkout components:
Externals.cfg:
```
[norcpm]
tag = n2.0.8.3
protocol = git
repo_url = https://github.com/NorESMhub/NorCPM_ESP
local_path = components/norcpm
required = true
```

## Add NorCPM settings in NorESM xml file
### Add NorCPM as an ESP component
NorESM/cime/config/cesm/config_files.xml:
```
  <entry id="COMP_ROOT_DIR_ESP">
    <type>char</type>
    <default_value>unset</default_value>
    <values>
      <value component="norcpm"    >$SRCROOT/components/norcpm</value> !! <--- add this line
      <value component="desp"      >$CIMEROOT/src/components/data_comps/desp</value>
      <value component="sesp"      >$CIMEROOT/src/components/stub_comps/sesp</value>
    </values>
    <group>case_comps</group>
    <file>env_case.xml</file>
    <desc>Root directory of the case external system processing (esp) component  </desc>
    <schema>$CIMEROOT/config/xml_schemas/config_compsets.xsd</schema>
  </entry>
```
### Add a compset conatains NorCPM
cime_config/config_compsets.xml:
```
  <!-- NorCPM as esp HIST compset. Add this compset -->
  <compset>
    <alias>NHISTNCPM</alias>
    <lname>HIST_CAM60%NORESM_CLM50%BGC-CROP_CICE%NORESM-CMIP6_BLOM%ECO_MOSART_SGLC_SWAV_NORCPM_BGC%BDRDDMS</lname>
  </compset>
```
