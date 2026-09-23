# Data from: Drivers of coastal benthic communities in a complex environmental setting

---

This README.md was generated on 2024-06-03 by Yuting Vicky Lin ([vicky.linyuting@gmail.com](mailto:vicky.linyuting@gmail.com)).

1.  **Author Information**
    *   First author
        *   Name: Dr. Yuting Vicky Lin
        *   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
    *   Second author
        *   Name: Prof. Pierre-Alexandre Château
        *   Institution: Department of Marine Environment and Engineering, National Sun Yat-Sen University, Kaohsiung 80420, Taiwan
    *   Third author
        *   Name: Prof. Yoko Nozawa
        *   Institution: Tropical Biosphere Research Center, University of the Ryukyus, Okinawa 905-0227, Japan
            Biodiversity Research Center, Academia Sinica, Taipei 11529, Taiwan
            Department of Marine Science, Faculty of Fisheries and Marine Science, Universitas Diponegoro, Semarang 50275, Indonesia
    *   Fourth author
        *   Name: Prof. Chih-Lin Wei
        *   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
    *   Fifth author
        *   Name: Dr. Rainer Ferdinand Wunderlich
        *   Institution: Department of Bioenvironmental Systems Engineering, National Taiwan University, Taipei 10617, Taiwan &
            5INRAE, UR EABX, 33612 Cestas, France
    *   Sixth & corresponding author
        *   Name: Prof. Vianney Denis
        *   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
        *   Email: [vianneydenis@ntu.edu.tw](mailto:vianneydenis@ntu.edu.tw)
2.  **Date of data collection**: 2016-2020
3.  **Geographic location of data collection**: Taiwan, West Pacific
4.  **Funding sources that supported the collection of the data**: Marine National Park of Taiwan, National Science and Technology Council of Taiwan, Ocean Affairs Council of Taiwan
5.  **Recommended citation for this dataset**
    Lin YT, Château P-A, Nozawa Y, Wei C-L, Wunderlich RF, Denis V (2024), Data from: Drivers of coastal benthic communities in a complex environmental setting, Dryad, Dataset, [https://doi.org/10.5061/dryad.j0zpc86nz](https://doi.org/10.5061/dryad.j0zpc86nz)

\[Data and R script are also available through the GitHub repository [https://github.com/NTU-FRELab/drivers-communities](https://github.com/NTU-FRELab/drivers-communities) in order to replicate our analyses]

## DESCRIPTION OF THE DATA AND FILE STRUCTURE

### DATA & FILE OVERVIEW

1.  **Description of dataset**

    These data were used to identify distinct benthic communities around Taiwan and the important environmental drivers causing the difference.
2.  **File List**

    File 1: Linetal_dataset_Benthic.csv: benthic cover estimates of the 433 sampling transects
    File 2: Linetal_dataset_Env.csv: data of 21 environmental drivers in the 433 sampling transects
    File 3: Linetal_dataset_Cor.csv: coordination and benthic information used to generate Figure 1b
    File 4: Linetal_dataset_TW.shp, .cpg, .dbf, .prj, .shx: a shapefile of Taiwan map used to generate Figure 1a

### METHODOLOGICAL INFORMATION

A detailed description of data acquisition and processing can be found in the published manuscript in the Marine Pollution Bulletin ([https://doi.org/10.1016/j.marpolbul.2024.116462](https://doi.org/10.1016/j.marpolbul.2024.116462)).

#### DATA-SPECIFIC INFORMATION

##### **Linetal\_dataset\_Benthic.csv**

1.  Number of variables/columns: 28
2.  Number of cases/rows: 434
3.  Missing data codes
    None
4.  Variable List
    Column A - transect
    Column B - ac_erect_single
    Column C - ag_articulated_calcareous
    Column D - ag_corticated_foliose
    Column E – ag_corticated_macrophyte
    Column F – cca_crustose
    Column G – cy_filamentous
    Column H – hc_arborescent
    Column I – hc_bushy
    Column J – hc_column
    Column K – hc_encrusting
    Column L – hc_foliose
    Column M – hc_massive
    Column N – hc_table
    Column O – oc_arborescent
    Column P – oc_bushy
    Column Q – oc_cluster
    Column R – oc_digitate
    Column S – oc_encrusting
    Column T – oc_lobate
    Column U – oc_massive
    Column V – sg
    Column W – sp_encrusting
    Column X – sp_massive
    Column Y – sp_repent
    Column Z – turf_filamentous
    Column AA – zo_encrusting
    Column AB – zo_massive
5.  Descriptions
    The values in this .csv file are benthic covers estimated based on the relative probability of the top five identifications suggested by a trained deep neural network in CoralNet (Global ID: 22264). The descriptions of each variables are as below:
    ac_erect_single: an ascidian individual that is in the erect morphology; ac denotes to ascidian
    ag_articulated_calcareous: Articulated calcareous alga; ag denotes to alga
    ag_corticated_foliose: Corticated foliose alga
    ag_corticated_macrophyte: Corticated macrophyte
    cca_crustose: Crustose coralline alga; cca denotes to crustose coralline alga
    cy_filamentous: Filamentous cyanobacterium; cy denotes to cyanobacteria
    hc_arborescent: A stony coral colony with arborescent morphology; hc denotes to stony coral
    hc_bushy: A stony coral colony with bushy morphology
    hc_column: A stony coral colony with columnar morphology
    hc_encrusting: A stony coral colony with encrusting morphology
    hc_foliose: A stony coral colony with foliose morphology
    hc_massive: A stony coral colony with massive morphology
    hc_table: A stony coral colony with table morphology
    oc_arborescent: A octocoral colony with arborescent morphology; oc denotes to octocoral
    oc_bushy: A octocoral colony with bushy morphology
    oc_cluster: A octocoral colony with cluster morphology
    oc_digitate: A octocoral colony with digitate morphology
    oc_encrusting: A octocoral colony with encrusting morphology
    oc_lobate: A octocoral colony with lobate morphology
    oc_massive: A octocoral colony with massive morphology
    sg: Seagrass; sg denotes to seagrass
    sp_encrusting: A sponge colony with encrusting morphology; sp denotes to sponge
    sp_massive: A sponge colony with massive morphology
    sp_repent: A sponge colony with repent morphology
    turf_filamentous: Filamentous turf alga; turf denotes to turf alga
    zo_encrusting: A zoanthid colony with encrusting morphology; zo denotes to zoanthid
    zo_massive: A zoanthid colony with massive morphology

##### **Linetal\_dataset\_Env.csv**

1.  Number of variables/columns: 26
2.  Number of cases/rows: 434
3.  Missing data codes
    None
4.  Variable List
    *   Column A - Region
    *   Column B - Site
    *   Column C - Depth
    *   Column D - Latitude
    *   Column E – Longitude
    *   Column F – Mean_SST
    *   Column G – SD_SST
    *   Column H – Light_intensity
    *   Column I – Wave_exposure
    *   Column J – Wave_height
    *   Column K – Mean_chl_a
    *   Column L – SD_chl_a
    *   Column M – Nitrate
    *   Column N – Nitrite
    *   Column O – Phosphate
    *   Column P – DHW
    *   Column Q – DHW_recovery
    *   Column R – Typhoon_disturbance
    *   Column S – Typhoon_recovery
    *   Column T – Typhoon_frequency
    *   Column U – Anthropogenic_land_use
    *   Column V – Forest_land_use
    *   Column W – Population_density
    *   Column X – Tourist_visitors
    *   Column Y – Unstable_substrate_cover
    *   Column Z – Management_status
5.  Descriptions
    The environmental data in 433 benthic sampling transects around the coastline of Taiwan. The descriptions of each variables are as below:
    Region: Region in Taiwan where the benthic sampling was conducted
    Site: Name of the sampling site
    Depth: Depth where the benthic sampling was conducted. Unit: meter
    Latitude: Latitude where the benthic sampling was conducted
    Longitude: Longitude where the benthic sampling was conducted
    Mean_SST: Mean sea surface temperature extracted from NOAA Global Coral Bleaching Monitoring. Unit: °C. Mean sea surface temperature is abbreviated to Mean_SST
    SD_SST: Standard deviation of sea surface temperature extracted from NOAA Global Coral Bleaching Monitoring. Unit: °C. Standard deviation of sea surface temperature is abbreviated to SD_SST
    Light_intensity: Light intensity calculated from Aqua MODIS, NPP, L3SMI, Global, 4km. Unit: E x m-2 x year-1
    Wave_exposure: Wave exposure extracted from Yeager et al. (2017). Unit: kW x m-1
    Wave_height: Wave height extracted from Open Whether Data, Central Weather Bureau, Taiwan. Unit: meter
    Mean_chl_a: Mean chlorophyll a extracted from Global Ocean Chlorophyll, PP and PFT from Satellite Observations: Monthly and Daily Interpolated. Unit: mg x m-3. Mean chlorophyll a is abbreviated to Mean_chl_a
    SD_chl_a: Standard deviation of chlorophyll a extracted from Global Ocean Chlorophyll, PP and PFT from Satellite Observations: Monthly and Daily Interpolated. Unit: mg x m-3. Standard deviation of chlorophyll a is abbreviated to SD_chl_a
    Nitrate: Nitrate concentration data downloaded from iOcean, Ocean Conservation Administration, Taiwan. Dongsha's data was provided by the laboratory of Prof. Fuh-Kwo Shiah in Academia Sinica, Taipei, Taiwan. Unit: μm
    Nitrite: Nitrite concentration data downloaded from iOcean, Ocean Conservation Administration, Taiwan. Dongsha's data was provided by the laboratory of Prof. Fuh-Kwo Shiah in Academia Sinica, Taipei, Taiwan. Unit: μm
    Phosphate: Phosphate concentration data downloaded from iOcean, Ocean Conservation Administration, Taiwan. Dongsha's data was provided by the laboratory of Prof. Fuh-Kwo Shiah in Academia Sinica, Taipei, Taiwan. Unit: μm
    DHW: Degree heating week extracted from NOAA Global Coral Bleaching Monitoring. Unit: °C-weeks. Degree heating week is abbreviated to DHW
    DHW_recovery: Recovery time from the last extreme heat event. The time gap between the year the maximum DHW occurred and the year of sampling, and the data of the maximum DHW were downloaded from NOAA Global Coral Bleaching Monitoring. Unit: year. Recovery time from the last extreme heat event is abbreviated to DHW_recovery
    Typhoon_disturbance: Maximum wind speed for the last severe typhoon, and the wind speed data were from CWB observation data inquire system, Central Weather Bureau, Taiwan. Dongsha's data was provided by Prof. Po-Hsiung Lin in the Department of Atmospheric Science, National Taiwan University, Taipei, Taiwan. Unit: m x s-1. Maximum wind speed for the last severe typhoon is abbreviated as Typhoon disturbance
    Typhoon_recovery: Recovery time from the last severe typhoon. It was calculated by the time gap between the year the nearest severe typhoon occurred and the year of sampling. Data were from CWB observation data inquire system, Central Weather Bureau, Taiwan. Unit: year. Recovery time from the last severe typhoon is abbreviated to Typhoon_recovery
    Typhoon_frequency: Number of severe typhoons in the last decade. Data were from CWB observation data inquire system, Central Weather Bureau, Taiwan. Unit: decade-1. Number of severe typhoons in the last decade is abbreviated to Typhoon frequency
    Anthropogenic_land_use: The area of the agricultural, aquacultural, and urban land use in the township located within a 10 km radius of sampling locations. Data were from MOI Open Data, Ministry of Interior, Taiwan. Unit: km2
    Forest_land_use: The area of the forest land use in the township located within a 10 km radius of sampling locations. Data were from MOI Open Data, Ministry of Interior, Taiwan. Unit: km2
    Population_density: The  average population density within a 10 km radius of sampling locations. Data were from Gridded Population of the World (GPW), v4, 2015, Socioeconomic data and application center, NASA. Unit: people x 10km-radius-1
    Tourist_visitors: The number of visit times in the scenic spots located within a 10 km radius of sampling locations. Data were from Tourism statistics database of the Taiwan Tourism Bureau, Taiwan. Unit: visits
    Unstable_substrate_cover: Covers of unstable subtract (i.e., sand, rubbles) on the canopy layer of reefs. Data were extracted from each of 433 transects. Unit: Percentage (absolute cover).
    Management_status: The marine protected status of sampling locations was classified into “fully no-take,” “open access,” and “restrictions on fishing gear, catch size, species, and/or any recreational activities." Data were from the marine protected area database of Functional Reef Ecology Laboratory, Taiwan

Yeager LA, Marchand P, Gill DA, Baum JK, McPherson JM (2017) Marine Socio-Environmental Covariates: queryable global layers of environmental and anthropogenic variables for marine ecosystem studies. Ecology 98:1976-1976

##### **Linetal\_dataset\_Cor.csv**

1.  Number of variables/columns: 14
2.  Number of cases/rows: 88
3.  Missing data codes
    None
4.  Variable List
    *   Column A - Region
    *   Column B - Site
    *   Column C - Depth
    *   Column D - Latitude
    *   Column E – Longitude
    *   Column F – 1
    *   Column G – 2
    *   Column H – 3
    *   Column I – 4
    *   Column J – 5
    *   Column K – y
    *   Column L – x
    *   Column M – Label1
    *   Column N – Label2
5.  Descriptions
    Coordination, depth and benthic composition of each sampling transect. All these information was used to generate the Figure 1. Abbreviations used for variables are as below:
    1: CCA Community
    2: Turf & Stony Coral Community
    3: Turf Community
    4: Digitate Community
    5: Bushy Community
    y: The vertical position of the pie chart
    x: The horizontal position of the pie chart
    Label1: The number of sampling sites used in Figure 1a
    Label2: The number of sampling sites used in Figure 1b

##### **Linetal\_dataset\_TW\.shp, .cpg, .dbf, .prj, .shx**

1.  Geometry type: MULTIPOLYGON
2.  Spatial boundary: latitude - 20.69750 to 25.63431; longitude - 116.71000 to 122.10850
3.  Missing data codes
    None
4.  Spatial reference: WGS 84
5.  Descriptions
    This shapefile contains spatial polygons within the specified spatial boundary, which is used for generating the map of Taiwan.

## Usage notes

R and QGIS can be used to view Linetal_dataset_TW.shp. Linetal_dataset_TW.cpg, Linetal_dataset_TW.dbf, Linetal_dataset_TW.prj, and Linetal_dataset_TW.shx are the companion files for Linetal_dataset_TW.shp that can not be opened individually but are part of the shapefile data.

Microsoft Excel can be used to view Linetal_dataset_Benthic.csv, Linetal_dataset_Env.csv, and Linetal_dataset_Cor.csv.
