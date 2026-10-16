#!/usr/bin/env python

import os
import sys
import re
import pandas as pd

# metadata_csv = "../ccir_data/NCBI_ST398_metadata_7flds_no_NA.csv"

metadata_csv = sys.argv[1]

# data = pd.read_csv(metadata_csv, dtype=str)
# # print(data["Location"].value_counts())
# country = data["Location"].copy()
#
# def consolidate(string, data, replacewith):
#     for i in range(len(data)):
#         if string in str(data[i]):
#             country.replace(data[i], replacewith, inplace=True)

location_to_country = {
    "Germany": "Germany",
    "China": "China",
    "Denmark": "Denmark",
    "Switzerland": "Switzerland",
    "Portugal": "Portugal",
    "Japan": "Japan",
    "Thailand": "Thailand",
    "United Kingdom": "United_Kingdom",
    "New Zealand": "New_Zealand",
    "Ecuador": "Ecuador",
    "Mexico": "Mexico",
    "Africa": "South_Africa",
    "Dominican Republic": "Dominican_Republic",
    "Australia": "Australia",
    "Korea": "South_Korea",
    "Brazil": "Brazil",
    "Argentina": "Argentina",
    "Russia": "Russia",
    "Canada": "Canada",
    "Cape Verde": "Cape_Verde",
    "Colombia": "Colombia",
    "Sweden": "Sweden",
    "USA": "USA",
    "France": "France",
    "Spain": "Spain",
}
# country = country.str.replace(" ", "_")
# print(country.value_counts())
# data["Country"] = country

country_to_ccode = {
    "Germany": "DEU",
    "China": "CHN",
    "France": "FRA",
    "Netherlands": "NLD",
    "Denmark": "DNK",
    "Switzerland": "CHE",
    "Spain": "ESP",
    "Portugal": "PRT",
    "Italy": "ITA",
    "Belgium": "BEL",
    "Hong_Kong": "HKG",
    "Hong Kong": "HKG",
    "HongKong": "HKG",
    "United_Kingdom": "GBR",
    "UnitedKingdom": "GBR",
    "Japan": "JPN",
    "Slovenia": "SVN",
    "Brazil": "BRA",
    "Thailand": "THA",
    "Australia": "AUS",
    "South_Korea": "KOR",
    "SouthKorea": "KOR",
    "Taiwan": "TWN",
    "Finland": "FIN",
    "New_Zealand": "NZL",
    "NewZealand": "NZL",
    "Ecuador": "ECU",
    "Russia": "RUS",
    "Poland": "POL",
    "Mexico": "MEX",
    "Hungary": "HUN",
    "South_Africa": "ZAF",
    "SouthAfrica": "ZAF",
    "Dominican_Republic": "DOM",
    "DominicanRepublic": "DOM",
    "Canada": "CAN",
    "Argentina": "ARG",
    "Algeria": "DZA",
    "Austria": "AUT",
    "Ireland": "IRL",
    "Dominica": "DMA",
    "Martinique": "MTQ",
    "Czech_Republic": "CZE",
    "CzechRepublic": "CZE",
    "Colombia": "COL",
    "Paraguay": "PRY",
    "Singapore": "SGP",
    "Sweden": "SWE",
    "Trinidad_and_Tobago": "TTO",
    "TrinidadandTobago": "TTO",
    "Turkey": "TUR",
    "Cape_Verde": "CPV",
    "CapeVerde": "CPV",
    "Luxembourg": "LUX",
    "Greece": "GRC",
    "Hong_Kong": "HKG",
    "Greenland": "GRL",
    "Czech Republic": "CZE",
}
# data["Country_Code"] = cCode
# data.to_csv("_added_columns.csv", index=False)

host_to_common = {
    "Homo_sapiens": "Human",
    "H_sapiens": "Human",
    "Sus_scrofa_domesticus": "Domestic_pig",
    "Sus_scrofa": "Wild_pig",
    "Sus": "pig",
    "Swine": "pig",
    "Equus_ferus_caballus": "Horse",
    "Equus_caballus": "Horse",
    "Equus_asinus_africanus": "African_wild_donkey",
    "Gallus_gallus_domesticus": "Chicken",
    "Gallus_gallus": "Chicken",
    "Meleagris_gallopavo": "Wild_turkey",
    "Canis_lupus_familiaris": "Dog",
    "Felis_catus": "Cat",
    "Rattus": "Rodent",
    "Rodentia": "Rodent",
    "Rat": "Rodent",
    "Ovis_aries": "Sheep",
    "Capra_aegagrus_hircus": "Goat",
    "Capra_hircus": "Goat",
    "Capra_pyrenaica": "Iberian_ibex",
    "Anas_platyrhynchos_domesticus": "Domestic_duck",
    "Cervus_elaphus": "Red_deer",
    "Gyps_fulvus": "Eurasian_griffon_vulture",
    "Delphinus_delphis": "Common dolphin",
    "Dama_dama": "Fallow_deer",
    "Cathartes_aura": "Turkey_vulture",
    "Erinaceus_europaeus": "Hedgehog",
    "environment_sample": "Environment",
}

def str_map(str_, map_):
    for k, v in map_.items():
        if k.lower() in str_.lower():
            return re.sub(" ", "", re.sub("_", "", v.lower()))
    return re.sub(" ", "", re.sub("_", "", str_.lower()))

file = open(metadata_csv, "r")
lines = file.readlines()
hdr_flds = []
out_flds = ['Accession_number', 'ST', 'ST_pattern', 'Biosample_ID',
            'Host', 'Collection_date', 'Location',
            'Comments', 'Country', 'Country_Code', 'Collection_year',
            'Host_common_name', 'Linkage_name']
for line in lines:
    line = re.sub("\s*$", "", line)
    # print(f"line= {line}")
    vec = line.split(',')
    if len(hdr_flds) == 0:
        hdr_flds = vec
        # print(f"hdr_flds= {hdr_flds}")
        print(','.join(out_flds))
    else:
        rec = {hdr_flds[i]: re.sub(' ', '_', vec[i]) for i in range(len(hdr_flds))}
        # print(f"rec= {rec}")
        rec['Accession_number'] = rec['NCBI_assembly_ID']
        rec['Comments'] = ""
        rec['Country'] = str_map(rec['Location'], location_to_country)
        rec['Country_Code'] = str_map(rec['Country'], country_to_ccode).upper()
        rec['Collection_year'] = rec['Collection_date'][:4]
        rec['Host_common_name'] = str_map(rec['Host'], host_to_common).replace(" ", "_")
        rec['Linkage_name'] = f'{rec["Biosample_ID"]}_{rec["Host_common_name"]}_{rec["Country_Code"]}_{rec["Collection_year"]}'

        out_vec = [rec[key] for key in out_flds]
        print(','.join(out_vec))

'''

PYTHONPATH=. metadata_add_columns.py \
 NCBI_ST398_metadata_7flds_no_NA.csv \
 | grep -f tsun_accessions.txt > NCBI_CC398_metadata_v3_Oct16.csv

'''
