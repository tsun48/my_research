#!/usr/bin/env python

import os
import sys
import re

#---------------------------------------------------------
import csv
def read_csv_to_dict_of_lists(filename):
    with open(filename, 'r') as file:
        reader = csv.DictReader(file)
        result = {}
        for column in reader.fieldnames:
            result[column] = []
        for row in reader:
            for column in reader.fieldnames:
                result[column].append(row[column])
    return result
#---------------------------------------------------------

# data = read_csv_to_dict_of_lists(quast_report_csv)
# print(data)

#-----------------------------------------------------------------------------------------------
host_map = {
    "africanwilddonkey": "other",
    "commondolphin": "other",
    "cowwithbovinemastitis": "other",
    "fallowdeer": "other",
    "hedgehog": "pig",
    "iberianibex": "other",
    "milk": "milk",
    "monkey": "other",
    "plant": "other",
    "turkeyvulture": "poultry",
    "bird": "poultry",
    "duck": "poultry",
    "environment": "other",
    "goat": "other",
    "bovine": "bovine",
    "sheep": "other",
    "cow": "other",
    "wildturkey": "poultry",
    "rodent": "other",
    "cat": "cat",
    "dog": "dog",
    "animal": "other",
    "chicken": "poultry",
    "cattle": "cattle",
    "horse": "horse",
    "bostaurus": "other",
    "notapplicable": "other",
    "domesticpig": "pig",
    "pig": "pig",
    "wildpig": "pig",
    "human": "human",
}
def map_Host_common_name_to_HostCategoty(Host_common_name):
    # print(f"[map_Host_common_name_to_HostCategoty] Host_common_name= {Host_common_name}")
    res = ""
    if Host_common_name in host_map:
        res = host_map[Host_common_name]
    elif re.search("human", Host_common_name, re.IGNORECASE):
        res = "human"
    elif re.search("wild", Host_common_name, re.IGNORECASE):
        res = "wild"
    elif re.search("environment", Host_common_name, re.IGNORECASE):
        res = "environment"
    else:
        res = "misc"
    return res
#-----------------------------------------------------------------------------------------------
