from cloudmesh.common.FlatDict import FlatDict, read_config_parameters, flatten
from cloudmesh.common.util import readfile

def read_config(filename):
    config_dict = read_config_parameters(filename)
    config = FlatDict(config_dict, sep=".")
    del config["sep"]
    return config

"""
config = read_config("./config.yaml")
print (config)

print(config["covid.ReadMay2022Covid"])

ReadMay2022Covid = config["covid.ReadMay2022Covid"]

print (ReadMay2022Covid)
"""