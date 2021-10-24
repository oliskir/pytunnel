import json
import argparse


def parse_json(path):
    fil = open(path, "r")
    data = json.load(fil)
    fil.close()
    return data


def parse_cla():
    parser = argparse.ArgumentParser(description='Numerical simulation of quantum mechanical tunneling')
    parser.add_argument('-c', '--config_file', type=str, help='json configuration file', default='config.json')
    parser.add_argument('-e', '--energy', type=str, help='kinetic energy of incident particle', default=None)
    args = parser.parse_args()
    return args
