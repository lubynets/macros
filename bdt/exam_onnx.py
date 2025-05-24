import argparse

import onnx

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--input-file', dest='input_file', help='input file in .onnx format.', default='')
args = parser.parse_args()
input_file = args.input_file

# Load the ONNX model
model = onnx.load(input_file)

feature_names = None
for kv in model.metadata_props:
    if kv.key == "feature_names":
        feature_names = kv.value.split(",")   # recover the list
        break

if feature_names:
    print("Feature names in order:")
    for i, name in enumerate(feature_names, 1):
        print(f"{i:2d}. {name}")
else:
    print("No feature_names metadata found.")
