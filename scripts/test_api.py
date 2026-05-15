import requests, json

# Get Drug fields
r = requests.post(
    "https://api.platform.opentargets.org/api/v4/graphql",
    json={"query": '{ __type(name: "Drug") { fields { name } } }'}
)
print("Drug fields:")
for f in r.json()["data"]["__type"]["fields"]:
    print(" ", f["name"])

# Get ClinicalDiseaseListItem fields
r = requests.post(
    "https://api.platform.opentargets.org/api/v4/graphql",
    json={"query": '{ __type(name: "ClinicalDiseaseListItem") { fields { name } } }'}
)
print("\nClinicalDiseaseListItem fields:")
for f in r.json()["data"]["__type"]["fields"]:
    print(" ", f["name"])
