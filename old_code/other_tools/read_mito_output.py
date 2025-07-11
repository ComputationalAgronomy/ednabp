import pandas as pd

xls = pd.ExcelFile("MiFish_test.xlsx")

spc_smpdata_df = pd.read_excel(xls, "List of Sample Details")
spc_smpdata_df["Species"] = spc_smpdata_df["Species"].ffill()

hap_id = spc_smpdata_df["Haploid ID"]
seq = spc_smpdata_df["Sequence"]
seq_size = spc_smpdata_df["Size"]
spc = spc_smpdata_df["Species"]
hap_seq = dict(zip(hap_id, seq, strict=False))
hap_size = dict(zip(hap_id, seq_size, strict=False))
hap2level = {
    key: {"species": value, "genus": value.split()[0]}
    for key, value in zip(hap_id, spc, strict=False)
}

spc_metadata_df = pd.read_excel(xls, "Comparison of Samples")
spc2level = (
    spc_metadata_df.filter(
        items=["Class", "Order", "Family", "Scientific Name"]
    )
    .set_index("Scientific Name")
    .T.to_dict("dict")
)
for hap, level in hap2level.copy().items():
    hap2level[hap].update(spc2level[level["species"]])

spc_metadata_df = spc_metadata_df.dropna(how="all")
spc_metadata_df = spc_metadata_df[
    spc_metadata_df["Class"] != "Followings are non-fish species"
]

spc_metadata_df = (
    spc_metadata_df.filter(
        items=["Class", "Order", "Family", "Scientific Name"]
    )
    .set_index("Scientific Name")
    .T
)

spc_info = (
    spc_metadata_df.filter(
        items=[
            "Scientific Name",
            "Water area",
            "Habitat",
            "DepthS",
            "DepthD",
            "IUCN Red List Status",
            "Importance in Fisheries",
        ]
    )
    .set_index("Scientific Name")
    .T.to_dict("dict")
)
