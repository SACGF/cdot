from bioutils.assemblies import make_ac_name_map, make_name_ac_map

# Generated via:
# import pandas as pd
# assembly_report = "./GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt"
# names = ["Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name"]
# df = pd.read_csv(assembly_report, comment='#', sep='\t', header=None, names=names)
# equals_mask = df["Relationship"] == '='  # MT not there
# df = df[equals_mask]
# ac_name_map = dict(df[["RefSeq-Accn", "Sequence-Name"]].values)

T2T_CHM13v2 = {
    'NC_060925.1': '1',
    'NC_060926.1': '2',
    'NC_060927.1': '3',
    'NC_060928.1': '4',
    'NC_060929.1': '5',
    'NC_060930.1': '6',
    'NC_060931.1': '7',
    'NC_060932.1': '8',
    'NC_060933.1': '9',
    'NC_060934.1': '10',
    'NC_060935.1': '11',
    'NC_060936.1': '12',
    'NC_060937.1': '13',
    'NC_060938.1': '14',
    'NC_060939.1': '15',
    'NC_060940.1': '16',
    'NC_060941.1': '17',
    'NC_060942.1': '18',
    'NC_060943.1': '19',
    'NC_060944.1': '20',
    'NC_060945.1': '21',
    'NC_060946.1': '22',
    'NC_060947.1': 'X',
    'NC_060948.1': 'Y'
}

def get_ac_name_map(assembly_name):
    if assembly_name == "CHM13v2.0":
        return T2T_CHM13v2
    return make_ac_name_map(assembly_name)


def get_name_ac_map(assembly_name):
    if assembly_name == "CHM13v2.0":
        return {name: ac for ac, name in T2T_CHM13v2.items()}
    return make_name_ac_map(assembly_name)
