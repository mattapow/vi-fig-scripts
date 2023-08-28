import os
import dendropy
import numpy as np
import pandas as pd

# Function to extract unique bit strings from a tree list
def extract_bit_strings(treelist):
    bit_strings = set()
    for tree in treelist:
        tree.encode_bipartitions()
        for edge in tree.edges():
            bit_strings.add(edge.split_as_bitstring())
    return bit_strings

# Function to calculate the mean edge lengths for each unique bit string in a tree list
def calculate_mean_edge_lengths(treelist, bit_strings):
    edge_lengths_dict = {bit_string: [] for bit_string in bit_strings}
    for tree in treelist:
        for edge in tree.edges():
            if edge.length is not None:
                bipartition = edge.split_as_bitstring()
                if bipartition in bit_strings:
                    edge_lengths_dict[bipartition].append(edge.length)
    
    # Calculate the mean edge length for each unique bit string
    mean_edge_lengths = {
        bit_string: np.mean(edge_lengths) if edge_lengths else None
        for bit_string, edge_lengths in edge_lengths_dict.items()
    }
    std_edge_lengths = {
        bit_string: np.std(edge_lengths) if edge_lengths else None
        for bit_string, edge_lengths in edge_lengths_dict.items()
    }
    n_edge_lengths = {
        bit_string: np.count_nonzero(edge_lengths) if edge_lengths else None
        for bit_string, edge_lengths in edge_lengths_dict.items()
    }
    return mean_edge_lengths, std_edge_lengths, n_edge_lengths


def main():
    ds=1
    ds_path = os.path.join("..", "analysis", f"ds{ds}")

    path_dodona = os.path.join(ds_path, "vi", "up_nj", "boosts", "d3_lr1_i3_b3_boosts")
    tree_files = [
        os.path.join(ds_path, "mb", "gtr", 'gtr.run1.t'),
        os.path.join(ds_path, "mb", "gtr", 'gtr.run2.t'),
        os.path.join(path_dodona, "samples.t")
    ]
    treelists_names = ["MB1", "MB2", "Dodonaphy"]
    treelists = []
    for tree_file in tree_files:
        print(f"reading {tree_file}")
        treelist = dendropy.TreeList.get(path=tree_file, schema='nexus')
        treelists.append(treelist)

    # Get a unified taxon namespace
    taxon_namespace = dendropy.TaxonNamespace(treelists[0][0].taxon_namespace)
    for treelist in treelists:
        for tree in treelist:
            tree.retain_taxa_with_labels(labels=taxon_namespace.labels())

    # Extract unique bit strings from all tree lists
    print("Extracting bit_strings")
    all_bit_strings = set()
    for treelist in treelists:
        all_bit_strings.update(extract_bit_strings(treelist))

    # Calculate mean edge lengths for each tree list and store in a dictionary
    print("calculating mean + std of edge lengths")
    mean_edge_lengths_by_treelist = {}
    std_edge_lengths_by_treelist = {}
    count_edge_lengths_by_treelist = {}
    for i, treelist in enumerate(treelists):
        mean_edge_lengths_treelist, std_edge_lengths_treelist, count_edge_lengths_treelist= calculate_mean_edge_lengths(treelist, all_bit_strings)
        mean_edge_lengths_by_treelist[treelists_names[i]] = mean_edge_lengths_treelist
        std_edge_lengths_by_treelist[treelists_names[i]] = std_edge_lengths_treelist
        count_edge_lengths_by_treelist[treelists_names[i]] = count_edge_lengths_treelist

    # Save the results to a CSV file
    data = []
    for bit_string in all_bit_strings:
        row = [bit_string]
        for i in range(len(treelists)):
            mean_length = mean_edge_lengths_by_treelist[treelists_names[i]].get(bit_string, np.nan)
            row.append(mean_length)
            std_length = std_edge_lengths_by_treelist[treelists_names[i]].get(bit_string, np.nan)
            row.append(std_length)
            count_length = count_edge_lengths_by_treelist[treelists_names[i]].get(bit_string, np.nan)
            row.append(count_length)
        data.append(row)

    columns = ['Bit_String']
    for i in range(len(treelists_names)):
        columns.append(f'Mean_Edge_Length_{treelists_names[i]}')
        columns.append(f'Std_Edge_Length_{treelists_names[i]}')
        columns.append(f'Count_Edge_Length_{treelists_names[i]}')

    df = pd.DataFrame(data, columns=columns)
    sorted_df = df.dropna().append(df[df.isna().any(axis=1)])
    sorted_df.to_csv('tmp/branch_lengths_by_treelist_sorted.csv', index=False)


if __name__ == "__main__":
    main()
