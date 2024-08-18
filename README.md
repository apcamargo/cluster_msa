# `cluster_msa.py`

This script clusters sequences from a multiple sequence alignment (MSA) into distinct groups, each of which is saved as a separate MSA file. Clustering is performed using HDBSCAN, based on Euclidean distances between the sequences in a high-dimensional embedding space.

This script is adapted from [AF-Cluster](https://github.com/HWaymentSteele/AF_Cluster)[^1] with the following major modifications:
- Uses the [GIANA](https://github.com/s175573/GIANA/blob/master/example_of_CDR3_encoding.ipynb)[^2] amino acid embedding, which approximates BLOSUM62, instead of one-hot encoding.
- Allows using [HDBSCAN](https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html) for clustering instead of DBSCAN, enabling the identification of clusters with variable densities.
- Adds support for compressed input files.
- Removes dependencies and simplifies the code.

## Command-line usage

```sh
$ ./cluster_msa.py -h
```

    usage: cluster_msa.py [-h] [--write-tabular-outputs]
                          [--remove-lowercase-columns]
                          [--max-gap-sequence MAX_GAP_SEQUENCE]
                          [--max-gap-column MAX_GAP_COLUMN] [--use-one-hot-encoding]
                          [--use-hdbscan] [--min-cluster-size MIN_CLUSTER_SIZE]
                          [--min-samples MIN_SAMPLES] [--min-eps MIN_EPS]
                          [--max-eps MAX_EPS] [--eps-step EPS_STEP]
                          [--epsilon-value EPSILON_VALUE]
                          [--cluster-selection-epsilon CLUSTER_SELECTION_EPSILON]
                          input output prefix

### Positional arguments

| Argument | Description                                |
| :------- | :----------------------------------------- |
| `input`  | Multiple sequence alignment (FASTA format) |
| `output` | Output directory                           |
| `prefix` | Prefix for the output files                |

### Options

| Option                                                  | Description                                                                                                                                      | Default |
| :------------------------------------------------------ | :----------------------------------------------------------------------------------------------------------------------------------------------- | :------ |
| `--write-tabular-outputs`                               | Write tabular outputs that describe cluster assignments and cluster statistics                                                                   | False   |
| `--remove-lowercase-columns`                            | Remove columns with at least one lowercase residue (before distance computation)                                                                 | False   |
| `--max-gap-sequence MAX_GAP_SEQUENCE`                   | Remove sequences with gaps representing more than this fraction of the sequence (before distance computation)                                    | 0.25    |
| `--max-gap-column MAX_GAP_COLUMN`                       | Remove columns with gaps representing more than this fraction of the column (after distance computation)                                         | 0.75    |
| `--use-one-hot-encoding`                                | Use one-hot encoding instead of the BLOSUM62-like embedding                                                                                      | False   |
| `--use-hdbscan`                                         | Use HDBSCAN instead of DBSCAN                                                                                                                    | False   |
| `--min-cluster-size MIN_CLUSTER_SIZE`                   | Minimum cluster size to consider                                                                                                                 | 3       |
| `--min-samples MIN_SAMPLES`                             | The minimum number of neighbors required for a sequence to be considered a core data point                                                       | 3       |
| `--min-eps MIN_EPS`                                     | Minimum epsilon value to be considered in the evaluation range (only for DBSCAN)                                                                 | 0.1     |
| `--max-eps MAX_EPS`                                     | Maximum epsilon value to be considered in the evaluation range (only for DBSCAN)                                                                 | 22.1    |
| `--eps-step EPS_STEP`                                   | Step size for incrementing epsilon values within the evaluation range (only for DBSCAN)                                                          | 0.5     |
| `--set-epsilon-value SET_EPSILON_VALUE`                 | Set the maximum distance between two sequences for them to be considered neighbors, disabling the automatic epsilon selection (only for DBSCAN). | None    |
| `--cluster-selection-epsilon CLUSTER_SELECTION_EPSILON` | Clusters with a distance below this threshold will be merged (only for HDBSCAN)                                                                  | 0.0     |

## Example

```sh
# Download the full alignment of the PF01719 family from InterPro and convert it to the FASTA format
$ curl -LJs https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF01719\?annotation\=alignment:full | gzip -dc | esl-reformat --gapsym="-" afa - > PF01719.afa
# Cluster the sequences in the MSA and write the output files to the PF01719_clusters directory
$ ./cluster_msa.py --remove-lowercase-columns --write-tabular-outputs PF01719.afa PF01719_clusters PF01719
```

The script will create a directory named `PF01719_clusters` containing the following files:
- `PF01719_cluster_assignments.tsv`: Tabular file with the cluster assignments for each sequence.
- `PF01719_cluster_data.tsv`: Tabular file with statistics for each cluster.
- `PF01719_cluster_*.afa`: FASTA format MSAs, with each file representing a different cluster.

> [!NOTE]
> The `.tsv` outputs are not written by default. The `--write-tabular-outputs` option is required to generate them.

```sh
$ head -n 5 PF01719_clusters/PF01719_cluster_assignments.tsv | csvtk pretty -t
```

    sequence_accession        cluster_index
    -----------------------   -------------
    A0A919WZ91_9BACI/1-122    -1
    A0A844QCC7_9BACL/2-146    0
    A0A091BW94_9ENTE/5-123    1
    A0A6N6NL64_9ACTN/16-135   -1

> [!NOTE]
> DBSCAN and HDBSCAN may leave some sequences unclustered, so not all sequences in the input MSA will necessarily appear in one of the output FASTA files. Unclustered sequences are indicated in the `<prefix>_cluster_assignments.tsv` output file by a value of `-1` in the `cluster_index` column.

```sh
$ head -n 5 PF01719_clusters/PF01719_cluster_data.tsv | csvtk pretty -t --max-width 35 --clip
```

    cluster_index   cluster_size   avg_identity_to_consensus   consensus_sequence
    -------------   ------------   -------------------------   -----------------------------------
    0               13             0.8099                      -AKDKARYFTFLLYPESIPEDWESKLELLGVP...
    1               7              0.7048                      -KDTRGRNWTFIVYPESAPENWREILDDLHIP...
    2               25             0.8258                      -KEQRSNKWAFLFYQESAPENYLDILEELHIP...
    3               4              0.7769                      ---KRTRNFATVVYPESAPEDWEDKLQEQCIP...

```sh
$ seqkit head -n 3 PF01719_clusters/PF01719_cluster_0.afa
```

    >A0A844QCC7_9BACL/2-146
    NKEKSRYFTFLLYPESIPNDWEMKLELLGLPIAISPLHDRDETEEEGVYKKAHYHVIYIA
    KNPVTVDSVRKRIQRTLGKQSVAMVQSVENVYLYLTHESKKKKHVYDKKDIKLLNNFDID
    >A0A7X6S4J6_9LACO/1-130
    AKDKSRYFTFLLYPESIPEDWKSKLELIGVPIAVSPLHDKDKSAVPGEFKKPHYHVVYVA
    KNPVTADSVRYKIKQLLGDQSIAKVQSMTSMFLYLTHESKKKKHKYNKQDITLINNFDID
    >A0A8J2YP40_9BACL/1-130
    AKDKARYFTFLLYPESIPEDWEMRLESLGVPIAISPLHDKDLSNVEGKYKKAHYHVIYVA
    KNPVTAESVRLKIKRCLGSRSVAMVQSMENIYLYLTHESKKNKHVYDKKDIKLLNNFDID

[^1]: Wayment-Steele, Hannah K., et al. "[Predicting multiple conformations via sequence clustering and AlphaFold2](https://doi.org/10.1038/s41586-023-06832-9)". *Nature* **625.7996** (2024): 832-839.
[^2]: Zhang, Hongyi, Xiaowei Zhan, and Bo Li. "[GIANA allows computationally-efficient TCR clustering and multi-disease repertoire classification by isometric transformation](https://doi.org/10.1038/s41467-021-25006-7)". *Nature Communications* **12.1** (2021): 4699.