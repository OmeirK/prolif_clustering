# prolif_clustering
Series of scripts to cluster protein–ligand crystal structures on the basis of ProLIF similarity.

Code assumes that all input proteins have the same residue numbering. Fragalysis directory format is used as the input. Examples of this format are in the `examples/aligned_files` directory.

## 1. Calculate ProLIF Fingerprints

Save fingerprints in a tsv file.

``` 
python3 01_Py_get_prolif_fps.py -fd=example_inputs/aligned_files/ -o=tsv_prolif_fps.tsv
```
Any failed calculations will be saved to `tsv_prolif_fps.err`

Next, convert the fingerprint information to a sparse bit fingerprint for further calculations:

```
python3 02_Py_convert_prolif_to_fp.py -i=tsv_prolif_fps.tsv -o=pkl_sparse_fps.pkl
```

## 2. Cluster ligands by fingerprint

Save a json file with clustered input targets. Users can set a distance threshold for the Butina clustering algorithm. Please note that the threshold represents distane, where `distance = (1 - tanimoto_similarity)`.

```
python3 03_Py_calc_tanimoto_similarity.py -i=pkl_sparse_fps.pkl -c=0.5 -o=pairwise_clusters.json
```
