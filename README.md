# OriC_Finder
Python scripts that predict and plot the location of the origin of replication (oriC) of circular bacterial genomes based on Z-curve and GC-skew analysis.
**NOTE**: This is a work in progress.

### Order of operation
The scripts are excecuted in a specific order to work properly, but also work independently, so that each script can serve as a checkpoint.

#### DoriC data prep
The DoriC data can be downloaded from http://tubic.tju.edu.cn/doric/public/index.php as a .RAR. Unpack this however you want and you'll be left with a CSV-file.

1. `data_prep_doric.py`: Creates three new CSV-files (only \_concat.csv works for now) that have each ordered the relevent DoriC data slightly differenly.

#### NCBI data prep
These scripts prepare the NCBI data for analysis. Each script has docs-strings for further information.

1. `ncbi_download.py`: Use this script to download a dataset of your choice. Documentation for the `ncbi-genome-download` package can be found here: https://github.com/kblin/ncbi-genome-download.
2. `ncbi_to_fasta.py`: Unzips and extracts the downloaded FASTA-files from the dataset. Multiple cleaning/filtering options available.
3. `fasta_to_oriC_csv.py`: Predicts the oriC(s) for the whole dataset.

#### Comparison
Once both the DoriC and NCBI datasets have been processed, they can be compared. This is done with `oirC_comparison.py`.

### `oriC_Finder.py`
This script predicts the origin of replication for circular bacterial DNA. It makes use of a combination of Z-curve and GC-skew analysis.
- https://en.wikipedia.org/wiki/Z_curve

### `Plotting_functions.py`
There are 3 general functions in this file which can be used to plot any generic 1D-np.array.
- `plot_Z_curve_3D`: Makes a 3D-plot of the Z-curve.
- `plot_Z_curve_2D`: Can plot a maximum of four axes in a single 2D-plot. This one is useful for plotting single or multiple Z-curve or GC-skew component agaist each other.
- `plot_GC_skew`: Does the same as `plot_Z_curve_2D`, except only takes one array.
