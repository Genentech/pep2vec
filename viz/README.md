# Setup Code
install mamba
mamba create -n pep2vec_viz python=3.11
mamba activate pep2vec_viz
mamba install nodejs
pip install -r requirements_frozen.txt

# Download Data
Download this file to the viz directory
https://zenodo.org/records/13932198/files/df_for_vizualization.parquet?download=1

# Running
Use the following command to run the server
```
mamba activate pep2vec_viz
panel serve start_server.py --autoreload --show
```
autoreload will auto restart the server when source code changes, this is optional.
Terminal will contain useful messages, and upon boot, will include the url needed to access.

# Notes:
    - If desired inside the update_df_selected function, you can add a dump of the selected dataframe to a temp file for further processing.
    - Edit start_server.py to change input dataframe as well as pre filtering or whatnot that is desired
