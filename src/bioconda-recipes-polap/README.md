# Upload the packages to my channel of anaconda.org

bioconda-utils build --packages cflye
anaconda upload cflye-2.9.5.1-py3*
anaconda upload dflye-2.9.5.0-py3*

bioconda-utils build --packages polap
bioconda-utils build --packages cflye
conda install anaconda-client
bioconda-utils build --packages cflye --python 3.12
cd ~/miniconda3/envs/bioconda/conda-bld/linux-64/
anaconda login
anaconda upload cflye-2.9.5.0-py3\*
