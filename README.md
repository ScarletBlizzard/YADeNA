## YADeNA
Yet Another De Novo Assembler

### Contributing
1. Fork and clone the repository
2. Do your changes
3. Install [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)
4. Add bioconda channel:
    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    ```
5. Create conda environment and install dependencies from the ```environment.yml``` file:
    ```
    conda env create -n YADeNA -f environment.yml
    ```
6. Extract ```out_dir.zip``` into ```test_data/``` OR
    1. Put FASTA file with reference genome into ```test_data/``` directory
    2. Install [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
    3. Run ```simulate_data.sh```
7. Activate the environment:
    ```
    conda activate YADeNA
    ```
8. Make your changes and, if necessary, install or remove packages using conda
9. Test the program:
    ```
    python test_assembler.py
    ```
    Or run it directly with necessary arguments:
    ```
    python assembler.py ...
    ```
10. Before committing, if you installed or removed packages, run the command
    ```
    conda env export > environment.yml
    ```
11. Commit, push and do a pull request