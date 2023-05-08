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
6. Put FASTA file with reference genome into ```test_data/``` directory
7. Install [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
8. Run ```simulate_data.sh```
9. Activate the environment:
    ```
    conda activate YADeNA
    ```
10. Test the program:
    ```
    python test_assembler.py
    ```
    Or run it directly with necessary arguments:
    ```
    python assembler.py ...
    ```
11. Commit, push and do a pull request