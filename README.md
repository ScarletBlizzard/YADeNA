## YADeNA
Yet Another De Novo Assembler

### Contributing
1. Fork and clone the repository
2. Do your changes
3. Install [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)
4. Create conda environment and install dependencies from the ```environment.yml``` file:
    ```
    conda env create -n YADeNA -f environment.yml
    ```
5. Put FASTA file with reference genome into ```test_data/``` directory
6. Install [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
7. Run ```simulate_data.sh```
8. Activate the environment:
    ```
    conda activate YADeNA
    ```
9. Test the program:
    ```
    python test_assembler.py
    ```
    Or run it directly with necessary arguments:
    ```
    python assembler.py ...
    ```
10. Commit, push and do a pull request