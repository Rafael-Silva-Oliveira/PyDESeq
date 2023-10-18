# PyDESeq

A DESeq Pipeline run in python for pre-processing, differential expression, GEAS, plotting and more.

# How is the project structured?

The main script to execute the pipeline can be found under the `src\scripts\main.py`. This file is associated with the `src\scripts\config.json` to change different settings on differential expression, GSEA and plotting. The core of the script is under the `src\scripts\deseq.py` script. This script contains a class object with differential expression and GSEA methods. An auxiliary script for plotting is located in the `src\scripts\plots.py`.

# How to run the pipeline?

First, install the `requirements.txt` file using the following CMD promot: `pip install -r requirements.txt`. You can also run the `setup.py` to install the github repository by running the command `python setup.py develop --user` being under the path that contains the setup.py script.

After the dependencies have been installed, run the pipeline by executing the script `src\scripts\main.py`. This script should be able to run without changing any paths. If path errors pop up, it might be due to paths on the `main.py` script (parameter `CONFIG_PATH` in the `run_pipeline` function) or the `count_matrix` parameter in the `config.json` file. Make sure to add the correct paths according to the sub-directories where the folder was installed to.

## Project Organization

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── scripts       <- Scripts to turn the PyDESeq pipeline
    │   │   ├── main.py   <- main script to run the pipeline
    │   │   ├── deseq.py  <- script with the DESeq pipeline as a class object
    │   │   ├── plots.py  <- script to create the plots
    │   │   ├── utils.py  <- helper methods and functions
    │   │   └── config.json   <- configuration file
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── plots  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io

---

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
