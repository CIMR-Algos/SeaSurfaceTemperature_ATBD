# CIMR Level-2 Sea Surface Temperature ATBD

## Fetch and read the book
If you just want to download the repo and read the book (not building it):
```
cd </path/to/where/you/want/the/repo> # call it /PATH/
git clone <repo>
```
In your web browser load `file:///PATH/book/_build/html/intro.html` in the url field.

## Instructions for building the book

### If needed, install miniconda
See https://docs.conda.io/en/latest/miniconda.html

### Create a conda environment and install the jupyter-book tools
```
conda create -c conda-forge --name JB python=3.8
conda activate JB
conda env list # check that JB is the current env
pip install --upgrade pip #(if needed several times)
pip install jupyter-book
```

### install packages needed for the book (this list might grow as we develop the books)
```
conda install -c conda-forge matplotlib numpy
conda install -c conda-forge sphinxcontrib-mermaid
```

### alternatively, install the necessary dependencies from devalgo-ENV.yml
```
conda env create -n devalgo --file devalgo-ENV.yml
```

### build the book
```
cd /path/to/local/repo/with/atbd/ # the root of the repo, with book/ algorithm/ data/ ...
jupyter-book build book/
jupyter-book build --all book/ #(force rebuild everything)
```

### develop the book and commit to GitHub
Edit the `.md` files of the book
Build the book (see above)
```
git add / commit #possibly in two commits, first the changes to the `.md` files, then the modified `_html` files
git push
```

## Other Resources:
* [Setting up repos on GitHub](https://kbroman.org/github_tutorial/pages/init.html)
* [Put your JupyterBook on Github with online webpage (requires Public repo)](https://github.com/pabloinsente/jupyter-book-tutorial)
