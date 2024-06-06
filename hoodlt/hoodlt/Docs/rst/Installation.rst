.. _getting-started:

Installation
===========================

How to Get the Program
----------------------

1. **Clone the Repository**:
    ::

        git clone https://your_user_name@bitbucket.org/trvsst/hoodlt.git

2. **Switch to the Appropriate Branch**:
    ::

        git checkout name_of_branch

3. **Ensure You Are in the Appropriate Directory**:
    ::

        ls -al pyproject.toml


Installation
------------

1. **Install Python with Miniconda**:

    Download Miniconda from the following link: https://www.anaconda.com/download/

2. **Create an Environment**:

    Replace `X` with the desired Python version.
    ::

        conda create -n your_environment_name python=3.X

3. **Install pip** (or ensure it is installed):
    ::

        conda install pip

4. **Install Hoomd**:

    Usually installed through conda-forge.
    ::

        conda install hoomd

5. **Install Required Packages**:

    Follow the given order to install the dependencies.
    ::

        conda install freud rowan gsd openpyxl pandas scipy matplotlib importlib_resources


Development Installation
------------------------

For development, use the following command:
    ::

        pip install -e .

For standard installation, use:
    ::

        pip install .


Generating Documentation
------------------------

Make sure LaTeX is installed. Then, navigate to the `docs` directory and generate the documentation using Sphinx:
    ::

        sphinx-build -b html ./rst /directory_where_you_want_your_documentation


Citing Hoodlt
-------------

When citing Hoodlt in your work, please refer to the official documentation on how to properly cite it.

.. note::

    This wiki uses the Markdown syntax. For more information, check out the MarkDownDemo tutorial on Bitbucket, which shows how various elements are rendered. Additionally, the Bitbucket documentation provides more details on using a wiki.
