from setuptools import setup, find_packages

requirements = [
    # numpy.
    "numpy>=1.15.0",
    # some plots use matplotlib
    "matplotlib >=3.0.0",
    # scipy is used in different places
    "scipy>=1.0.0",
    # freud
    "freud-analysis>=2.0.0",
    # pandas
    "pandas>= 1.0.0",
    # rowan
    "rowan>=1.2.1",
]

description = "Suite of Packages For All Things Nanoparticles"

setup(
      name='hoodlt',
      version='0.8',
      packages=find_packages(),
      description=description,
      url='https://bitbucket.org/trvsst/hoodlt',
      author='Alex Travesset group',
      author_email='trvsst@ameslab.gov',
      license='AmesLab',
      include_package_data=True,
      package_data={'hoodlt/Data': ['Data/nanoparticle_list.xlsx'],
                    'hoodlt/Data/ForceField': ['Data/Forcefield/params_forcefield.xlsx']},
      install_requires=requirements, python_requires=">=3.8",
      zip_safe=False
      )
