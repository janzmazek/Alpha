from setuptools import setup, find_packages

setup(
      name = "CellModel",
      version = "1.0.0",
      description = "Metabolic, Signaling and Electrophysiological Model of Gluucagon and Insulin Secretion",
      author = "Jan Zmazek",
      url = "https://github.com/janzmazek/Alpha-Beta-Model",
      license = "MIT License",
      packages = find_packages(exclude=['*test']),
      python_requires='>3.6.0',
      install_requires = ['numpy', 'scipy', 'matplotlib', 'pyyaml']
)
