from setuptools import setup, find_packages

setup(
    name="h-perforatum-network-tox",
    version="2.1.0",
    description="Network pharmacology analysis of H. perforatum hepatotoxicity",
    author="Antony Bevan",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires='>=3.10',
    install_requires=[
        'pandas>=2.3.0',
        "numpy>=1.24.0",
        "networkx>=3.0",
        "scipy>=1.10.0",
        "statsmodels>=0.14.0",
        "pyarrow>=12.0.0",
        "matplotlib>=3.7.0",
        "tqdm>=4.65.0",
        "requests>=2.28.0",
        "rdkit>=2023.3.1",
    ],
)
