from setuptools import setup, find_packages

setup(
    name="h-perforatum-network-tox",
    version="1.0.0",
    description="Network pharmacology analysis of H. perforatum hepatotoxicity",
    author="Antony Bevan",
    python_requires=">=3.9",
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
    ],
    entry_points={
        "console_scripts": [
            "run-analysis=network_tox.pipeline:main",
        ],
    },
)
