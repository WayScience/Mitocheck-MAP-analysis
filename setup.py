from setuptools import setup, find_packages

setup(
    name="mitocheck-map-analysi",
    version="0.0.1",
    packages=find_packages(),  # Automatically discover and include all packages
    author="Erik Serrano",
    author_email="erik.serrano@cuanschutz.edu",
    description="Utilizing the mean average precision (MAP) metric to assess"
                "reproducibility and perturbation effect on single-cell profiles"
                "in the MitoCheck dataset.",
    url="https://github.com/WayScience/Mitocheck-MAP-analysis",
)

