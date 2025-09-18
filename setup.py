from setuptools import setup, find_packages
import os

setup(
    name="HLA-PEPCLUST",
    version="1.1.1-dev",
    author="Sanjay Krishna",
    author_email="sanjay.sondekoppagopalakrishna@mail.com",
    packages=find_packages(),
    long_description=open("README.md").read(),
    install_requires=[
        "rich==13.9.4",
        "click==8.1.8",
        "pandas==2.2.3",
        "seaborn==0.13.2",
        "matplotlib==3.9.4",
        "rich-argparse==1.7.0",
        "opencv-python",
        "altair==5.5.0",
        "vl-convert-python==1.7.0",
        "ipykernel==6.30.0",
        "scipy==1.16.1",
        "beautifulsoup4==4.13.4",
        "requests==2.32.5",
        "numba==0.61.2"
    ],
    entry_points={
        "console_scripts": [
            "clust-search=cli.main:main",
        ],
    },
)
data_path = os.getenv('GLODBA_DATA_PATH', 'data/')