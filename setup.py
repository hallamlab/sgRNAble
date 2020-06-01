import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sgRNAble",
    version="0.0.1",
    author="Avery Noonan, Siddarth Raghuvanshi, Ahmed Abdelmoneim",
    author_email="sidr97@gmail.com",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Siddarth-Raghuvanshi/CRISPR-Guide-RNA",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        "console_scripts": ['sgrnable = optimal_guide_finder.guide_finder:main']
    },
    install_requires=[
        'matplotlib==3.2.1',
        'pandas==0.24.2',
        'XlsxWriter==1.1.5',
        'numpy==1.16.2',
        'scipy==1.4.1',
        'biopython==1.76',
        'GPy==1.9.9',
        'ipython==7.10.1',
        'azimuth==2.0',
        'elevation==1.0.6',
        'hyperopt==0.2.2',
        'ipdb==0.12.3',
        'mkl==2019.0',
        'scikit_learn==0.20.2',
        'theanets==0.7.3',
        'numba==0.48.0',
    ],
)
