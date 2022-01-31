"""
Setup file for making the package pip installable
"""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sgRNAble",
    version="1.0.11",
    author="Avery Noonan, Siddarth Raghuvanshi, Ahmed Abdelmoneim",
    author_email="sidr97@gmail.com",
    description="CRISPR-Cas9 single guide RNA generation tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hallamlab/sgRNAble",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    python_requires='==3.7.*',
    entry_points={
        "console_scripts": ['sgrnable = optimal_guide_finder.guide_finder:main']
    },
    install_requires=[
        "biopython",
        "matplotlib",
        "numba",
        "numpy==1.16.6",
        "pandas",
        "psutil",
        "scikit-learn==0.20.4",
        "scipy",
        "tqdm",
    ],
)
