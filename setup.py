from setuptools import setup, find_packages

setup(
    name="datafunk",
    version="0.0.5",
    packages=find_packages(),
    url="https://github.com/cov-ert/datafunk",
    license="MIT",
    entry_points={"console_scripts": ["datafunk = datafunk.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=[
        "biopython>=1.70",
        "pandas>=0.25.0",
        "pycountry>=pycountry-19.8.18",
        "pysam>=0.15.4",
        "datapackage"
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
