from setuptools import setup, find_packages

setup(
    name="datafunk",
    version="0.1.1",
    packages=find_packages(),
    package_data={'datafunk':['datafunk/resources/*']},
    include_package_data=True,
    url="https://github.com/cov-ert/datafunk",
    license="MIT",
    entry_points={"console_scripts": ["datafunk = datafunk.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=[
        "biopython>=1.70",
        "pandas>=0.25.0",
        "pycountry",
        "pysam>=0.17.0",
        "datapackage",
        "epiweeks",
        "unidecode",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
