# datafunk
Miscellaneous data manipulation tools

## Install
Either pip install using command
```
pip install .
```
or use
```
python setup.py install
```
and test with
```
python setup.py test
```

## Adding functions to this suite
Ideally the function name should be used as the filename and the argparse group. In the following example,
the function name is `remove_dat_junk`.
1. Add your script to directory `datafunk` e.g. `datafunk/remove_dat_junk.py`
2. Update `datafunk/__init__.py` and `datafunk/subcommands/__init__.py` by adding the command name to the `all` lists
3. Add a new command line parameter section to `datafunk/__main__.py`
This should start by defining a new argparse group, e.g.
```
subparser_remove_dat_junk = subparsers.add_parser(
        "remove_dat_junk",
        usage="datafunk remove_dat_junk -i <input>",
        help="Example command",
    )
```
then include all the arguments, e.g.
```
    subparser_remove_dat_junk.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        action="store",
        type=str,
        help="Input file: something about the input file format",
    )
```
and end with the entry point
```
    subparser_remove_dat_junk.set_defaults(func=datafunk.subcommands.remove_dat_junk.run)
```
4. Create file `datafunk/subcommands/remove_dat_junk.py` which defines how to `run` given the command line parameters. Alternatively, specify the entrypoint within the main script file.
5. If you have tests, add the test data to a subdirectory e.g. `tests/data/remove_dat_junk`, and add the test file
`tests/remove_dat_junk_test.py`. This file should contain unit tests which have names `test_*` and ideally be informative about which function they test/the result.
6. If the script has any new dependencies, update `install_requires` section of `setup.py` - this means that it can be pip installed from a conda environment file without a hitch.



Newly added functions in version 0.00000001 (Prototype)

       - clean_names (clean_dat_name)
       - merge_fasta (merge_dat_fasta)
       - remove_fasta (remove_dat_fasta)

Examples to use for comparison in test_files (will be deleted soon)
readme_test for commands

Slack me for bugs (there will be many)
Will update for them in v0.000000002 (Prototype2)
