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
4. Create file `datafunk/subcommands/remove_dat_junk.py` which defines how to `run` given the command line parameters
5. If you have tests, add the test data to a subdirectory e.g. `tests/data/remove_dat_junk`, and add the test file
`tests/remove_dat_junk.py`
