# rdsss

CLI tool of RDKit rdSubstructLib

## Description

CLI tool of substructure search

## Getting Started

### Dependencies

* click
* rdkit

### Installing

```
$ gh repo clone iwatobipen/rdsss
$ cd rdsss
$ pip install -e .
```

### Executing program / basic usage

```
# make ssslib from sdf.gz
$ make_rdssslib <input.sdf.gz> <output.sslib.pkl>

# search with ssslib from CLI
$ run_rdsss 'SMARTS query' <output.ssslib.pkl>
```

## Help

## Authors

* iwatobipen

## Vesion History

* 0.1
  * Initial Release

## License

This project is licensed under the [MIT] License - see the LICENSE.md file for details

## Acknowledgements

