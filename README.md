# Spectral Proper Orthogonal Decomposition


## Table of contents

  * [Description](#description)
  * [Installation](#installation)
  * [License](#license)
  * [Contact us](#contact-us)
  * [Contributors](#contributors)


## Description

__SPyOD__ is the python implementation of the Spectral Proper Orthogonal Decomposition published by [Sieber et al. in 2016](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/spectral-proper-orthogonal-decomposition/DCD8A6EDEFD56F5A9715DBAD38BD461A).

It includes two __.py-files__:

- `spod.py` - Includes the function `spod` which calculates the SPOD and

- `findpairs.py` - Includes the post-processing of the SPOD in the function `findpairs` which finds linked modes as described in [Sieber et al. in 2016](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/spectral-proper-orthogonal-decomposition/DCD8A6EDEFD56F5A9715DBAD38BD461A)


and one __jupyter notebook__ example `example_SPOD_free_jet.ipynb` of the SPOD of experimental PIV data from a free jet flow. The data are stored in `PIV_jext_example.mat`.

The paper describing the SPOD method is made publicly available by the TU Berlin at the following link: [https://doi.org/10.14279/depositonce-6377](https://doi.org/10.14279/depositonce-6377).

## Installation 

The __SPyOD__ package can be installed using the following command:
```bash
$ pip install SPyOD
```
The package can be imported by
```bash
$ from spyod.spod import spod
$ from spyod.findpairs import findpairs
```

## License

__SPyOD__ is released under the MIT License. Please have a look at [LICENSE.md](LICENSE.md) for more details.

## Contact us
The best way to contribute is to report any issues you mention for improvement and extension of the code or documentation as well as fixing bugs etc. Feel free to contact us. 


## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center"><a href="https://github.com/grigorishat"><img src="https://avatars.githubusercontent.com/u/114856563?s=400&u=9eea6aaba80fe841c18c8a621111e2d9f3da63ed&v=4" width="100px;" alt="Grigorios Hatzissawidis"/><br /><sub><b>Grigorios Hatzissawidis</b></sub></td>
      <td align="center"><a href="https://github.com/morsieber"><img src="https://avatars.githubusercontent.com/u/116639701?v=4" width="100px;" alt="Moritz Sieber"/><br /><sub><b>Moritz Sieber</b></sub></td>

  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

