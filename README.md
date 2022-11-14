# Spectral Proper Orthogonal Decomposition


## Table of contents

  * [Description](#description)
  * [Installation](#installation)
  * [License](#license)
  * [Contributors](#contributors)


## Description

__SPyOD__ is the python implementation of the Spectral Proper Orthogonal Decomposition published by [Sieber et al. in 2016](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/spectral-proper-orthogonal-decomposition/DCD8A6EDEFD56F5A9715DBAD38BD461A).

It includes two __.py-files__:

- calc_SPOD.py - Includes the function `spod` which calculates the SPOD and

- find_mode_pairs_from_SPOD.py - Includes the post-processing of the SPOD in the function `find_spod_pairs` which finds linked modes as described in [Sieber et al. in 2016](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/spectral-proper-orthogonal-decomposition/DCD8A6EDEFD56F5A9715DBAD38BD461A)


and one __jupyter notebook__ example `example_SPOD_free_jet.ipynb` of the SPOD of experimental PIV data from a free jet flow. The data are stored in `PIV_jext_example.mat`.

## Installation 

The SPyOD package can be installed using the following command:
```bash
$ pip install spyod
```
## Licence

__SPyOD__ is released under the MIT License. Please have a look at [LICENSE.md](LICENSE.md) for more details.

## Contributors

Thanks go to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

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

