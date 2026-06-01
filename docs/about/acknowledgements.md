# Acknowledgements

ICE is built upon several open-source projects.

## FiNeR

**Fortran INI ParseR and generator**

- **Repository:** [github.com/szaghi/FiNeR](https://github.com/szaghi/FiNeR)
- **License:** GPL v3.0

FiNeR is a pure Fortran 2003+ OOP library for reading and writing INI configuration files. ICE uses FiNeR to parse the `input.ini` parameter file.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

ORION provides built-in functions to read and write files in different formats. ICE uses ORION for solution output in Tecplot and VTK formats and for reading the structured-grid mesh.

## Documentation Tools

### MkDocs

**Static Site Generator**

- **Website:** [mkdocs.org](https://www.mkdocs.org)
- **License:** BSD-2-Clause

MkDocs transforms ICE's documentation into a searchable website.

### Material for MkDocs

**Modern Documentation Theme**

- **Website:** [squidfunk.github.io/mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
- **License:** MIT

Material for MkDocs provides the interface for ICE's documentation.

## License Compliance

ICE respects all licenses of its dependencies:

- **GPL v3.0** — ORION, FiNeR

See the [License](license.md) page for ICE's full license text.

---
