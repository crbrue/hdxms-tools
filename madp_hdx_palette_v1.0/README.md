# MADP HDX Palette v1.0 bundle

From the root of `hdx-ms-tools`, run:

```bash
bash /path/to/install_madp_hdx_palette.sh "$(pwd)"
```

This installs the palette resources under:

```text
resources/styles/
```

The bundle does not automatically rewrite your current coloring code. Use
`resources/styles/hdx_palette.py` as the canonical Python helper when you
consolidate the PyMOL and ChimeraX exporters.
