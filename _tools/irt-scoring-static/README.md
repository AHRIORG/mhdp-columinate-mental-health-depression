# IRT Scoring Static Prototype

This folder contains a browser-only proof-of-feasibility for GitHub Pages
hosting.

The prototype:

- loads a browser-safe JSON export of the production 1-factor engines
- scores manual response patterns fully in the browser
- supports browser-local batch CSV scoring
- keeps uploads local to the user session

It is not yet the final production tool. The browser scorer still needs formal
row-by-row validation against the R implementation before it should replace the
hosted Shiny runtime.
