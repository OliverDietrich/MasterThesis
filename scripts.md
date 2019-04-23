# Automating the computational workflow using R scripts

## Patching
Patches are used to update the commonly used version of a script (e.g. dropseQC.R) with changes from the latest draft (e.g. dropseQC-1.7) without replacing the file.

The code used for patching such a script is:

> diff -c dropseQC.R dropseQC-1.x.R | patch

Alternatively, the patch (file containing only the changes) can be distributed and applied to the older script version.

> diff -c dropseQC.R dropseQC-1.x.R > patch-1.x

> patch -i patch-1.x

