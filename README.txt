Before running lp constructors, ensure that proper alias files for
prize assignment and conversion between Reactome and UniprotKB identifiers
have been set up:

TOPFOLDER
  programs
    output
  datasets
    aliases
      *-conversion.txt
      *-name2prize.txt
    hypergraphs
      *-edges.txt
      *-nodes.txt
    reactome
      converted
        *-hyperedges.txt
        *-hypernodes.txt
      parsed
        *-complexes.txt
        *-controls.txt
        *-elements.txt
        *-entitysets.txt
        *-reactions.txt
        *-subpathways.txt

These directories can be initialized with initialize_directories.sh (eventually
will be master_script)

Sample lp file construction
$ build_lp.py Singaling-by-Hedgehog
