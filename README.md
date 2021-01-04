# Primer_Blast_like

**Aim:**

To carry out *[NCBI Primer-Blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)-like* function locally with unpublished genome database, and maybe for further use.

**Current plans:**

0. Local blast with NCBI Blast+ through Biopython.
1. Regular design of PCR primers and locally blast for specificity.  
↘︎ Regardless of the position of the primers
2. Design primers for qPCR in case of needing Exon junction span or Intron inclusion.
3. Add a visualize function to visualize every primer pair's position on the given gene sequence with each use.
4. While dealing with species like *B. napus*(2n = 4x = 38, AACC), each gene could have multiple copies(or homologous genes). For initial expression analysis, input multiple gene data and get consensus primers.
5. For detecting SNP of the same gene of different cultivars(same species), design Allele-specific PCR primers.

**Required environment:**

1. [Python3]()

2. [NCBI Blast+]()

3. Database created with NCBI Blast+

4. Certain directory in order for outputting

5. Python Packages: 

    * [Biopython](https://biopython.org/)

    * os

    * [Pandas](https://pandas.pydata.org/)

    * [Primer3-py](https://pypi.org/project/primer3-py/)

**Required parameters:**

↘︎ Check [Primer3 Documentation](http://primer3.org/manual.html) for detail

1. Email (in case of contact from NCBI)

2. Gene id, Database (to download sequence from NCBI GenBank)

   ↘︎[NCBI efetch database](https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly)

3. Output address, Database address (to save the result)

4. E-value, Identity (to check specificity)


**Optimal parameters: (for designing primers)**

1. Gene name (to name saved files)

2. Other parameters (for particular need in designing primers)
