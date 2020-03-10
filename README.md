README
================
Megan Barkdull

## Changes from Ben Rubin’s Comparative Genomics Pipeline

This set of scripts is a slightly modified version of [Ben Rubin’s
Comparative Genomics
Pipeline](https://github.com/berrubin/ComparativeGenomics). I have:

  - Reformatted all python scripts to be compatible with python3.
  - Changed hard-coded paths so that the scripts should be portable
    between different environments.
  - Re-written the ReadMe to reflect how I am using this workflow.

## Using This Workflow:

Note that this workflow must be run from the beginning; even if you
**only** want to convert Orthofinder outputs to RERconverge inputs, all
steps depend on the output of steps before them.

### Getting Started:

Download this repository into your working directory.

Make sure you have the required inputs. To begin using this workflow,
you will need:

  - A file that identifies orthogroups.
      - The default is to read the output format of
        [Orthofinder](https://github.com/davidemms/OrthoFinder).
      - When identifying orthogroups, it is necessary that all gene
        names be prefaced by a four-character taxon designation code.
      - See the Example Orthogroups file for an example.
  - The coding sequences for each gene listed in the orthogroups file.
      - The coding sequence files are specified in a two-column
        parameters file; the first column gives the four-character taxon
        designation code, and the second column gives the path to the
        coding sequence fasta file for that taxon.
      - See the Example Coding Sequence Parameters file for an example.

Note that the gene names in the orthogroup file **must exactly match**
the gene names in the coding sequence fasta files.

#### Potential problems:

This step of the workflow relies on several tools which you may need to
install, and which then may need to be [added to your path
variable](https://superuser.com/questions/284342/what-are-path-and-other-environment-variables-and-how-can-i-set-or-use-them),
including:

  - [BioPython](https://biopython.org/wiki/Download)
  - [pyvcf](https://anaconda.org/bioconda/pyvcf)
  - [PySAL](https://pysal.org/install)
  - [PAML](https://bioconda.github.io/recipes/paml/README.html)
  - [statsmodels](https://www.statsmodels.org/stable/install.html)
  - [FSA](http://fsa.sourceforge.net/)
  - [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html)
  - [trimAl](http://trimal.cgenomics.org/)

### 1\. Gathering and Initially Filtering Orthogroups

#### What does this step do?

This first step will create .fasta files for each orthogroup that meet
the filtering requirements. Several filtering steps are hard-coded here:

  - Any sequence that is less than half of the median length of all the
    sequences in an orthogroup is removed from that orthogroup.
  - Any species with more than three sequences in an orthogroup is
    removed from that orthogroup.
  - Orthogroups where the total number of sequences is more than 1.5
    times the number of species in the orthogroup are completely
    removed.

You also specify the minimum number of taxa that must be represented in
an orthogroup in order for that orthogroup to be retained.

#### How do I do it?

To gather the orthogroups, use the command (may need to be prefaced by
`python`):

`selection_pipeline.py -a write_orthos -b [output directory] -o [output
prefix] -r [orthogroup file] -t [min. taxa to include locus] -d [params
file]`

   `-a` is the action to perform; in this case, we are doing
`write_orthos`  
   `-b` is the output directory, which will be created in this step; all
of your outputs, for every step of this workflow, should be written to
this same output directory.  
   `-o` is the prefix that will be used for labelling all output in your
project. Make this something logical and informative\!  
   `-r` is the path to the orthogroup file, as described under
**Required Inputs**.  
   `-t` is the minimum number of taxa that must be represented in an
orthogroup in order for that orthogroup to be retained through
filtering. Choose a logical proportion of your total number of taxa
here.  
   `-d` is the path to the parameters file that gives the location of
the coding sequence file for each of your taxa, as described under
**Required Inputs**.

#### What are the outputs?

This command creates:

  - A directory `./[output directory]/[output_prefix]_orthos`, which
    contains fasta files of the coding sequence for every orthogroup.
  - A directory `./[output directory]/[output_prefix]_orthos_prots`,
    which contains fasta files of the translated amino acid sequence for
    every orthogroup.
  - An index file, `[output_prefix]_ortho.index`, which lists the number
    of taxa and number of sequences present for each orthogroup.
      - This is a critical file that acts as a reference for future
        steps which process and choose orthogroups.

### 2\. Aligning Orthogroups

#### What does this step do?

Now we must align the written orthogroups, so that homologous
nucleotides are at the same position in every sequence. This step uses
[FSA](https://anaconda.org/bioconda/fsa) with the `--nucprot` option to
align the coding sequences.

#### How do I do it?

To align the written orthogroups, use the command (may need to preface
with `python`):

`selection_pipeline.py -a align_coding -p [number of threads to use] -b
[output directory] -o [output prefix] -r [orthogroup file] -t [min. taxa
to include locus] -d [params file]`

   `-a` is the action to perform; in this case, we are doing
`align_coding`  
   `-p` is the number of threads you want to use to run the alignment.
This step can take a while, so using multiple threads is a good idea.  
   `-b` is the output directory, which was created in the previous step;
all of your outputs, for every step of this workflow, should be written
to this same output directory.  
   `-o` is the prefix that you are using to label all output in your
project.  
   `-r` is the path to the orthogroup file, as described under
**Required Inputs**.  
   `-t` is the minimum number of taxa that must be represented in an
orthogroup in order for that orthogroup to be retained through
filtering. Choose a logical proportion of your total number of taxa
here.  
   `-d` is the path to the parameters file that gives the location of
the coding sequence file for each of your taxa, as described under
**Required Inputs**.

#### What are the outputs?

This command creates:

  - A directory, `./[output directory]/[output_prefix]_fsa_coding`, that
    contains unaligned (\*.fa) and aligned (\*.afa) fasta files for each
    orthogroup.
  - a file, `[output_prefix].afa`, containing an aligned, concatenated
    matrix of all proteins from all orthogroups with a single sequence
    per species.
      - The protein sequences are from the
        [trimAl](http://trimal.cgenomics.org/) alignments, and so are
        fairly conservatively filtered.
      - This matrix can be used to infer a phylogeny, and is formatted
        to work with
        [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html).

### 3\. Filtering the Alignments

#### What does this step do?

Next, we will filter the alignments, to remove low confidence positions
and low information or low confidence sequences. By default, Gblocks and
trimAl are run on each orthogroup; however, both of these are fairly
stringent methods that may remove potentially valuable informative
sequence. Therefore, this step implements several other post-alignment
filters, which can be adjusted to suit your data:

1.  Columns in the alignment with fewer than a specified number of
    non-gap characters are removed.
2.  Columns in the alignment with less than a specified proportion of of
    known nucleotides (anything that is not “-” or “N”) are removed.
3.  Columns in the alignment that contain sequences from fewer than a
    specified number of taxa can be removed.
4.  Then the amino acid alignments are run through the [Jarvis et al.
    (Science 346: 1320-1331) Avian Phylogenomics
    Project](https://science.sciencemag.org/content/346/6215/1320.full.pdf+html)
    scripts for filtering amino acid alignments, which mask over poorly
    aligning regions of individual sequences (rather than omitting
    entire alignment columns).
      - The scripts spotProblematicSeqsModules.py and
        spotProblematicSeqsModules-W12S4.py were downloaded from
        <ftp://parrot.genomics.cn/gigadb/pub/10.5524/101001_102000/101041/Scripts.tar.gz>
        on 1/31/2019.
      - These scripts were incorporated into the selection\_pipeline.py
        pipeline through the jarvis\_filtering() function.
5.  The outputted, masked sequences are passed through the first three
    filters again.
6.  Sequences that now contain fewer than a specified proportion of
    their original, known nucleotides are removed entirely.
7.  Sequences that contain more than a specified proportion of unknown
    sequence are removed entirely.
8.  Sequences with fewer than a specified total number of known
    nucleotides are removed entirely.

#### How do I do it?

To filter the aligned orthogroups, use the command (may need to preface
with `python`):

`selection_pipeline.py -a alignment_filter -b [output directory] -o
[output prefix] -r [orthogroup file] -p [number of threads to use] -t
[min. taxa to include locus] -d [parameters file] --nogap_min_count
[filtering step 1] --nogap_min_prop [filtering step 2]
--nogap_min_species [filtering step 3] --min_seq_prop_kept [filtering
step 6] --max_seq_prop_gap [filtering step 7] --min_cds_len [filtering
step 8]`

   `-a` is the action to perform; in this case, we are doing
`alignment_filter`  
   `-p` is the number of threads you want to use to run the filtering.
The Jarvis filter can take a while, so you should use multiple threads
here.  
   `-b` is the output directory, which was created in the first step;
all of your outputs, for every step of this workflow, should be written
to this same output directory.  
   `-o` is the prefix that you are using to label all output in your
project.  
   `-r` is the path to the orthogroup file, as described under
**Required Inputs**.  
   `-t` is the minimum number of taxa that must be represented in an
orthogroup in order for that orthogroup to be retained through
filtering. Choose a logical proportion of your total number of taxa
here.  
   `-d` is the path to the parameters file that gives the location of
the coding sequence file for each of your taxa, as described under
**Required Inputs**.  
   `--nogap_min_count` is the minimum number of non-gap characters that
must be present in a column, and is used in filtering step 1.  
   `--nogap_min_prop` is the required proportion of known nucleotides in
a column to retain that column, and is used in filtering step 2.  
   `nongap_min_species` is the required number of species with a non-gap
character in a column, and is used in filtering step 3.  
   `--min_seq_prop_kept` is the proportion of original known nucleotides
that must be retained in a sequence, and is used in filtering step 6.  
   `--max_seq_prop_gap` is the maximum allowed proportion of unknown
sequence in each sequence, and is used in filtering step 7.  
   `--min_cds_len` is the minimum known sequence length required after
filtering, and is used in filtering step 8.

#### What are the outputs?

This command creates:

  - A directory, `./[output
    directory]/[output_prefix]_fsa_coding_columnfilt`, that contains the
    aligned fasta files produced by the first three filtering steps.
  - A directory, `./[output directory]/[output_prefix]_coding_jarvis`,
    that contains the masked output from filtering step 4.
  - A directory, `./[output
    directory]/[output_prefix]_coding_jarvis_columnfilt`, containing the
    masked, column-filtered output from step 5.
  - A directory, `./[output
    directory]/[output_prefix]_coding_jarvis_columnfilt_seqfilt`,
    containing the outputs from filtering steps 6-8.

### 4\. Compiling Data for RERconverge

#### What does this step do?

[RERconverge](https://github.com/nclark-lab/RERconverge) is a tool that
estimates correlations between rates of molecular evolution and the
evolution of a trait. RERconverge requires as input a phylogeny for each
gene of interest, with branch lengths inferred from degree of protein
sequence divergence.

This step of the workflow creates those phylogenies, using
[AAML](http://web.mit.edu/6.891/www/lab/paml.html) to infer branch
lengths for the orthogroups identified by Orthofinder.

#### How do I do it?

To filter the aligned orthogroups, use the command (may need to preface
with `python`):

`selection_pipeline.py -a rer_converge -p [number of threads] -b [output
directory] -o [output prefix] -t [min. taxa to include locus]
--outputfile [output file name] --taxa_inclusion [taxon requirement
file] -e [newick species tree file]`

   `-a` is the action to perform; in this case, we are doing
`rer_converge`  
   `-p` is the number of threads you want to use.  
   `-b` is the output directory, which was created in the first step;
all of your outputs, for every step of this workflow, should be written
to this same output directory.  
   `-o` is the prefix that you are using to label all output in your
project.  
   `-t` is the minimum number of taxa that must be represented in an
orthogroup in order for that orthogroup to be used in this step. Choose
a logical proportion of your total number of taxa here.  
   `--outputfile` is the name you want your output file to have.  
   `--taxa_inclusion` is the path to a tab-delimited file that specifies
the patterns of included taxa you require be present in the orthogroups
passed to RERconverge. See [Ben Rubin’s
ReadMe](https://github.com/berrubin/ComparativeGenomics) for details.  
   `-e` is the path to a Newick format species tree; this is one of the
outputs of Orthofinder, so you can go ahead and use that.

#### What are the outputs?
