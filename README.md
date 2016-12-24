README.txt for DarkHorse version 2.0
December 2016
Author: Sheila Podell

----------------------------------------
CONTENTS

 - DESCRIPTION
   
 - SYSTEM REQUIREMENTS/DEPENDENCIES
 
 - INSTALLATION
 
 - TESTING

 - USAGE INSTRUCTIONS

 - VERSION HISTORY

 - CONTACT INFORMATION

----------------------------------------
DESCRIPTION

    Darkhorse is an experimental program that defines phylogenetic
    relatedness of blastp hits for a set of proteins against a taxonomically
    diverse reference database (e.g. NCBI Genbank nr) using a lineage 
    probability index (LPI) score. The basic lgorithm used to calculate 
    LPI scores and its application in predicting
    horizontal gene transfer are described in the following publications. These
    publications should be cited in any work that uses the algorithm or results
    obtained using this software.
    
        Podell, S. and Gaasterland, T. (2007). DarkHorse: A method for 
        genome-wide prediction of horizontal gene transfer. Genome Biology 
        8(2):R16. 
        
        Podell, S Gaasterland, T, and Allen, EE (2008). A database of
        phylogentically atypical genes in archaeal and bacterial
        genomes, identified using the DarkHorse algorithm. BMC
        Bioinformatics 9:419          
    
    Those desiring to incorporate the DarkHorse software, the associated
    DarkHorse HGT Candidate database, or information downloaded from the
    database into commercial products, or to use any of these materials for
    commercial purposes, should contact the Technology Transfer &
    Intellectual Property Services, University of California, San Diego,
    9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, Ph: (858)
    534-5815, FAX: (858) 534-7345, E-MAIL: invent@ucsd.edu.
    
----------------------------------------
SYSTEM REQUIREMENTS/DEPENDENCIES
   
   Hardware requirements include sufficient disk space to store reference protein
   sequences, large MySQL relational database tables based on these sequences, 
   and tab-delimited BLASTP sequence comparisons for query sequences.  
   Total disk space required is roughly triple the size of the decompressed protein 
   reference database. Program performance depends primarily on MySQL database 
   access efficiency. A multi-processor CPU with at least 16 GB of RAM is recommended. 
     
    DarkHorse is a command-line only program, and requires the Unix OS. 
    It has been installed, tested, and run on multiple versions of Linux and 
    Macintosh OS X platforms. Additional software pre-requisites include: 
    
        PERL version 5.8.1 or later, with the following modules installed 
        (http://www.cpan.org/):
            DBI
            DBD::mysql
        Mysql Relational Database software version 5.5 or later (or MariaDB equivalent)
        	(https://dev.mysql.com/downloads/mysql/ or https://downloads.mariadb.org)     	
        Diamond software for performing protein sequence similarity searches:
       		(https://github.com/bbuchfink/diamond)   
            
    Prior to DarkHorse program installation, a database must be created 
    from within the Mysql program to accept DarkHorse input, e.g.
    
        mysql> create database darkhorse_01;
        
    The following public database information must be downloaded,  
    de-compressed, and locally available:
    
        Genbank nr database raw fasta file
            wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
        Genbank Taxonomy database 
	    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
            wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
        	      
    Users must be able to obtain protein BLAST (BLASTP) search data for query
    sequences versus a large, taxonomically informative set of reference sequences, 
    (e.g. Genbank nr) with output in the same  tab-delimited format as 
    NCBI BLAST+ outfmt 6 (equivalent to old NCBI blastall format -m8). The program 
    Diamond is highly recommended for this purpose, although NCBI blast or any 
    alternative (e.g. cluster-accelerated) software may also be used instead, as long 
    as the output ends up in the same tab-delimited format. These tab-delimited output 
    files are used as DarkHorse program input.  
    
    The DarkHorse program is critically dependent on pre-computed associations of
    BLAST search output with taxonomic classifications that are stored in a MySQL relational
    database. The program will not work unless all database sequence ID numbers from
    the blast search match ID numbers in the MySQL database EXACTLY.
	
    To make this easier, the DarkHorse installer program now creates a fasta-format file 
    of taxonomically informative protein sequences suitable as blast search references. 
    This file (or a subset of its sequences) must be used as the reference database 
    for BLAST search input to the DarkHorse program. It is possible to add custom, unpublished 
    sequences to the reference database, but this must be done using the specific scripts 
    and procedures described below.
	     
----------------------------------------
INSTALLATION
       
    Installation is performed using a script that tests availability
    of pre-requisite programs and databases and creates a customized
    configuration file, based on interactive user input. Optionally, a
    template configuration file is provided (templates/config_template),
    which can be manually edited and specified as an additional command 
    line parameter to the installation program:
       
        ./install.pl  -p program_directory_path 
            
            OR
            
         ./install.pl  -p program_directory_path -c config_filename   
    
    Note that the installation script may take a substantial amount of time to run 
    (e.g many hours). If more than 16GB of RAM are available at the time of program execution, 
    performance can be accelerated by adjusting the "max_lines_per_packet" parameter 
    in the configuration file. However, if this parameter is set too high, MySQL operations
    will fail due to insufficient memory. Recommended settings are as follows, based on total
    system memory:
    
    	16 GB RAM  [max_lines_per_packet]=4000  (default value)
    	32 GB RAM  [max_lines_per_packet]=8000
    	64 GB RAM  [max_lines_per_packet]=16000
	
   In addition to creating and populating required MySQL tables, the installation script produces
   a new configuration file (if none was specified on the command line), a log file containing 
   detailed notes on any errors that may have occurred, and new fasta formatted sequence file
   of taxonomically informative sequences. The taxonomically informative fasta file will be 
   named using the following convention:
            
                databasename_informative_ref_seqs.fasta    
   
   Users will need to format this file of protein sequences as a reference database for subsequent 
   blastp searches using Diamond or NCBI-blast software.   	

----------------------------------------

INSTALLATION TESTING

    NOTE: Revised installation testing software for DarkHorse version 2 is still under development,
    and is not yet included in the product package.
    
    After DarkHorse installation is complete, it can be tested using the
    sample data and test script provided:
    
         ./test_dh.pl  -c configuration file
             
    Results of this test are written to a new directory called
    "test_results". Previous test results are overwritten unless this
    directory is re-named. The testing script may attempt to install sample
    sequences in the DarkHorse MySQL database specified by the configuration
    file. If any of these sequences have been entered previously, the
    command line may show error messages of the following type, which can be
    safely ignored:
     
     DBD::mysql::db do failed: Duplicate entry'gi|110807720|ref|YP_691240.1|' for 
     key 1 at darkHorse-0.8/bin/installation_scripts/load_lineage_index_table.pl 
     line 184, <INPUT> line 7.     
     
----------------------------------------

USAGE INSTRUCTIONS

    A. General Procedure
        
        1. Run a BLAST query of all protein sequences in a genome (or metagenome)
           against the informative sequence set output by the installer program.
            
            Example Diamond BLASTP parameters:
	    	diamond makedb --in dh2_informative.fasta -d dh2_informative --threads 2
            	diamond blastp -d dh2_informative.dmnd -q genome.faa -a genome_dh2_informative.daa -t . --max-target-seqs 100 -p 4
            	diamond view -a genome_dh2_informative.daa -f tab  -o genome_dh2_informativep.m8
		
          Don't forget - if you update your DarkHorse database version, you must also 
	  update the informative fasta sequences used in the blast search. 
	                                
        2. Run darkhorse.pl from the unix command line, e.g.
            
            darkhorse.pl -c config file -t blast_tab_infile -e exclude list -g self.fasta 
               
            Output from the program will be written to a new directory, named
            using the following convention:
            
                calcs_processID_filt_setting 
                (e.g. calcs_13299_filt_0.1)           
            
      B. Customizing parameter settings
      
        DarkHorse provides a number of command line options, as described below.
        Options specified at the command line over-ride default parameters
        specified in the configuration file. File names should include complete
        paths, if not in current working directory.
        
        Usage: ./darkhorse.pl -c config file -t blast_tab_infile -e exclude list -g self.fasta 

          Required parameters
            -c  <configuration file name>   
            
            -g  <query sequence source file> 
                Fasta format query sequence file used for BLAST search
            
            -t  <tab delimited BLAST input file>
                BLAST results for query sequences versus Genbank nr database.
                
            -e  <exclude list>
                File containing a list of exclusion terms, one to a line. Can be any
                combination of words or NCBI taxonomy id numbers, in any order. BLAST
                matches containing any of these terms in are discarded. Used to exclude
                uninformative sequences (e.g. "contaminant") and/or matches to "self",
                defined at desired level of granularity (e.g. strain, species, or genus).
                       
          Optional parameters:
          
            -f  <filter threshold> [default = 0.1]
                Can be used to tune program sensitivity according to relative abundance
                of sequence data for species related to the test genome. Default value
                (0.1) requires all candidate BLAST matches for a given query to have
                bit scores in a range within 10% of the highest scoring non-self match
                for the same query. Lower values may increase sensitivity to potential
                horizontal gene transfer in query sets with few database relatives, but
                tend to increase false positive HGT predictions for query sets with
                abundant database relatives. 
            
            -n  <minimum number lineage terms> [default = 3]
                Sets the minimum number of lineage terms required for a match to be
                included in the taxonomy determination. For horizontal gene transfer
                studies, setting this value lower than 3 is not recommended, because
                it may cause matches with poorly characterized taxonomy to obscure more
                informative matches. For metagenomic studies, users willing to accept
                reduced specificity may want to try lowering this value to boost
                sensitivity.
		
    C. Supplementing a previously installed MySQL database 
    	
	Installed DarkHorse databases can be supplemented with additional sequences, including
	private, unpublished data, using the following script
	
		Darkhorse2/bin/installation_scripts/update_dh_database.pl
		
        In addition to a DarkHorse configuration file and fasta format file of input amino acid 
	sequences, this script also requires a tab-delimited index file with one line per sequence, 
	divided into three columns (note that columns 2 and 3 will be repeated on multiple lines).
	
		sequence_id_number    species_name     lineage
		
	The lineage column is a semicolon(;)-delimited lineage string, with taxonomic groups ranging
	from superkingdom to genus level listed in left-to-right order, e.g.
	
		Bacteria;Cyanobacteria;Gloeobacteria;Gloeobacterales;Gloeobacteraceae;Gloeobacter
    	
	If the species of interest is already included in NCBI's taxonomy database, it's lineage
	string can be found by searching the NCBI taxonomy website (https://www.ncbi.nlm.nih.gov/taxonomy).
	For novel organisms where no exact lineage is available from NCBI, users should just provide 
	their best guess, using lineages for closest available relatives as a model.
		
----------------------------------------

VERSION HISTORY
    
  Version 2.0 (beta) December 24, 2016
     Removed reliance on NCBI Genbank gi numbers. Provides a faster, more efficient 
     database installation process, and allows use of custom reference data sets, 
     including private and/or unpublished sequences.
  Version 1.5 October 2, 2013
     fixed bug causing LPI mis-calculation for organisms whose NCBI taxonomy lineage string
     contains two identical terms, for example Actinobacteria, which is both a phylum name 
     and a class name in the following example:
     Bacteria;Actinobacteria;Actinobacteria;Actinobacteridae;Bifidobacteriales;Bifidobacteriaceae;Gardnerella)
  Version 1.4 July 19, 2010
     added workaround for for NCBI BLAST bug, that gives inconsistent, truncated ID numbers
     for EMBL database entries if Genbank nr sequences are extracted from pre-formatted
     blast files using fastacmd, instead of downloaded directly from NCBI as as 
     unformatted fasta files.
  Version 1.3 December 14, 2009
     fixed incompatibility with some versions of Linux operating system (e.g. Centos and Ubuntu)
     due to syntax differences in unix sort command.
  Version 1.2 August 19, 2009
     fixed problem recognizing high-level self-id keywords occurring in lineage string only
     (e.g. kingdom, phylum, class, order)
  Version 1.1 April 6, 2009
     revised program installation script to handle larger data sets without memory errors
     added new requirement for NCBI Taxonomy reference file gi_taxid_prot.dmp
     removed requirement for NCBI accessory program fastacmd
  Version 1.0 October 24, 2008
     added error check for Genbank ID number mismatch between database and blast hits
  Version 0.9 May 9, 2008 
     fixed problem parsing extra space characters in NCBI blast output
     fixed problem parsing non-self matches if number exceeds 999
     added filters to to remove "subg", and "candidate division" from lineage
    Version 0.8 August 15, 2007
     initial beta distribution.
  	
---------------------------------------- 

 CONTACT INFORMATION
 
    Please send bug reports to: spodell@ucsd.edu
    Copyright 2017 S. Podell, University of California, San Diego. 
    All rights reserved.
  
----------------------------------------

