README.txt for DarkHorse version 2.0
May 2017
Author: Sheila Podell

----------------------------------------
CONTENTS

 - DESCRIPTION
   
 - SYSTEM REQUIREMENTS/DEPENDENCIES
 
 - INSTALLATION

 - USAGE INSTRUCTIONS

 - VERSION HISTORY

 - CONTACT INFORMATION

----------------------------------------
DESCRIPTION

   Darkhorse is an experimental program that defines phylogenetic
   relatedness of blastp hits for a set of proteins against a
   taxonomically diverse reference database (e.g. NCBI Genbank nr) using
   a taxonomically-weighted distance algorithm, with results expressed as
   a lineage probability index (LPI) score. The basic algorithm used to
   calculate LPI scores and its application in predicting horizontal gene
   transfer are described in the following publications. These
   publications should be cited in any work that uses the algorithm or
   results obtained using this software.
    
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
   access efficiency. A multi-processor CPU with at least 32 GB of RAM is recommended. 
     
   DarkHorse is a command-line only program, and requires the Unix OS. 
   It has been installed, tested, and run on multiple versions of Linux and 
   Macintosh OS X platforms. Additional software pre-requisites include: 
    
        PERL version 5.8.1 or later, with the following modules installed 
        (http://www.cpan.org/):
            DBI
            DBD::mysql
        MySQL Relational Database software version 5.5 or later (or MariaDB equivalent)
        	(https://dev.mysql.com/downloads/mysql/ or https://downloads.mariadb.org)     	
        Diamond software for performing protein sequence similarity searches:
       		(https://github.com/bbuchfink/diamond)   
            
    Prior to DarkHorse program installation, a database must be created 
    from within the MySQL program to accept DarkHorse input, e.g.
    
        mysql -u username -p -e "create database db_name;"
        
    The following public database information must be downloaded, de-compressed, 
    and locally available:
    
        Genbank nr database raw fasta file
            wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
	    gunzip nr.gz
        Genbank Taxonomy database 
            wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	    tar -xzvf taxdump.tar.gz
	Protein accession-tax_id index file
            wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
	    gunzip prot.accession2taxid.gz
        	      
    Users must be able to obtain protein BLAST (BLASTP) search data for query
    sequences versus a large, taxonomically informative set of reference sequences, 
    (e.g. Genbank nr) with output in the same tab-delimited format as 
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
       
   Installation is performed using a script that tests availability of
   pre-requisite programs and databases based on a user-customized configuration
   file. Note that some distributions may require program permissions to be
   re-set as executable after decompressing the downloaded files, e.g.
   
   	tar -xzvf DarkHorse-2.0_revXX.tar.gz
	chmod -R 755 DarkHorse-2.0_revXX
  
   A template configuration file is provided (templates/config_template), 
   which must be manually edited to include appropriate local information, 
   including paths to downloaded files, MySQL database name, MySQL user name,
   and MySQL user password. 
   
   If more than 16GB of RAM are available at the time of program execution,
   performance can be accelerated by increasing the "max_lines_per_packet"
   parameter in the configuration file. However, if this parameter is set too
   high, MySQL operations will fail due to insufficient memory (typical error
   message is something like DBD::mysql::db do failed: MySQL server has gone
   away.) Suggested settings based on total system memory (with no programs 
   other than DarkHorse being run concurrently), are as follows:
    
    	16 GB RAM  [max_lines_per_packet]=4000  (default value)
    	32 GB RAM  [max_lines_per_packet]=8000
    	64 GB RAM  [max_lines_per_packet]=16000
   
   The user-edited version of the configuration file is passed as a command 
   line parameter to the installation program:
            
         ./install_darkhorse2.pl -c config_filename        	
  
   Users of Oracle-MySQL version 5.6 or later should disable "strict mode" 
   in their MySQL program configuration prior to DarkHorse installation. 
   After installation is complete, this parameter can optionally be turned 
   back on, if desired. Strict mode is turned off by default in MariaDB 
   and earlier versions of MySQL. The reason why strict mode is a problem 
   is that database loading is performed in very large batches (packets), 
   containing thousands of sequences at a time. Every time there is an entry 
   with a field that fails to meet table format specifications (e.g. too many 
   characters) MySQL in strict mode discards the entire packet, including all 
   of the entries that were actually good. When strict mode is disabled, lines 
   exceeding the maximum number of characters are fixed by truncation, and then 
   everything else proceeds normally. 
   
   Note that the installation script may take a substantial amount of time to
   run (e.g many hours). In addition to creating and populating required MySQL
   tables, the installation script produces a log file containing detailed
   notes on any errors that may have occurred, and new fasta formatted
   sequence file of taxonomically informative sequences. The taxonomically
   informative fasta file will be named using the following convention:
            
         databasename_informative_ref_seqs.fasta    
   
   Users will need to format this file of protein sequences as a reference database for subsequent 
   blastp searches using Diamond or NCBI-blast software. 
   
----------------------------------------

USAGE INSTRUCTIONS

    A. General Procedure
        
       1. Run a BLAST query of all protein sequences in a genome (or metagenome)
          against the informative sequence set output from the installer program.
	  Recommended e-value setting for HGT determinations is 1e-5. For metagenomic
          binning or highly fragmented draft sequences, a less stringent e-value,
          (e.g. 1e-3) may improve sensitivity.
            
          Example Diamond BLASTP parameters:
             diamond makedb --in dh2_informative.fasta -d dh2_informative --threads 2
             diamond blastp -d dh2_informative.dmnd -q genome.faa -a genome_dh2.daa -e 1e-5 -t . --max-target-seqs 100 -p 4
             diamond view -a genome_dh2.daa -f tab  -o genome_dh2.m8          
                                            
          Don't forget - if you update your DarkHorse database version, you must also 
          update the informative fasta sequences used in the blast search. 
      
       2. Create an exclude_list file to remove self-matches from blast search results,
          if necessary. A template file is provided (templates/exclude_list_template),
          which can be used as is for classification of metagenomic sequences. For HGT
          analysis of individual genomes, an optional script has been provided to assist
          users in creating exclusion lists at different levels of taxonomic granularity
          (e.g. genes acquired since divergence from other strains, species, or genera
          sharing common ancestors with the test genome). However, these lists may need
	  to be manually edited and adjusted in cases where taxonomy identification 
	  numbers are duplicated or inconsistently applied in the NCBI taxonomy database.
	
	    Darkhorse2/bin/accessory_scripts/generate_dh_self_keywords.pl    	 	  
                              
        3. Run darkhorse.pl from the unix command line, e.g.
            
            darkhorse2.pl -c config_file -t blast_tab_infile -e exclude list -g self.fasta 
               
            Output from the program will be written to a new directory, named
            using the following convention:
            
                calcs_processID_filt_setting 
                (e.g. calcs_13299_filt_0.1)           
            
      B. Customizing parameter settings
      
        DarkHorse provides a number of command line options, as described below.
        Options specified at the command line over-ride default parameters
        specified in the configuration file. File names should include complete
        paths, if not in current working directory.
        
        Usage: ./darkhorse2.pl -c config_file -t blast_tab_infile -e exclude list -g self.fasta 

          Required parameters
            -c  <configuration file name>   
            
            -g  <query sequence source file> 
                Fasta format query sequence file used for BLAST search
            
            -t  <tab delimited BLAST input file>
                BLAST results for query sequences versus informative reference sequences
		obtained during program installation.
                
            -e  <exclude list>
                File containing a list of exclusion terms, one to a line. Can be any
                combination of words or NCBI taxonomy id numbers, in any order. BLAST
                matches containing any of these terms in are discarded. Used to exclude
                uninformative sequences (e.g. "contaminant") and/or matches to "self",
                defined at desired level of granularity (e.g. strain, species, or genus).
		An optional script has been provided to assist users in creating exclusion 
		lists for species that have been assigned multiple taxonomy ids in 
		the NCBI database:
	
	    		Darkhorse2/bin/accessory_scripts/generate_dh_self_keywords.pl
                       
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
		reduced specificity and phylum-level granularity for taxonomic assignment
		may want to try lowering this value to boost sensitivity, especially for
		poorly characterized taxa with few database representatives.
		
	   -b <blast filter> [default = 0.7]
	   	Values should be between 0-1 (Recommended 0.7). Sets the minimum alignment 
		length for blast matches as a percentage of the query and match sequence
		lengths (e.g. for a setting of 0.7, the blast alignment must cover at least
		70% of both query and database match sequences). Slows performance significantly 
		for large metagenomic data sets, but greatly improves accuracy and reduces 
		spurious matches to low-complexity regions. Set to value of zero to turn off.
     
     C. Supplementing a previously installed MySQL database 
    	
        Installed DarkHorse databases can be supplemented with additional sequences, 
        including private, unpublished data, using the following script:
	
		Darkhorse2/bin/installation_scripts/update_dh_database.pl
		
        This script will update a pre-existing MySQL DarkHorse database 
	and create a file of additional taxonomically informative protein sequences
	that can be used in BLASTP searches to supplement those previously 
	associated with the database. 
	
	In addition to a DarkHorse configuration file and fasta format file of
	input amino acid sequences, update_dh_database.pl also requires users to 
	prepare a tab-delimited text file (unix line endings please) with one 
	line per sequence. Each line should contain three columns (note that 
	columns 2 and 3 may be repeated on multiple lines).
	
			sequence_id_number    species_name     lineage
		
        The lineage column is a semicolon(;)-delimited lineage string, with taxonomic  
        group names ranging from superkingdom to genus listed in descending
	left-to-right order, e.g.
	
	    Bacteria;Cyanobacteria;Gloeobacteria;Gloeobacterales;Gloeobacteraceae;Gloeobacter
    	
        If the species of interest is already included in NCBI's taxonomy
 	database, it's lineage string can be found by searching the NCBI taxonomy
 	website (https://www.ncbi.nlm.nih.gov/taxonomy). For novel organisms where
 	no exact lineage is available from NCBI, users should just provide their
 	best guess, using the lineages of closest available relatives as a model. 
	
     D. Obtaining a taxonomically selected subset of reference fasta sequences
     
     	Taxonomically selected subsets can be extracted from the database-matched  
     	informative sequence fasta file generated by the DarkHorse installer  
	by using the following script:
	
	    Darkhorse2/bin/installation_scripts/subset_informative_fasta.pl
     	
     	This command-line perl script finds fasta format reference sequences 
     	from lineages containing one or more keywords entered by the user as 
     	a separate file (one keyword per line), and writes these selected 
     	sequences to a new file.
   
----------------------------------------

VERSION HISTORY
    
     Version 2.0 (beta) December 24, 2016
     	Removed reliance on NCBI Genbank gi numbers. Provides a faster, more efficient 
     	database installation process, and allows the use of custom reference data sets, 
     	including private and/or unpublished sequences. The installation software now 
	furnishes users with a verified, database-matched set of informative reference   
	sequences for use in subsequent BLAST searches, along with optional tools  
	for subdividing these sequences into smaller, taxonomically focused subsets.
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
