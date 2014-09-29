#################
#
# SCARecROW
# in Silico Consensus sequence Alignment from RECent maRkers Of Whole-population admixture
# Ryan Neff
# December 21, 2012
#
# Whole genome alignment tool which tolerates large differences between
# the individual's genome and the consensus sequence. Requires the 
# individual's subpopulation to be known, as well as 1000 Genomes data 
# of unmapped contigs in the individual's subpopulation.
#
# The program runs in two stages. First, the sub-population consensus
# sequence for the individual is created from 1000 Genomes' data, 
# accounting for SNPs, indels, and large structural variation using 
# available genotype data. Unmapped contigs from the 1000 Genomes project
# are assembled and by de novo sequencing and mapped to the human genome
# by admixture mapping of linkage-associated markers present in the 
# de novo contigs and available genotype data. Next, the sub-population
# consensus sequence is used to align the reads from the individual. 
#
# Requires BWA, SAMtools, Cortex_var, MALDsoft, Structure on the system path.
#
##################

#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

# declare variables for input
my $file_1000G_SNPs_vcf;
my $file_1000G_indels_vcf;
my $file_1000G_SVs_vcf;
my $file_consensus_fa;
my $file_reads_fastq;
my $file_1000G_samples_txt;
my $file_error;
my $threads;

my $path_unmapped;
my $path_cortexvar;
my $path_gatk;
my $path_structure;
my $path_maldsoft;
my $path_output;

# init program and validate program options
GetOptions( 'I=s' => \$file_reads_fastq,
			'R=s' => \$file_consensus_fa,
			'G_snps=s' => \$file_1000G_SNPs_vcf,
			'G_indels=s' => \$file_1000G_indels_vcf,
			'G_SVs=s' => \$file_1000G_SVs_vcf,
			'O=s' => \$path_output,
			'U=s' => \$path_unmapped;
			'E=s' => \$file_error,
			'T=d' => \$threads,
			'sample_file=s' => \$file_1000G_samples_txt,
			'cortex=s' => \$path_cortexvar,
			'gatk' => \$path_gatk,
			'structure_path=s' => \$path_structure,
			'maldsoft=s' => \$path_maldsoft,
		    'includeinfo' => \$includeinfo );

# set up log for output
open LOG, ">$file_error" || die "$0: Can't open error file!";
open ( STDERR, ">>$file_error" );
select( LOG );
$| = 1; # Turn on buffer autoflush for log output
select( STDOUT );

# step 1 - split 1000 Genomes genotypes by sample
&handleBash("java -Xmx4g -jar ". $path_gatk . "--unsafe ALL" . 
				" -R " . $file_consensus_fa . 
				" -T " . "SelectVariants" . 
				" -V " . $file_1000G_SNPs_vcf . 
				" -o " . $path_output . "/subpop_SNPs.vcf" . 
				" --sample_file " . $file_1000G_samples_txt );

&handleBash("java -Xmx4g -jar ". $path_gatk . "--unsafe ALL" . 
				" -R " . $file_consensus_fa . 
				" -T " . "SelectVariants" . 
				" -V " . $file_1000G_indels_vcf . 
				" -o " . $path_output . "/subpop_indels.vcf" . 
				" --sample_file " . $file_1000G_samples_txt );
				
&handleBash("java -Xmx4g -jar ". $path_gatk . "--unsafe ALL" . 
				" -R " . $file_consensus_fa . 
				" -T " . "SelectVariants" . 
				" -V " . $file_1000G_SVs_vcf . 
				" -o " . $path_output . "/subpop_SVs.vcf" . 
				" --sample_file " . $file_1000G_samples_txt );

&handleBash("java -Xmx4g -jar ". $path_gatk . "--unsafe ALL" . 
				" -R " . $file_consensus_fa . 
				" -T " . "CombineVariants" . 
				" -V:VCF " . $path_output . "/subpop_SNPs.vcf" .
				" -V:VCF " . $path_output . "/subpop_indels.vcf" .
				" -V:VCF " . $path_output . "/subpop_SVs.vcf" . 
				" -o " . $path_output . "/subpop_all_genotypes.vcf" );

# step 2 - filter samples by allele frequency >50%

&FilterAlleleFrequencies( $path_output . "/subpop_all_genotypes.vcf", 
						  0.5,   # minimum AF
						  20,	 # minimum number of samples (always)
						  $path_output . "/subpop_1000G_consensus.vcf");

# step 3 - assemble unmapped contigs from sample file
	# 3a. count the number of samples in path_unmapped
	my $samplecount = `wc -l \$file_1000G_samples_txt`;
	
	# 3b. compile cortex_var for that number of unmapped samples
	
	
	# 3c. run de novo mapping

# step 4 - place unmapped contigs on the current reference sequence
	# 4a. Structure v. 2.1
	# 4b. MALDsoft
	
# step 5 - create VCF file of unmapped contigs from MALDsoft

# step 6 - combine VCF of population changes and unmapped contigs

&handleBash("java -Xmx4g -jar ". $path_gatk . "--unsafe ALL" . 
				" -R " . $file_consensus_fa . 
				" -T " . "CombineVariants" . 
				" -V:VCF " . $path_output . "/subpop_all_genotypes.vcf"
				" -V:VCF " . $path_output . "/subpop_MALD_mapped_contigs.vcf"
				" -o " . $path_output . "/subpopulation_all_changes.vcf" );

# step 7 - make population-specific reference sequence

&handleBash("java -Xmx4g -jar ". $path_gatk . "--unsafe ALL" . 
				" -R " . $file_consensus_fa . 
				" -T " . "FastaAlternateReferenceMaker" . 
				" -V:VCF " . $path_output . "/subpopulation_all_changes.vcf" . 
				" -o " . $path_output . "/subpopulation_consensus.fa" );

# step 8 - align sample reads to new reference

&handleBash("bwa aln -R 100 -t " . $threads . " " .
			. $path_output . "/subpopulation_consensus.fa " . 
			. $file_reads_fastq . 
			" > " . $path_output . "/" . $file_reads_fastq . ".sai");

&handleBash("bwa samse -a 486 " .
			. $path_output . "/subpopulation_consensus.fa " . 
			. $path_output . "/" . $file_reads_fastq . ".sai " . 
			. $file_reads_fastq . 
			" > " . $path_output . "/" . $file_reads_fastq . ".sai.sam");

# step 9 - output alignment BAMs

&handleBash("samtools view -bS " . 
			$path_output . "/" . $file_reads_fastq . ".sai.sam" . 
			" > " $path_output . "/" . $file_reads_fastq . ".sai.sam.bam");

exit 0;


# ------------------------------- subroutines ------------------------

sub FilterAlleleFrequencies {
	# declare variables
	my $inputfilename = $_[0];
	my $afscutoff = $_[1];
	my $samplecutoff = $_[2];
	my $output_VCF = $_[3];
	my $inverse;
	my $line = "";

	# open input file
	open FILE, '<'.$inputfilename or die $!;
	open OUTVCF, '>>'.$output_VCF or die $!;
	print STDERR ("Starting analysis...\n");

	# stats
	my $numlinesprinted = 0;
	my $numlinesomitted = 0;
	my $numlinestotal = 0;

	# read in each line of VCF
	while (<FILE>){
		
		$numlinestotal++; # count lines 
		# print friendly updates
		my $rem = $numlinestotal % 10000;
		if($rem == 0){
			print STDERR ("Read $numlinestotal lines\n");
		}
		
		# start processing the line, determine type of line
		my $line  = $_; # collect line variable
		my @linearray = split(/\t/,$line,9);
		my $firstchar = (split "",$_,2)[0]; # split it into characters
		my $secondchar = (split "",$_,3)[1]; # split it into characters
		
		if($firstchar eq "#" && $secondchar eq "#"){ 
			# if we're in the header, just print the header
			print OUTVCF $_;
			next;
		}
		elsif($firstchar eq "#"){ # when we've reached the labels, reformat
			print OUTVCF ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
			next;
		}
		my $info = (split /\t/, $line)[7]; # split out INFO field
		my @afsfieldarray = (split ';', $info); # collect all properties
		
		# get AF value	
		my $afsfield = $afsfieldarray[1]; # get AF field
		my $afslabel = (split '=', $afsfield)[0]; # get AF label
		my $afsvalue = (split '=', $afsfield)[1]; # get AF value
		
		# in case we didn't get the AF field in the predicted spot, loop
		# over all values in the INFO section
		my $i = 0;
		while (($afslabel ne "AF") && ($i < $#afsfieldarray)){
			$afsfield = $afsfieldarray[$i];
			$afsvalue = (split '=', $afsfield)[1];
			$afslabel = (split '=', $afsfield)[0];
			$i++;
		}
		
		# get AC value
		my $acfield = $afsfieldarray[0]; # get AC field
		my $aclabel = (split '=', $acfield)[0]; # get AC label
		my $acvalue = (split '=', $acfield)[1]; # get AC value
		
		# in case we didn't get the AC field in the predicted spot, loop
		# over all values in the INFO section
		$i = 0;
		while (($aclabel ne "AC") && ($i < $#afsfieldarray)){
			$acfield = $afsfieldarray[$i];
			$acvalue = (split '=', $acfield)[1];
			$aclabel = (split '=', $acfield)[0];
			$i++;
		}
		
		# check for multiallelic - pass to subroutine to handle this case 
		if(! looks_like_number($afsvalue)){
			&handleMultiAllelic($line, $afsvalue, $acvalue);
			next;
		}
		
		# check it actually got set to AF, else print line information to STDERR
		if($afslabel ne 'AF'){
			print STDERR ("WARNING: Line $numlinestotal does not have an AF tag (omitted).\n");
			$numlinesomitted++;
		}
		elsif($aclabel ne 'AC'){
			print STDERR ("WARNING: Line $numlinestotal does not have an AC tag (omitted).\n");
			$numlinesomitted++;
		}
		else{
			# print all lines over the cutoff - under if inverse is set
			if (((($afscutoff < $afsvalue) && ! $inverse) || ($afscutoff > $afsvalue && $inverse)) && $acvalue >= $samplecutoff){
				print OUTVCF ("$linearray[0]\t$linearray[1]\t$linearray[2]\t$linearray[3]\t$linearray[4]\t$linearray[5]\t$linearray[6]\tZZ=Single;AF=$afsvalue;AC=$acvalue\n");
				$numlinesprinted++;
			}
			else{
				$numlinesomitted++;
			}
		}
	}
	# print run statistics
	print STDERR ("Filtering complete for $numlinestotal input lines.\n");
	print STDERR ("\tNumber of lines printed: $numlinesprinted\n");
	print STDERR ("\tNumber of lines discarded: $numlinesomitted\n");

	# end of program
	close(FILE);
	close(OUTVCF);
	return 0;
}

sub handleMultiAllelic {
	my $inputline = $_[0];
	my @linearray = split(/\t/,$inputline,9);
	my @allelefreq = split(',',$_[1]);
	my @acfreq = split(',',$_[2]);
	my @alleles = split(',', $linearray[4]);
	
	# loop over all alleles at that site
	my $int = 0;
	foreach(@alleles){
		if(((($allelefreq[$int] > $afscutoff) && ! $inverse) || ($afscutoff > $allelefreq[$int] && $inverse)) && $acfreq[$int] >= $samplecutoff){
			print OUTVCF ("$linearray[0]\t$linearray[1]\t$linearray[2]\t$linearray[3]\t$alleles[$int]\t$linearray[5]\t$linearray[6]\tZZ=Multiallelic;AF=$allelefreq[$int];AC=$acfreq[$int]\n");
			$numlinesprinted++;
		}
		else{
			$numlinesomitted++;
		}
		$int++;
	}
	return 0;
}

sub handleBash {
	my $command = $_[0];
	
	system($command);
	
	if ($? == -1) {
		print STDERR ("failed to execute: $!\n");
		exit 1;
	  }
	elsif ($? & 127) {
		print STDERR ("child died with signal %d, %s coredump\n"),
		($? & 127),  ($? & 128) ? 'with' : 'without';
		exit 1;
	 }
	else {
        print STDERR "child exited with value %d\n", $? >> 8;
        exit 1;
     }
	return 0;
}
