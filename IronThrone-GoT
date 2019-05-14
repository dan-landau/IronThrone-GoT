#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use threads;
use List::Util qw(sum);
use Compress::Zlib;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);


my $what_is_this_code = <<'__conquest_GoT__';

     ##### created by the Landau lab (@New York Genome Center)
     ##### for processing GoT amplicon data
__conquest_GoT__


# Global variables, they could be overwritten at init()
my $fastqR1    = undef;
my $fastqR2    = undef;
my $dupcut     = undef;
my $my_postP   = undef;
my $conf_file  = undef;
my $PRIMER_SEQ      = undef;
my $PRIMER_ALW      = undef;
my $GENERAL_SEQ     = undef;
my $GENERAL_ALW     = undef;
my $WT_SEQ          = undef;
my $WT_ALW          = undef;
my $pcr2_fw_wt_SEQ  = undef;
my $pcr2_fw_wt_ALW  = undef;
my $pcr2_fw_mut_SEQ = undef;
my $pcr2_fw_mut_ALW = undef;
my $pcr2_rv_wt_SEQ  = undef;
my $pcr2_rv_wt_ALW  = undef;
my $pcr2_rv_mut_SEQ = undef;
my $pcr2_rv_mut_ALW = undef;
my %MATRIX_HEADER;
my @WhitelistBarcodes;
my %WhitelistBarcodes_hash;
my %BARCODE_stringdist_hamm;


##################################################
# Parse paramters then assign provided global variables
##################################################

# default values
my $run_mode = "linear";
my $MMTCH = 0.2 ;
$my_postP = 0.99 ;
$dupcut = 1 ;
my $thread_max = 4;
my $len_barcode = 16;
my $len_umi = 10;
my $ALLOW_ERR = 0;
my $VERBOSE = 0 ;
my $keep_outs = 0 ;
my $out_dir = "./out";
my $samplename = "myGoT";
my $log_file="$samplename.log";
my $dir=dirname($0);
my $whitelist_file="$dir/barcodes10X/737K-august-2016.txt"; # provided from 10x


##################################################
# initiation of pipeline
##################################################

init();

open LOG, ">$out_dir/$log_file" || die $!; close LOG; # clear old file if exists
open LOG, ">>$out_dir/$log_file" || die $!; # open log file for append mode


##################################################
# main function
##################################################

pre_process()  ;
post_process() ;

##################################################

# Get option parameters
sub init {
	GetOptions (
		"run|r=s"      => \$run_mode,
		"fastqR1|f1=s" => \$fastqR1,
		"fastqR2|f2=s" => \$fastqR2,
		"config=s"     => \$conf_file,
		"whitelist=s"  => \$whitelist_file,
		"mmtch=s"      => \$MMTCH,
		"postP=s"      => \$my_postP,
		"dupcut=s"     => \$dupcut,
		"bclen=s"      => \$len_barcode,
		"umilen=s"     => \$len_umi,
		"sample=s"     => \$samplename,
		"outdir=s"     => \$out_dir,
		"log=s"        => \$log_file,
		"thread|t=s"   => \$thread_max,
		"keepouts=s"   => \$keep_outs,
		"vervose=s"    => \$VERBOSE,
		"help"         => sub { usage() },
		);

	if(!defined $fastqR1)  {usage();	exit 1	}
	if(!defined $fastqR2)  {usage();	exit 1	}
	if(!defined $conf_file){usage();	exit 1	}

	unless(-e $out_dir) {
		`mkdir -p $out_dir`;
	}
}


# Show option parameters
sub usage {
	print <<"USAGE";
     $what_is_this_code
     USAGE =================================================================================================
           IronThrone-GoT [options] --run linear --fastqR1 <in.R1.fastq> --fastqR2 <in.R2.fastq> \\
                                    --config <in.config> --sample <out.prefix> --outdir <out.path>

     REQUIRED INPUTS =======================================================================================
           -r/--run        run module for processing 'linear'/'circ' GoT (default: linear)
           -f1/--fastqR1   input R1.FASTQ FILE (input file can be in GZip format with .gz extension)
           -f2/--fastqR2   input R2.FASTQ FILE (input file can be in GZip format with .gz extension)
           -c/--config     input CONFIG FILE   (input file should be in tab-separated)
                           |
                           | .config file (for 'linear' GoT) should be prepared as follows
                           | to scan expected sequences at specific positions in amplicon reads
                           |
                           | (1st row)  primer.sequence	start.pos end.pos
                           | (2nd row)  shared.sequence	start.pos end.pos
                           | (3rd row)      WT.sequence	start.pos end.pos
                           | (4th row)     MUT.sequence	start.pos end.pos
                           |
                           | as an example, see testdata/linear_example.config
                           | or README file to prepare circ_example.config

     OPTIONS ===============================================================================================
           -m/--mmtch      allowed mismatch ratio to grep the expected sequences (default: 0.2)
           -p/--postP      cutoff for the posterior probability in the barcode replacement (default: 0.99)
           -d/--dupcut     cutoff for the total number of duplication (default: 1)
           -b/--bclen      length of barcode (default: 16)
           -u/--umilen     length of UMI (default: 10)
           -w/--whitelist  file for whitelisted barcodes (737K-august-2016.txt)
           -s/--sample     prefix for outputs (myGoT)
           -o/--outdir     path for outputs (./out)
           -l/--log        logfile name (myGoT.log)
           -t/--thread     number of threads to run in parallel (default: 4)
           -k/--keepouts   if set to 1, keeping intermediate files (default: 0)
           -v/--verbose    if set to 1, returning more logs (default: 0)
           -h/--help

USAGE

	exit ;
}


##################################################

# Parse necessary information from paired fastq files
sub pre_process {

	($fastqR1, $fastqR2) = un_gzip($fastqR1, $fastqR2); # unzip if .gz

	sed_repeat($fastqR1, "$out_dir/$samplename.BARCODE", 2, 4,  1, $len_barcode);
	sed_repeat($fastqR1, "$out_dir/$samplename.BARCODE_Q", 4, 4,  1, $len_barcode);
	sed_repeat($fastqR1, "$out_dir/$samplename.UMI", 2, 4,  $len_barcode+1, $len_barcode+$len_umi);
	sed_repeat($fastqR2, "$out_dir/$samplename.SEQ", 2, 4);
	sed_repeat($fastqR2, "$out_dir/$samplename.SEQ_Q", 4, 4);

	convert_phred_BQ("$out_dir/$samplename.BARCODE_Q" , "$out_dir/$samplename.BARCODE_Qnum");
	convert_numeric_BQ("$out_dir/$samplename.SEQ_Q", "$out_dir/$samplename.SEQ_Qnum");
	sum_of_BQ("$out_dir/$samplename.SEQ_Q");

	paste_cols(
		"$out_dir/$samplename.BARCODE",
		"$out_dir/$samplename.BARCODE_Qnum",
		"$out_dir/$samplename.UMI",
		"$out_dir/$samplename.SEQ",
		"$out_dir/$samplename.SEQ_Qnum",
		"$out_dir/$samplename.SEQ_Q.P_errors_SUM",

		"$out_dir/$samplename.SEQUENCE",
	);

  if ($run_mode eq "linear") {
    seqs_looking("$out_dir/$samplename.SEQUENCE") ;
  } elsif ($run_mode eq "circ") {
    seqs_looking2("$out_dir/$samplename.SEQUENCE") ;
  }

}


##################################################

# Perform additional processing on parsed amplicon data
sub post_process {
	##########################
	# load .config file
	##########################
  if ($run_mode eq "linear") {

    open CONF, $conf_file or die "[$conf_file] $!";
    my @conf_lines = <CONF>;
    close CONF;
    for(@conf_lines) { chomp; }

    $PRIMER_SEQ=(split "\t", $conf_lines[0])[0];
    $PRIMER_ALW=sprintf("%.0f", length($PRIMER_SEQ)*$MMTCH);

    $GENERAL_SEQ=(split "\t", $conf_lines[1])[0];
    $GENERAL_ALW=sprintf("%.0f", length($GENERAL_SEQ)*$MMTCH);

    $WT_SEQ=(split "\t", $conf_lines[2])[0];
    $WT_ALW=sprintf("%.0f", length($WT_SEQ)*$MMTCH);

  } elsif ($run_mode eq "circ") {

    open CONF, $conf_file or die "[$conf_file] $!";
    my @conf_lines = <CONF>;
    close CONF;
    for(@conf_lines) { chomp; }

    $PRIMER_SEQ=(split "\t", $conf_lines[0])[0];
    $PRIMER_ALW=sprintf("%.0f", length($PRIMER_SEQ)*$MMTCH);

    $GENERAL_SEQ=(split "\t", $conf_lines[1])[0];
    $GENERAL_ALW=sprintf("%.0f", length($GENERAL_SEQ)*$MMTCH);

    $WT_SEQ=(split "\t", $conf_lines[2])[0];
    $WT_ALW=sprintf("%.0f", length($WT_SEQ)*$MMTCH);

    $pcr2_fw_wt_SEQ=(split "\t", $conf_lines[4])[0];
    $pcr2_fw_wt_ALW=sprintf("%.0f", length($pcr2_fw_wt_SEQ)*$ALLOW_ERR);

    $pcr2_rv_wt_SEQ=(split "\t", $conf_lines[6])[0];
    $pcr2_rv_wt_ALW=sprintf("%.0f", length($pcr2_rv_wt_SEQ)*$ALLOW_ERR);
  }


	## Filter primed reads by mismatch allowance
	my @primed_looked;
	my @primed_looked_barcode;
	open LOOKED, "$out_dir/$samplename.looked" or die $!;
	my @LOOKED_RAW = <LOOKED>;
	close LOOKED;
	my $file_count_looked = @LOOKED_RAW;

  if ($run_mode eq "linear") {

		logging("seeking for mismatches from $file_count_looked lines | Primer<= $PRIMER_ALW | Shared<= $GENERAL_ALW ") if $VERBOSE;

		for my $line (@LOOKED_RAW) {
			my ($barcode, $mm_primer, $mm_share) = (split "\t", $line)[0, 3, 4];
			if ($mm_primer <= $PRIMER_ALW && $mm_share <= $GENERAL_ALW) {
				push @primed_looked, $line;
				push @primed_looked_barcode, $barcode;
			}
		}
		my $primed_looked_barcode_cnt = @primed_looked_barcode;
		logging("   -> $primed_looked_barcode_cnt found from $file_count_looked lines") if $VERBOSE;

	} elsif ($run_mode eq "circ") {

    logging("seeking for mismatches from $file_count_looked lines | Primer<= $PRIMER_ALW | Shared<= $GENERAL_ALW | PCR2FW_WT/PCR2FW_MUT<= $pcr2_fw_wt_ALW | PCR2RV_WT/PCR2RV_MUT<= $pcr2_rv_wt_ALW ") if $VERBOSE;

  	for my $line (@LOOKED_RAW) {
  		my ($barcode, $mm_primer, $mm_share, $mm_pcr2fw_wt, $mm_pcr2fw_mut, $mm_pcr2rv_wt, $mm_pcr2rv_mut) = (split "\t", $line)[0, 3, 4, 7, 8, 9, 10];
  		if ($mm_primer <= $PRIMER_ALW && $mm_share <= $GENERAL_ALW && ($mm_pcr2fw_wt <= $pcr2_fw_wt_ALW || $mm_pcr2fw_mut <= $pcr2_fw_wt_ALW) && ($mm_pcr2rv_wt <= $pcr2_rv_wt_ALW || $mm_pcr2rv_mut <= $pcr2_rv_wt_ALW)) {
  			push @primed_looked, $line;
  			push @primed_looked_barcode, $barcode;
  		}
  	}
  	my $primed_looked_barcode_cnt = @primed_looked_barcode;
  	logging("   -> $primed_looked_barcode_cnt found from $file_count_looked lines") if $VERBOSE;

  }


	# Load barcodes whitelist
	my $arr_tmp = file_read($whitelist_file);
	@WhitelistBarcodes = @$arr_tmp;

	my @obs_uniq_barcodes = unique(@primed_looked_barcode);

	$arr_tmp = intersect(\@WhitelistBarcodes, \@obs_uniq_barcodes);
	my @obs_barcodes_in_whitelist = @$arr_tmp;

	my @obs_barcodes_NOTin_whitelist = @obs_uniq_barcodes[ find_in_array("not in", \@obs_uniq_barcodes, \@obs_barcodes_in_whitelist) ];
	my @MTX_done = array_to_matrix( @primed_looked[ find_in_array("in", \@primed_looked_barcode, \@obs_barcodes_in_whitelist) ] );
	my @MTX_ND = array_to_matrix( @primed_looked[ find_in_array("in", \@primed_looked_barcode, \@obs_barcodes_NOTin_whitelist) ] );

	matrix_header_create("MTX_done", $MTX_done[0] );
	matrix_header_create("MTX_ND", $MTX_ND[0] );

	matrix_put_values(\@MTX_done, "MTX_done", "whitelist", "Y");
	matrix_put_values(\@MTX_ND, "MTX_ND", "whitelist", "N");

	my @MTX_sub = @MTX_ND;
	matrix_header_create("MTX_sub", undef, matrix_header_get("MTX_ND" ));
	matrix_header_rename("MTX_sub", "V1", "original_BC");
	matrix_put_values(\@MTX_sub, "MTX_sub", "replaced_BC", "NA");

	stringdist_hamm_pre_calc(\@MTX_sub); # pre-calculate stringdist_hamm() by multi-thread

  # Test barcode replacement probability
	my @tested;
	for (my $r=0; $r < @MTX_sub; $r++) {

		my $ori_BC = $MTX_sub[$r][0];
		my @BQ = split ";", $MTX_sub[$r][1];

		my @Hamm_1 = ();
		my $ref = $BARCODE_stringdist_hamm{ $ori_BC };
		if ($ref) {
		@Hamm_1 = @$ref;
		}

		if (@Hamm_1 < 1) {
			$MTX_sub[$r][ matrix_header_get_index("MTX_sub", "replaced_BC") ] = "No_Hamm_1";
		} else {
			my %dat_n_i;
			my @dat_n;
			my @dat_n_header = ("class", "bc", "qv", "At", "readN", "priorP", "p_edit", "likelihood", "postP");
			matrix_header_create("dat_n", undef, @dat_n_header );
			matrix_header_create("maxdat_n", undef, @dat_n_header );

			for (my $i=0; $i < @Hamm_1; $i++) {
				my $Hamm_1_i  = $Hamm_1[$i];

				my $diff_pos = str_diff_pos_1($ori_BC, $Hamm_1_i);

				my @tmp_data_test = ($Hamm_1_i);
				my @tmp_data_barcode = matrix_get_array_by_col(\@MTX_done, 0);
				my @inMTX_done_i =  matrix_get_matrix_by_row(\@MTX_done, find_in_array("in", \@tmp_data_barcode, \@tmp_data_test)  );
				my $NinMTX_done_i = @inMTX_done_i;

				if (@inMTX_done_i ==0) {
					$MTX_sub[$r][ matrix_header_get_index("MTX_sub", "replaced_BC") ] = "NO_in_SAVED";
					my @tmp = ($i, $Hamm_1_i, 0, $diff_pos, 0);
					push @dat_n, \@tmp;
				} else {
					my $Errs_n_i = (split ";", $MTX_sub[$r][1])[$diff_pos];
					my $ErrsEdit_n_i = min_value(33, $Errs_n_i);

					my @tmp = ($i, $Hamm_1_i, $ErrsEdit_n_i, $diff_pos, $NinMTX_done_i);
					push @dat_n, \@tmp;
				}
			}

			matrix_calc_priorP(\@dat_n, "dat_n");
			matrix_calc_p_edit(\@dat_n, "dat_n");
			matrix_calc_likelihood(\@dat_n, "dat_n");
			matrix_calc_postP(\@dat_n, "dat_n");

			my @Maxdat_n = matrix_calc_maxdat_n(\@dat_n, "dat_n", "postP");

			if ( matrix_sum_cols(\@dat_n, matrix_header_get_index("dat_n", "readN") ) == 0 && @Maxdat_n >= 2) {
				@Maxdat_n = ($Maxdat_n[0]);
				$MTX_sub[$r][ matrix_header_get_index("MTX_sub", "replaced_BC") ] = "NO_in_SAVED";
			}

			if ($Maxdat_n[0][ matrix_header_get_index("maxdat_n", "postP") ] ne "NA" && $Maxdat_n[0][ matrix_header_get_index("maxdat_n", "postP") ] >= $my_postP) {
				$MTX_sub[$r][ matrix_header_get_index("MTX_sub", "replaced_BC") ] = $Maxdat_n[0][ matrix_header_get_index("maxdat_n", "bc") ];
			} else {
				$MTX_sub[$r][ matrix_header_get_index("MTX_sub", "replaced_BC") ] = "Not_Sig";
			}
			if (@Maxdat_n == 2 && $Maxdat_n[0][ matrix_header_get_index("maxdat_n", "postP") ] == 0.5) {
				$MTX_sub[$r][ matrix_header_get_index("MTX_sub", "replaced_BC") ] = "Not_Sig";
				@Maxdat_n = ($Maxdat_n[0]);
			}
		}

		logging("processing row $r") if $VERBOSE;
		push @tested, $MTX_sub[$r];
	}

	my @MTX_changed = @tested;
	matrix_header_create("MTX_changed", undef, matrix_header_get("MTX_sub") );
	matrix_save(\@MTX_changed, "$out_dir/$samplename.MTX_changed.txt", "MTX_changed") if $keep_outs;


	# Replace barcode and merge with rescued reads
	my @MTX_changed_extracted = matrix_calc_mtx_extract(\@MTX_changed);
	matrix_header_create("MTX_changed_extracted", undef, matrix_header_get("MTX_sub") );
	matrix_calc_MTX_changed_extracted_put_bc(\@MTX_changed_extracted);

  my @MTX_RE;
  if ($run_mode eq "linear") {
    matrix_remove_col(\@MTX_changed_extracted, "MTX_changed_extracted", 13);
    matrix_header_rename("MTX_changed_extracted", "original_BC", "V1");
    matrix_header_rename("MTX_changed_extracted", "V13", "whitelist");

    @MTX_RE = (@MTX_done, @MTX_changed_extracted);
    matrix_header_create("MTX_RE", undef, matrix_header_get("MTX_changed_extracted") );
    matrix_header_rename("MTX_RE", "V6", "mismatch.WT");
    matrix_header_rename("MTX_RE", "V7", "mismatch.MUT");
    matrix_header_rename("MTX_RE", "V11", "BQ.in.WT");
    matrix_header_rename("MTX_RE", "V12", "BQ.in.MUT");

  } elsif ($run_mode eq "circ") {
    matrix_remove_col(\@MTX_changed_extracted, "MTX_changed_extracted", 21);
    matrix_header_rename("MTX_changed_extracted", "original_BC", "V1");
    matrix_header_rename("MTX_changed_extracted", "V21", "whitelist");

    @MTX_RE = (@MTX_done, @MTX_changed_extracted);
    matrix_header_create("MTX_RE", undef, matrix_header_get("MTX_changed_extracted") );
    matrix_header_rename("MTX_RE", "V6", "mismatch.WT");
    matrix_header_rename("MTX_RE", "V7", "mismatch.MUT");
    matrix_header_rename("MTX_RE", "V15", "BQ.in.WT");
    matrix_header_rename("MTX_RE", "V16", "BQ.in.MUT");

  }

	matrix_header_insert("MTX_RE", "BARCODE_UMI");
	matrix_calc_MTX_RE_BARCODE_UMI(\@MTX_RE, "MTX_RE");

	my @mtx_re_barcode_umi_list = matrix_get_array_by_col(\@MTX_RE, matrix_header_get_index("MTX_RE", "BARCODE_UMI"));
	my @uniq_BC_UMI = unique(@mtx_re_barcode_umi_list); #BARCODE_UMI
	my @MTX_RE_dedup;

	for (my $d=0; $d <@uniq_BC_UMI; $d++) {
		my $uniq_BC_UMI_d = $uniq_BC_UMI[$d];

		my @uniq_BC_UMI_d_arr_tmp = ($uniq_BC_UMI_d);
		my @MTX_RE_d = @MTX_RE[ find_in_array("in", \@mtx_re_barcode_umi_list, \@uniq_BC_UMI_d_arr_tmp) ];
		my $nrow_MTX_RE_d = @MTX_RE_d;
		matrix_header_create("MTX_RE_d", undef, matrix_header_get("MTX_RE") );
		matrix_put_values(\@MTX_RE_d, "MTX_RE_d", "total_dups", $nrow_MTX_RE_d);
		matrix_put_values(\@MTX_RE_d, "MTX_RE_d", "mismatch.WT.call", "");
		matrix_put_values(\@MTX_RE_d, "MTX_RE_d", "mismatch.MUT.call", "");
		matrix_put_values(\@MTX_RE_d, "MTX_RE_d", "mismatch.call", "amb");

		matrix_calc_MTX_RE_d_mismatch_WT_assign(\@MTX_RE_d, "MTX_RE_d");
		matrix_calc_MTX_RE_d_mismatch_MUT_assign(\@MTX_RE_d, "MTX_RE_d");
		matrix_calc_MTX_RE_d_mismatch_WT_call_assign(\@MTX_RE_d, "MTX_RE_d");
		matrix_calc_MTX_RE_d_mismatch_MUT_call_assign(\@MTX_RE_d, "MTX_RE_d");
		matrix_calc_MTX_RE_d_mismatch_amb_assign(\@MTX_RE_d, "MTX_RE_d");

		matrix_put_values(\@MTX_RE_d, "MTX_RE_d", "WT.inDups", scalar find_in_matrix("in", \@MTX_RE_d, "MTX_RE_d", "mismatch.call", "WT") );
		matrix_put_values(\@MTX_RE_d, "MTX_RE_d", "MUT.inDups", scalar find_in_matrix("in", \@MTX_RE_d, "MTX_RE_d", "mismatch.call", "MUT") );
		matrix_put_values(\@MTX_RE_d, "MTX_RE_d", "amb.inDups", scalar find_in_matrix("in", \@MTX_RE_d, "MTX_RE_d", "mismatch.call", "amb") );

		#1) Select highest frequency read type (i.e. mutant or wild-type) and obtain quality score
		my @MTX_wrap_d;
		my @MTX_wrap_d_BEST;
		matrix_header_create("MTX_wrap_d", undef, matrix_header_get("MTX_RE_d") );
		matrix_header_create("MTX_wrap_d_BEST", undef, matrix_header_get("MTX_RE_d") );

		if ($MTX_RE_d[0][ matrix_header_get_index("MTX_RE_d", "WT.inDups") ] > $MTX_RE_d[0][ matrix_header_get_index("MTX_RE_d", "MUT.inDups") ] ) {
			@MTX_wrap_d = @MTX_RE_d[ find_in_matrix("in", \@MTX_RE_d, "MTX_RE_d", "mismatch.call", "WT") ];
			@MTX_wrap_d_BEST = @MTX_RE_d[ find_in_matrix("in", \@MTX_RE_d, "MTX_RE_d", "BQ.in.WT", matrix_stat_col("min", \@MTX_RE_d, matrix_header_get_index("MTX_RE_d", "BQ.in.WT") ) ) ];

			if (@MTX_wrap_d_BEST ==1) {
			} else {
				#2) If the target quality of the remaining reads is equivalent then select the highest quality barcode read
				if (@MTX_wrap_d_BEST) {
					@MTX_wrap_d_BEST = @MTX_wrap_d_BEST[ matrix_calc_best_quality(\@MTX_wrap_d_BEST, matrix_header_get_index("MTX_wrap_d_BEST", "V2") ) ];
				}
			}
		}

		if ($MTX_RE_d[0][ matrix_header_get_index("MTX_RE_d", "WT.inDups") ] < $MTX_RE_d[0][ matrix_header_get_index("MTX_RE_d", "MUT.inDups") ]) {
			@MTX_wrap_d = @MTX_RE_d[ find_in_matrix("in", \@MTX_RE_d, "MTX_RE_d", "mismatch.call", "MUT") ];
			@MTX_wrap_d_BEST = @MTX_wrap_d[ find_in_matrix("in", \@MTX_wrap_d, "MTX_wrap_d", "BQ.in.WT", matrix_stat_col("min", \@MTX_wrap_d, matrix_header_get_index("MTX_wrap_d", "BQ.in.WT") ) ) ];

			if (@MTX_wrap_d_BEST ==1) {

			} else {
				if (@MTX_wrap_d_BEST) {
					@MTX_wrap_d_BEST = @MTX_wrap_d_BEST[ matrix_calc_best_quality(\@MTX_wrap_d_BEST, matrix_header_get_index("MTX_wrap_d_BEST", "V2") ) ];
				}
			}
		}

		#3) If both read types have equivalent frequencies, choose the type with the highest target quality score
		if ($MTX_RE_d[0][ matrix_header_get_index("MTX_RE_d", "WT.inDups") ] == $MTX_RE_d[0][ matrix_header_get_index("MTX_RE_d", "MUT.inDups") ]) {
			@MTX_wrap_d_BEST = @MTX_RE_d[ find_in_matrix("in", \@MTX_RE_d, "MTX_RE_d", "BQ.in.WT", matrix_stat_col("min", \@MTX_RE_d, matrix_header_get_index("MTX_RE_d", "BQ.in.WT") ) ) ];
			if (@MTX_wrap_d_BEST ==1) {

			} else {
				if (@MTX_wrap_d_BEST) {
					@MTX_wrap_d_BEST = @MTX_wrap_d_BEST[ matrix_calc_best_quality(\@MTX_wrap_d_BEST, matrix_header_get_index("MTX_wrap_d_BEST", "V2") ) ];
				}
			}
		}

    if ($run_mode eq "linear") {

      @MTX_wrap_d_BEST = matrix_remove_col(\@MTX_wrap_d_BEST, "MTX_wrap_d_BEST", 16);
  		@MTX_wrap_d_BEST = matrix_remove_col(\@MTX_wrap_d_BEST, "MTX_wrap_d_BEST", 15);
  		@MTX_wrap_d_BEST = matrix_remove_col(\@MTX_wrap_d_BEST, "MTX_wrap_d_BEST", 14);

    } elsif ($run_mode eq "circ") {

      @MTX_wrap_d_BEST = matrix_remove_col(\@MTX_wrap_d_BEST, "MTX_wrap_d_BEST", 25);
  		@MTX_wrap_d_BEST = matrix_remove_col(\@MTX_wrap_d_BEST, "MTX_wrap_d_BEST", 24);
  		@MTX_wrap_d_BEST = matrix_remove_col(\@MTX_wrap_d_BEST, "MTX_wrap_d_BEST", 23);

    }

		if (@MTX_wrap_d_BEST) {
			push @MTX_RE_dedup, @MTX_wrap_d_BEST;
		} else {
		}

		my $uniq_BC_UMI_len = @uniq_BC_UMI;
		my $d_no = $d + 1;
		logging("processing $d_no out of $uniq_BC_UMI_len") if $VERBOSE;
	}

	matrix_header_create("MTX_RE_dedup", undef, matrix_header_get("MTX_wrap_d_BEST") );

	matrix_save(\@MTX_RE_dedup, "$out_dir/$samplename.MTX_RE.dedup.txt","MTX_RE_dedup") if $keep_outs ;


	my @MTX_MERGE = @MTX_RE_dedup;
	matrix_header_create("MTX_MERGE", undef, matrix_header_get("MTX_RE_dedup") );


  if ($run_mode eq "linear") {

    @MTX_MERGE = matrix_remove_col(\@MTX_MERGE, "MTX_MERGE", 13);
    my @col_select = (1,13,3,15,16,17,14,  8,9,10,11,12,  4,5,6,7);
    for my $e (@col_select) { $e--; }

    @MTX_MERGE = matrix_get_col(\@MTX_MERGE, "MTX_MERGE", @col_select);

    matrix_header_create("MTX_MERGE", undef, ("BC","whitelist",
           "UMI","num.WT.in.dups","num.MUT.in.dups","num.amb.in.dups","call.in.dups",
           "avg.base_error.R2",
           "avg.base_error.primer","avg.base_error.shared","avg.base_error.WT","avg.base_error.MUT",
           "mismatch.primer","mismatch.shared","mismatch.WT","mismatch.MUT") );

  } elsif ($run_mode eq "circ") {

    @MTX_MERGE = matrix_remove_col(\@MTX_MERGE, "MTX_MERGE", 21);
    my @col_select = (1,21,3,23,24,25,22,  12,13,14,15,16, 4,5,6,7,  8,9,10,11);
    for my $e (@col_select) { $e--; }

    @MTX_MERGE = matrix_get_col(\@MTX_MERGE, "MTX_MERGE", @col_select);

    matrix_header_create("MTX_MERGE", undef, ("BC","whitelist",
           "UMI","num.WT.in.dups","num.MUT.in.dups","num.amb.in.dups","call.in.dups",
           "avg.base_error.R2",
           "avg.base_error.primer","avg.base_error.shared","avg.base_error.WT","avg.base_error.MUT",
           "mismatch.primer","mismatch.shared","mismatch.WT","mismatch.MUT",
           "mismatchPCR2FW_WT", "mismatchPCR2FW_MUT", "mismatchPCR2RV_WT", "mismatchPCR2RV_MUT") );

  }

	matrix_save(\@MTX_MERGE, "$out_dir/$samplename.MTX_MERGE.txt", "MTX_MERGE") if $keep_outs ;

	# Summarize accroding to barcode (for each cell)
	my @uniq_BC = unique( matrix_get_array_by_col(\@MTX_MERGE, matrix_header_get_index("MTX_MERGE", "BC") ) );
	my @summing_u;
	my @summing; # final result table;

	for (my $u=0; $u < @uniq_BC; $u++) {
		my $uniq_BC_u = $uniq_BC[$u] ;

		my @tmp_data_test = ($uniq_BC_u);
		my @tmp_data_barcode = matrix_get_array_by_col(\@MTX_MERGE,  matrix_header_get_index("MTX_MERGE", "BC") );
		my @MTX_MERGE_u =  matrix_get_matrix_by_row(\@MTX_MERGE, find_in_array("in", \@tmp_data_barcode, \@tmp_data_test)  );
		matrix_header_create("MTX_MERGE_u", undef, matrix_header_get("MTX_MERGE") );

		matrix_put_values(\@MTX_MERGE_u, "MTX_MERGE_u", "mismatch.WT.call", "" );
		matrix_put_values(\@MTX_MERGE_u, "MTX_MERGE_u", "mismatch.MUT.call", "" );
		matrix_put_values(\@MTX_MERGE_u, "MTX_MERGE_u", "mismatch.call", "" );

		matrix_calc_MTX_MERGE_u_WT_assign(\@MTX_MERGE_u, "MTX_MERGE_u");

		matrix_calc_MTX_MERGE_u_MUT_assign(\@MTX_MERGE_u, "MTX_MERGE_u");

		matrix_calc_MTX_MERGE_u_WT_call_assign(\@MTX_MERGE_u, "MTX_MERGE_u");

		matrix_calc_MTX_MERGE_u_MUT_call_assign(\@MTX_MERGE_u, "MTX_MERGE_u");

		matrix_calc_MTX_MERGE_u_amb_assign(\@MTX_MERGE_u, "MTX_MERGE_u");

		@MTX_MERGE_u = matrix_remove_col(\@MTX_MERGE_u, "MTX_MERGE_u", matrix_header_get_index("MTX_MERGE_u", "mismatch.MUT.call") );
		@MTX_MERGE_u = matrix_remove_col(\@MTX_MERGE_u, "MTX_MERGE_u", matrix_header_get_index("MTX_MERGE_u", "mismatch.WT.call") );

		my $arr_ref0 = $MTX_MERGE_u[0];
		my $arr_ref0_cnt = @$arr_ref0;
		my @arr_ref0_tmp = ( ("0") x  ($arr_ref0_cnt-1+3) );
		@summing_u = (\@arr_ref0_tmp);
		my @summing_u_hd = matrix_header_get("MTX_MERGE_u");
		@summing_u_hd = (@summing_u_hd[0..$#summing_u_hd - 1]);
		@summing_u_hd = (@summing_u_hd, "WT.calls", "MUT.calls", "amb.calls");
		matrix_header_create("summing_u", undef, @summing_u_hd );

		$summing_u[0][0] = $MTX_MERGE_u[0][0];

		my ($cnt_wt, $cnt_mut, $cnt_amb) = matrix_calc_get_mismatch_call(\@MTX_MERGE_u, "MTX_MERGE_u");
		matrix_put_values(\@summing_u, "summing_u", "WT.calls", $cnt_wt );
		matrix_put_values(\@summing_u, "summing_u", "MUT.calls", $cnt_mut );
		matrix_put_values(\@summing_u, "summing_u", "amb.calls", $cnt_amb );


		for (my $c =0; $c < matrix_header_get_ncol("MTX_MERGE_u") -1; $c++) {
			my $tmp = join ";", matrix_get_array_by_col(\@MTX_MERGE_u, $c);
			matrix_put_values(\@summing_u, "summing_u", matrix_header_get_name("summing_u", $c), $tmp);
		}
		push @summing, @summing_u;
	}

  # Filter by total number of duplicates
  my @final_sum;
  for my $line_ref (@summing) {
	  my @line_col = @$line_ref;

	  my $dup_WT = sum(split(";", $line_col[3]));
	  my $dup_MUT= sum(split(";", $line_col[4]));
	  my $dup_amb= sum(split(";", $line_col[5]));

      push @final_sum, $line_ref if($dup_WT + $dup_MUT + $dup_amb >= $dupcut);
  }
	unless ($keep_outs) {`rm $out_dir/$samplename.*`}
	matrix_save(\@final_sum, "$out_dir/$samplename.summTable.txt", "summing_u");
	print "______________________________END____________________________\n";

}



###############
# sub-functions
###############
sub find_mm_1 {
	my ($ref, $str, $str_exclude) = @_;
	my $mm_cnt;
	my @out;
	for my $e (@$ref) {
		next if $e eq $str_exclude;
		$mm_cnt = ($str ^ $e) =~ tr/\001-\255//; # XOR
		push @out, $e if $mm_cnt == 1;
	}
	@out;
}


sub sum_of_BQ {
	my ($file_in) = @_;

	if (!-e $file_in || -s $file_in <= 0) {
		logging("$file_in does not exist or 0 byte!") if $VERBOSE;
		exit;
	}

	my %phred_scores = load_phred_scores();

	open (FILE, $file_in) or die $!;
	my @MY_input =<FILE>;
	close FILE;

	open OUT, ">$file_in.P_errors_SUM" or die $!;
	my $cnt = 1;
	for my $line (@MY_input) {
	  chomp $line;
	  my $sum = 0;
	  my @bases = split(//,$line);
		for (my $i=0; $i < @bases; $i++){
			$sum += $phred_scores{ $bases[$i] }->{'p_err'};
	   }
	  print OUT "$sum\n";
	  $cnt++;
	}
	close OUT;
}


sub convert_phred_BQ {
	my ($file_in, $file_out) = @_;

	if (!-e $file_in || -s $file_in <= 0) {
		logging("$file_in does not exist or 0 byte!");
		exit;
	}

	my %phred_scores = load_phred_scores();

	open OUT, ">$file_out" || die $!;
	open FILE, $file_in || die $!;
	while(<FILE>){
		chomp;
		my @base = split(//,$_);
		for (my $i=0; $i <= (length($_)-1); $i++){
			print OUT $phred_scores{ $base[$i] }->{'q_score'}.";";
		}
		print OUT "\n";
	}

}


sub convert_numeric_BQ {
	my ($file_in, $file_out) = @_;

	if (!-e $file_in || -s $file_in <= 0) {
		logging("$file_in does not exist or 0 byte!");
		exit;
	}

	my %phred_scores = load_phred_scores();

	open OUT, ">$file_out" || die $!;
	open FILE, $file_in || die $!;
	while(<FILE>){
		chomp;
		my @base = split(//,$_);
		for (my $i=0; $i <= (length($_)-1); $i++){
			print OUT $phred_scores{ $base[$i] }->{'p_err'}.";";
		}
		print OUT "\n";
	}

	logging("[$file_out] created");
}


sub seqs_looking {
	my ($file_in) = @_;

	if (!-e $file_in || -s $file_in <= 0) {
		logging("$file_in does not exist or 0 byte!");
		exit;
	}

	open CONF, $conf_file or die "[$file_in] $!";
	my @conf_lines = <CONF>;
	close CONF;
	for(@conf_lines) { chomp; }

	my $PRIMER_SEQ=(split "\t", $conf_lines[0])[0];
	my $PRIMER_LEN=length($PRIMER_SEQ);
	my $PRIMER_ALW=sprintf("%.0f", $PRIMER_LEN*$ALLOW_ERR);
	my $PRIMER_START=(split "\t", $conf_lines[0])[1];
	my $PRIMER_END=(split "\t", $conf_lines[0])[2];
	my $PRIMER_MMTCH=sprintf("%.0f", $PRIMER_LEN*$MMTCH);

	my $GENERAL_SEQ=(split "\t", $conf_lines[1])[0];
	my $GENERAL_START=(split "\t", $conf_lines[1])[1];
	my $GENERAL_END=(split "\t", $conf_lines[1])[2];
	my $GENERAL_ALW=sprintf("%.0f", length($GENERAL_SEQ)*$ALLOW_ERR);

	my $WT_SEQ=(split "\t", $conf_lines[2])[0];
	my $WT_START=(split "\t", $conf_lines[2])[1];
	my $WT_END=(split "\t", $conf_lines[2])[2];
	my $WT_ALW=sprintf("%.0f", length($WT_SEQ)*$ALLOW_ERR);

	my $MUT_SEQ=(split "\t", $conf_lines[3])[0];
	my $MUT_START=(split "\t", $conf_lines[3])[1];
	my $MUT_END=(split "\t", $conf_lines[3])[2];
	my $MUT_ALW=sprintf("%.0f", length($MUT_SEQ)*$ALLOW_ERR);
	##########################

	my $SUM_SEQ_ref = file_read($file_in);
	my @SUM_SEQ = @$SUM_SEQ_ref;

	my $cnt_tmp = 0;
	open (my $look_seq, '>', "$out_dir/$samplename.looked");
	my $line_cnt = 0;

	for (my $i=0; $i < @SUM_SEQ; $i++) {
		$line_cnt++;
		my $line_raw = $SUM_SEQ[$i];
		chomp $line_raw;
		next if $line_raw =~ /^\s*$/;

		my @line_col = split("\t", $line_raw);
		my $line = $line_col[3];

		my $PRIMER_FOUND;
		my $PRIMER_STARTz;
		my $GENERAL_STARTz;
		my $GENERAL_ENDz;
		my $WT_STARTz;
		my $MUT_STARTz;
		my $MUT_ENDz;
		my $WT_ENDz;
		my $PRIMER_ENDz;

		$PRIMER_FOUND = "";
		for my $loc (0..3) { # consider four shifted locations
			my $seq_in_line = substr($line, $PRIMER_START-1+$loc, length($PRIMER_SEQ)+$loc);
			my $PRIMER_FOUND_check = mismatch_check($PRIMER_SEQ, $seq_in_line, $PRIMER_ALW);

			if ($PRIMER_FOUND_check <= $PRIMER_MMTCH) {
				$PRIMER_FOUND = $PRIMER_FOUND_check ;
				$GENERAL_STARTz = $GENERAL_START+$loc ; $GENERAL_ENDz = $GENERAL_END+$loc ;
				$WT_STARTz = $WT_START+$loc ; $WT_ENDz = $WT_END+$loc ;
				$MUT_STARTz = $MUT_START+$loc ; $MUT_ENDz = $MUT_END+$loc ;
				last;
			}
		}
		if ($PRIMER_FOUND eq "") {
			$PRIMER_FOUND =  $PRIMER_LEN;
			$GENERAL_STARTz = $GENERAL_START ; $GENERAL_ENDz = $GENERAL_END ;
			$WT_STARTz = $WT_START ; $WT_ENDz = $WT_END ;
			$MUT_STARTz = $MUT_START ; $MUT_ENDz = $MUT_END ;
		}

	    my $GENERAL_FOUND = mismatch_check($GENERAL_SEQ, substr($line, $GENERAL_STARTz-1, length($GENERAL_SEQ)), 0);
	    my $WT_FOUND = mismatch_check($WT_SEQ, substr($line, $WT_STARTz-1, length($WT_SEQ)), 0);
	    my $MUT_FOUND = mismatch_check($MUT_SEQ, substr($line, $MUT_STARTz-1, length($MUT_SEQ)), 0);

if (!$line or !$line_col[5] or !length($GENERAL_SEQ) or
	!sum((split(";", $line_col[4]))[$GENERAL_STARTz-1..$GENERAL_ENDz-1])) {
	print "line: $line\n";
	exit;
}

	    my $avg_R2_BQ = $line_col[5]/length($line);
	    my @R2_BQ = split(";", $line_col[4]);

	    my @BQ_primer = @R2_BQ[1-1..$PRIMER_LEN-1];
	    my $BQ_primer_avg = sum(@BQ_primer)/$PRIMER_LEN ;

	    my @BQ_general = @R2_BQ[$GENERAL_STARTz-1..$GENERAL_ENDz-1];
	    my $BQ_general_avg = sum(@BQ_general)/length($GENERAL_SEQ) ;

	    my @BQ_WT = @R2_BQ[$WT_STARTz-1..$WT_ENDz-1];
	    my $BQ_WT_avg = sum(@BQ_WT)/length($WT_SEQ) ;

	    my @BQ_MUT = @R2_BQ[$MUT_STARTz-1..$MUT_ENDz-1];
	    my $BQ_MUT_avg = sum(@BQ_MUT)/length($MUT_SEQ) ;

		print $look_seq "$line_col[0]\t$line_col[1]\t$line_col[2]\t$PRIMER_FOUND\t$GENERAL_FOUND\t$WT_FOUND\t$MUT_FOUND\t$avg_R2_BQ\t$BQ_primer_avg\t$BQ_general_avg\t$BQ_WT_avg\t$BQ_MUT_avg\n";
  }
  close $look_seq;
}

sub seqs_looking2 {
	my ($file_in) = @_;

	if (!-e $file_in || -s $file_in <= 0) {
		logging("$file_in does not exist or 0 byte!");
		exit;
	}

	open CONF, $conf_file or die "[$file_in] $!";
	my @conf_lines = <CONF>;
	close CONF;
	for(@conf_lines) { chomp; }

	my $PRIMER_SEQ=(split "\t", $conf_lines[0])[0];
	my $PRIMER_LEN=length($PRIMER_SEQ);
	my $PRIMER_ALW=sprintf("%.0f", $PRIMER_LEN*$ALLOW_ERR);
	my $PRIMER_START=(split "\t", $conf_lines[0])[1];
	my $PRIMER_END=(split "\t", $conf_lines[0])[2];
	my $PRIMER_MMTCH=sprintf("%.0f", $PRIMER_LEN*$MMTCH);

	my $GENERAL_SEQ=(split "\t", $conf_lines[1])[0];
	my $GENERAL_START=(split "\t", $conf_lines[1])[1];
	my $GENERAL_END=(split "\t", $conf_lines[1])[2];
	my $GENERAL_ALW=sprintf("%.0f", length($GENERAL_SEQ)*$ALLOW_ERR);

	my $WT_SEQ=(split "\t", $conf_lines[2])[0];
	my $WT_START=(split "\t", $conf_lines[2])[1];
	my $WT_END=(split "\t", $conf_lines[2])[2];
	my $WT_ALW=sprintf("%.0f", length($WT_SEQ)*$ALLOW_ERR);

	my $MUT_SEQ=(split "\t", $conf_lines[3])[0];
	my $MUT_START=(split "\t", $conf_lines[3])[1];
	my $MUT_END=(split "\t", $conf_lines[3])[2];
	my $MUT_ALW=sprintf("%.0f", length($MUT_SEQ)*$ALLOW_ERR);

	my $pcr2_fw_wt_SEQ=(split "\t", $conf_lines[4])[0];
	my $pcr2_fw_wt_START=(split "\t", $conf_lines[4])[1];
	my $pcr2_fw_wt_END=(split "\t", $conf_lines[4])[2];
	my $pcr2_fw_wt_ALW=sprintf("%.0f", length($pcr2_fw_wt_SEQ)*$ALLOW_ERR);

	my $pcr2_fw_mut_SEQ=(split "\t", $conf_lines[5])[0];
	my $pcr2_fw_mut_START=(split "\t", $conf_lines[5])[1];
	my $pcr2_fw_mut_END=(split "\t", $conf_lines[5])[2];
	my $pcr2_fw_mut_ALW=sprintf("%.0f", length($pcr2_fw_mut_SEQ)*$ALLOW_ERR);

	my $pcr2_rv_wt_SEQ=(split "\t", $conf_lines[6])[0];
	my $pcr2_rv_wt_START=(split "\t", $conf_lines[6])[1];
	my $pcr2_rv_wt_END=(split "\t", $conf_lines[6])[2];
	my $pcr2_rv_wt_ALW=sprintf("%.0f", length($pcr2_rv_wt_SEQ)*$ALLOW_ERR);

	my $pcr2_rv_mut_SEQ=(split "\t", $conf_lines[7])[0];
	my $pcr2_rv_mut_START=(split "\t", $conf_lines[7])[1];
	my $pcr2_rv_mut_END=(split "\t", $conf_lines[7])[2];
	my $pcr2_rv_mut_ALW=sprintf("%.0f", length($pcr2_rv_mut_SEQ)*$ALLOW_ERR);


	##########################

	my $SUM_SEQ_ref = file_read($file_in);
	my @SUM_SEQ = @$SUM_SEQ_ref;

	my $cnt_tmp = 0;
	open (my $look_seq, '>', "$out_dir/$samplename.looked");
	my $line_cnt = 0;

	for (my $i=0; $i < @SUM_SEQ; $i++) {
		$line_cnt++;
		my $line_raw = $SUM_SEQ[$i];
		chomp $line_raw;
		next if $line_raw =~ /^\s*$/;

		my @line_col = split("\t", $line_raw);
		my $line = $line_col[3];

		my $PRIMER_FOUND;
		my $PRIMER_STARTz; my $PRIMER_ENDz;
		my $GENERAL_STARTz; my $GENERAL_ENDz;
		my $WT_STARTz; my $WT_ENDz;
		my $MUT_STARTz; my $MUT_ENDz;
		my $pcr2_fw_wt_STARTz ; my $pcr2_fw_wt_ENDz;
		my $pcr2_fw_mut_STARTz; my $pcr2_fw_mut_ENDz;
		my $pcr2_rv_wt_STARTz ; my $pcr2_rv_wt_ENDz;
		my $pcr2_rv_mut_STARTz; my $pcr2_rv_mut_ENDz;


		$PRIMER_FOUND = "";
		for my $loc (0..3) { # Consider four shifted locations in the case of using staggered primer sets
			my $seq_in_line = substr($line, $PRIMER_START-1+$loc, length($PRIMER_SEQ)+$loc);
			my $PRIMER_FOUND_check = mismatch_check($PRIMER_SEQ, $seq_in_line, $PRIMER_ALW);

			if ($PRIMER_FOUND_check <= $PRIMER_MMTCH) {
				$PRIMER_FOUND = $PRIMER_FOUND_check ;
				$GENERAL_STARTz = $GENERAL_START+$loc ; $GENERAL_ENDz = $GENERAL_END+$loc ;
				$WT_STARTz = $WT_START+$loc ; $WT_ENDz = $WT_END+$loc ;
				$MUT_STARTz = $MUT_START+$loc ; $MUT_ENDz = $MUT_END+$loc ;
				$pcr2_fw_wt_STARTz = $pcr2_fw_wt_START+$loc ; $pcr2_fw_wt_ENDz = $pcr2_fw_wt_END+$loc ;
				$pcr2_fw_mut_STARTz = $pcr2_fw_mut_START+$loc ; $pcr2_fw_mut_ENDz = $pcr2_fw_mut_END+$loc ;
				$pcr2_rv_wt_STARTz = $pcr2_rv_wt_START+$loc ; $pcr2_rv_wt_ENDz = $pcr2_rv_wt_END+$loc ;
				$pcr2_rv_mut_STARTz = $pcr2_rv_mut_START+$loc ; $pcr2_rv_mut_ENDz = $pcr2_rv_mut_END+$loc ;
				last;
			}
		}
		if ($PRIMER_FOUND eq "") {
			$PRIMER_FOUND =  $PRIMER_LEN;
			$GENERAL_STARTz = $GENERAL_START ; $GENERAL_ENDz = $GENERAL_END ;
			$WT_STARTz = $WT_START ; $WT_ENDz = $WT_END ;
			$MUT_STARTz = $MUT_START ; $MUT_ENDz = $MUT_END ;
			$pcr2_fw_wt_STARTz = $pcr2_fw_wt_START ; $pcr2_fw_wt_ENDz = $pcr2_fw_wt_END ;
			$pcr2_fw_mut_STARTz = $pcr2_fw_mut_START ; $pcr2_fw_mut_ENDz = $pcr2_fw_mut_END ;
			$pcr2_rv_wt_STARTz = $pcr2_rv_wt_START ; $pcr2_rv_wt_ENDz = $pcr2_rv_wt_END ;
			$pcr2_rv_mut_STARTz = $pcr2_rv_mut_START ; $pcr2_rv_mut_ENDz = $pcr2_rv_mut_END ;
		}

	    my $GENERAL_FOUND = mismatch_check($GENERAL_SEQ, substr($line, $GENERAL_STARTz-1, length($GENERAL_SEQ)), 0);
	    my $WT_FOUND = mismatch_check($WT_SEQ, substr($line, $WT_STARTz-1, length($WT_SEQ)), 0);
	    my $MUT_FOUND = mismatch_check($MUT_SEQ, substr($line, $MUT_STARTz-1, length($MUT_SEQ)), 0);
	    my $pcr2_fw_wt_FOUND = mismatch_check($pcr2_fw_wt_SEQ, substr($line, $pcr2_fw_wt_STARTz-1, length($pcr2_fw_wt_SEQ)), 0);
	    my $pcr2_fw_mut_FOUND = mismatch_check($pcr2_fw_mut_SEQ, substr($line, $pcr2_fw_mut_STARTz-1, length($pcr2_fw_mut_SEQ)), 0);
	    my $pcr2_rv_wt_FOUND = mismatch_check($pcr2_rv_wt_SEQ, substr($line, $pcr2_rv_wt_STARTz-1, length($pcr2_rv_wt_SEQ)), 0);
	    my $pcr2_rv_mut_FOUND = mismatch_check($pcr2_rv_mut_SEQ, substr($line, $pcr2_rv_mut_STARTz-1, length($pcr2_rv_mut_SEQ)), 0);


if (!$line or !$line_col[5] or !length($GENERAL_SEQ) or
	!sum((split(";", $line_col[4]))[$GENERAL_STARTz-1..$GENERAL_ENDz-1])) {
	print "line: $line\n";
	exit;
}

	    my $avg_R2_BQ = $line_col[5]/length($line);
	    my @R2_BQ = split(";", $line_col[4]);

	    my @BQ_primer = @R2_BQ[1-1..$PRIMER_LEN-1];
	    my $BQ_primer_avg = sum(@BQ_primer)/$PRIMER_LEN ;

	    my @BQ_general = @R2_BQ[$GENERAL_STARTz-1..$GENERAL_ENDz-1];
	    my $BQ_general_avg = sum(@BQ_general)/length($GENERAL_SEQ) ;

	    my @BQ_WT = @R2_BQ[$WT_STARTz-1..$WT_ENDz-1];
	    my $BQ_WT_avg = sum(@BQ_WT)/length($WT_SEQ) ;

	    my @BQ_MUT = @R2_BQ[$MUT_STARTz-1..$MUT_ENDz-1];
	    my $BQ_MUT_avg = sum(@BQ_MUT)/length($MUT_SEQ) ;

	    my @BQ_pcr2_fw_wt = @R2_BQ[$pcr2_fw_wt_STARTz-1..$pcr2_fw_wt_ENDz-1];
	    my $BQ_pcr2_fw_wt_avg = sum(@BQ_pcr2_fw_wt)/length($pcr2_fw_wt_SEQ) ;

	    my @BQ_pcr2_fw_mut = @R2_BQ[$pcr2_fw_mut_STARTz-1..$pcr2_fw_mut_ENDz-1];
	    my $BQ_pcr2_fw_mut_avg = sum(@BQ_pcr2_fw_mut)/length($pcr2_fw_mut_SEQ) ;

	    my @BQ_pcr2_rv_wt = @R2_BQ[$pcr2_rv_wt_STARTz-1..$pcr2_rv_wt_ENDz-1];
	    my $BQ_pcr2_rv_wt_avg = sum(@BQ_pcr2_rv_wt)/length($pcr2_rv_wt_SEQ) ;

	    my @BQ_pcr2_rv_mut = @R2_BQ[$pcr2_rv_mut_STARTz-1..$pcr2_rv_mut_ENDz-1];
	    my $BQ_pcr2_rv_mut_avg = sum(@BQ_pcr2_rv_mut)/length($pcr2_rv_mut_SEQ) ;

		print $look_seq "$line_col[0]\t$line_col[1]\t$line_col[2]\t$PRIMER_FOUND\t$GENERAL_FOUND\t$WT_FOUND\t$MUT_FOUND\t$pcr2_fw_wt_FOUND\t$pcr2_fw_mut_FOUND\t$pcr2_rv_wt_FOUND\t$pcr2_rv_mut_FOUND\t$avg_R2_BQ\t$BQ_primer_avg\t$BQ_general_avg\t$BQ_WT_avg\t$BQ_MUT_avg\t$BQ_pcr2_fw_wt_avg\t$BQ_pcr2_fw_mut_avg\t$BQ_pcr2_rv_wt_avg\t$BQ_pcr2_rv_wt_avg\n";

  }
  close $look_seq;
}


sub mismatch_check {
	my ($str_ref, $str_test, $mismatch_allow_count) = @_;
	chomp $str_ref;
	chomp $str_test;

	my $len = length $str_test;
	my $mismatched = 0;

	for (my $i =0; $i < length($str_ref); $i++) {
		my $letter_ref = substr($str_ref, $i, 1);
		my $letter_test = substr($str_test, $i, 1);
		if ($letter_ref ne $letter_test) { $mismatched++; }
	}

	if ($mismatched > $mismatch_allow_count) {
		return "$mismatched";
	} else {
		return "$mismatched";
	}
}


sub logging {
	my $msg = shift;
	my $time = localtime;
	print LOG "[$time] $msg\n";
	print "[$time] $msg\n";
}


sub file_read {
	my ($file) = @_;
	open F, $file or die $!;
	my @tmp = <F>;
	close F;
	my $file_size = -s $file;
	my $lines = @tmp;
	logging("File read: $file (size: $file_size,   lines: $lines)"); # if $VERBOSE;

	return \@tmp;
}


sub unique {
	my (@arr) = @_;
	my %seen;
	my @out;
	my $cnt_in = @arr;
	for my $e (@arr) {
		if (!exists $seen{$e}) { push @out, $e; $seen{$e}++; }
	}
	my $cnt_out = @out;
	logging("Unique done: $cnt_in -> $cnt_out") if $VERBOSE;
	return @out;
}


sub intersect {
	my ($arr_ref1, $arr_ref2) = @_;
	my %seen;
	my @out;
	my $cnt_in1 = @$arr_ref1;
	my $cnt_in2 = @$arr_ref2;

	for my $e (@$arr_ref2) { chomp $e; $seen{$e}++; }

	for my $e (@$arr_ref1) {
		chomp $e;
		if (exists $seen{$e}) { push @out, $e; }
	}
	my $cnt_out = @out;
	logging("Intersect done: input count($cnt_in1, $cnt_in2),  intersect count: $cnt_out") if $VERBOSE;

	return \@out;
}


sub find_in_array {
	my ($mode, $arr_ref1, $arr_ref2) = @_;
	my %seen;
	my @out;
	my $cnt_in1 = @$arr_ref1;
	my $cnt_in2 = @$arr_ref2;

	for my $e (@$arr_ref2) { chomp $e; $seen{$e}++; }

	for (my $i = 0; $i < @$arr_ref1; $i++) {
		my $e = $$arr_ref1[$i];
		chomp $e;

		if ($mode eq "in") {
			if (exists $seen{$e}) { push @out, $i; }
		} elsif ($mode eq "not in") {
			if (!exists $seen{$e}) { push @out, $i; }
		}
	}
	my $cnt_out = @out;
	logging("find $mode done: input count($cnt_in1, $cnt_in2),  found count: $cnt_out") if $VERBOSE;

	return @out;
}


sub matrix_get_array_by_col {
	my ($mtx_ref, $idx) = @_;
	my @out;
	for my $e_ref (@$mtx_ref) {
		push @out, ${$e_ref}[$idx];
	}
	@out;
}


sub matrix_get_matrix_by_row {
	my ($mtx_ref, @rows) = @_;
	my @out;
	my %rows_seen;
	for my $e (@rows) { $rows_seen{$e}++; }

	for (my $i=0; $i <@$mtx_ref; $i++) {
		if ($rows_seen{$i}) {
			push @out, $$mtx_ref[$i];
		}
	}
	@out;
}


sub matrix_get_col {
	my ($mat_ref, $mat_name, @cols) = @_;
	my @out;
	my %cols_seen;
	for my $e (@cols) { $cols_seen{$e}++; }

	# reset header
	my @hd_old = matrix_header_get($mat_name);
	$MATRIX_HEADER{$mat_name} = undef;
	my @new_hd;
	my $i = 0;
	for my $idx (@cols) {
		$MATRIX_HEADER{$mat_name}{$i++} = $hd_old[$idx];
	}

	# assign new data
	my @new;
	for (my $i=0; $i <@$mat_ref; $i++) {
		my @new_row;
		for my $idx (@cols) {
			push @new_row, $$mat_ref[$i][$idx];
		}
		push @new, \@new_row;
	}

	return @new;
}


sub find_in_matrix {
	my ($mode, $mat_ref, $mat_name, $col_name, $find_value) = @_;
	my @out;
	my $idx = matrix_header_get_index($mat_name, $col_name);
	my $cnt_in1 = @$mat_ref;

	for (my $i = 0; $i < @$mat_ref; $i++) {
		if ($mode eq "in") {
			if ($$mat_ref[$i][$idx] eq $find_value) { push @out, $i; }
		} elsif ($mode eq "not in") {
			if ($$mat_ref[$i][$idx] ne $find_value) { push @out, $i; }
		}
	}
	my $cnt_out = @out;
	logging("find $mode done: input count($cnt_in1),  found count: $cnt_out")  if $VERBOSE;

	return @out;
}


sub matrix_header_get {
	my ($matrix_name) = @_;
	my $ref_hash = \%MATRIX_HEADER;
	my $ref = $$ref_hash{$matrix_name};
	my @out;
	for my $i ( sort { $a <=> $b} keys %$ref ) {
		push @out, $$ref{$i};
	}
	return @out;
}


sub matrix_header_create {
	my ($matrix_name, $matrix_first_row_array_ref, @names) = @_;

	my @matrix_first_row_cols;
	if(defined $matrix_first_row_array_ref) {
		@matrix_first_row_cols = @$matrix_first_row_array_ref;
	}

	my $ref_hash = \%MATRIX_HEADER;
	$MATRIX_HEADER{$matrix_name} = undef; # clear;

	my $col_v = 1;
	if (! @names || @names < 0) {
		for (my $i=0; $i<@matrix_first_row_cols; $i++) {
			push @names, "V$i";
			$$ref_hash{$matrix_name}{$i} = "V$col_v";
			$col_v++;
		}
	} else {
		for (my $i=0; $i<@names; $i++) {
			$$ref_hash{$matrix_name}{$i} = $names[$i];
			$col_v++;
		}
	}

	$col_v--;
	logging("Matrix($matrix_name) header created: total $col_v columns.")  if $VERBOSE;
}


sub matrix_header_insert {
	my ($matrix_name, $name, $idx) = @_;

	if (! exists $MATRIX_HEADER{$matrix_name}) {
		$MATRIX_HEADER{$matrix_name} = {};
	}

	my $ref = $MATRIX_HEADER{$matrix_name};

	if (!defined $idx || $idx eq "") {
		my @tmp = keys %$ref;
		$idx = @tmp;
	}

	$$ref{$idx} = $name;
	logging("Matrix($matrix_name) index added: index($idx), name($name)") if $VERBOSE;
	return $idx;
}


sub matrix_header_remove {
	my ($matrix_name, $idx) = @_;
	my $ref = $MATRIX_HEADER{$matrix_name};
my $cnt1 = scalar keys %$ref;

	for my $i (keys %$ref) {
		if ($i == $idx) {
			delete $$ref{$i};
			last;
		}
	}
my $cnt2 = scalar keys %$ref;
	logging("Matrix($matrix_name) index removed: index($idx)") if $VERBOSE;
}


sub matrix_header_rename {
	my ($matrix_name, $name_old, $name_new) = @_;
	my $ref = $MATRIX_HEADER{$matrix_name};

	for my $i (keys %$ref) {
		if ( $$ref{$i} eq $name_old) {
			$$ref{$i} = $name_new;
			last;
		}
	}
	logging("Matrix($matrix_name) header renamed: old($name_old) new($name_new)") if $VERBOSE;
}


sub array_to_matrix {
	my (@arr) = @_;
	my @matrix;

	for (my $i=0; $i < @arr; $i++) {
		my $line = $arr[$i];
		chomp $line;
		next if $line =~ /^\s*$/;

		my @cols = (split "\t", $line);
		push @matrix, \@cols;
	}
	return @matrix;
}


sub matrix_header_get_index {
	my ($mat_name, $header_name) = @_;

	my $ref = $MATRIX_HEADER{$mat_name};
	my @idx = sort { $a <=> $b } keys %$ref;

	if (@idx) {
		for (my $i=0; $i < @idx; $i++) {
			if (! defined $$ref{ $i }) { next; }
			if ( $$ref{ $i } eq $header_name) {
				return $i;
			}
		}
	}

	return -1;
}


sub matrix_header_get_name {
	my ($mat_name, $header_idx) = @_;

	my $ref = $MATRIX_HEADER{$mat_name};
	my @idx = sort { $a <=> $b } keys %$ref;

	if (@idx) {
		for (my $i=0; $i < @idx; $i++) {
			if (! defined $$ref{ $i }) { next; }
			if ( $i eq $header_idx) {
	#			logging("'$header_name' index is $i");
				return $$ref{ $i };
			}
		}
	}

	return -1;
}


sub matrix_header_get_ncol {
	my ($mat_name) = @_;

	my $ref = $MATRIX_HEADER{$mat_name};
	my @idx = sort { $a <=> $b } keys %$ref;

	if (@idx) {
		return scalar @idx;
	}

	return -1;
}


sub matrix_put_values {
	my ($ref_matrix, $matrix_name, $header_name, @values) = @_;

	my $idx_to_put = matrix_header_get_index($matrix_name, $header_name);
	if ($idx_to_put == -1) {
		$idx_to_put = matrix_header_insert($matrix_name, $header_name);
	}
	my $m_row = @$ref_matrix;

	if ($m_row == 0) {
		$m_row = @values;
	}

	my $last_value = $values[ $#values ];
	my $cnt_values = @values;
	if ($cnt_values < $m_row) {
		my $cnt_to_fill = $m_row - $cnt_values;

		my $to_fill_value = $values[$#values];
		my $start = $#values + 1;
		my $end = $#values + 1 + $cnt_to_fill;

		for (my $i=$start; $i < $end; $i++) {
			push @values, $last_value;
		}
	}

	for (my $i=0; $i < $m_row; $i++) {
		my $row_ref = $$ref_matrix[$i];
		$$ref_matrix[$i][$idx_to_put] = $values[$i];
	}
}


sub matrix_sum_cols {
	my ($mat_ref, $idx_col) = @_;

	my $sum = 0;
	for (my $i=0; $i < @$mat_ref; $i++) {
		$sum += $$mat_ref[$i][$idx_col];
	}

	return $sum;
}


sub matrix_stat_col {
	my ($mode, $mat_ref, $idx_col) = @_;

	my $sum = 0;
	my $min;
	my $max;
	my $mean = 0;

	for (my $i=0; $i < @$mat_ref; $i++) {
		$sum += $$mat_ref[$i][$idx_col];

		if(!defined $min) { $min = $$mat_ref[$i][$idx_col]; }
		else { if ($$mat_ref[$i][$idx_col] < $min) { $min = $$mat_ref[$i][$idx_col]; } }

		if(!defined $max) { $max = $$mat_ref[$i][$idx_col]; }
		else { if ($$mat_ref[$i][$idx_col] > $max) { $max = $$mat_ref[$i][$idx_col]; } }

	}

	$mean = $sum / @$mat_ref;

	if ($mode eq "min") { return $min; }
	if ($mode eq "max") { return $max; }
	if ($mode eq "sum") { return $sum; }
	if ($mode eq "mean") { return $mean; }
}


sub matrix_calc_best_quality {
	my ($mat_ref, $idx_col) = @_;

	my $sum = 0;
	my $max;
	my $max_idx;

	for (my $i=0; $i < @$mat_ref; $i++) {
		my $sum = sum_value( split ";", $$mat_ref[$i][$idx_col] );

		if(!defined $max) {
			$max = $sum;
			$max_idx=$i;
		}
		else {
			if ($sum > $max) {
				$max = $sum;
				$max_idx=$i;
			}
		}
	}
	return $max_idx;
}


sub matrix_calc_priorP {
	my ($mat_ref, $header_name) = @_;

	my $idx_priorP = matrix_header_get_index("dat_n", "priorP");
	my $idx_readN = matrix_header_get_index("dat_n", "readN");
	my $sum_readN = matrix_sum_cols($mat_ref, $idx_readN);

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($sum_readN) {
			my $value = $$mat_ref[$i][$idx_readN] / $sum_readN;
			$$mat_ref[$i][$idx_priorP] = $value;
		} else {
			$$mat_ref[$i][$idx_priorP] = 0;
		}
	}
}


sub matrix_calc_p_edit {
	my ($mat_ref, $header_name) = @_;

	my $idx_p_edit = matrix_header_get_index("dat_n", "p_edit");
	my $idx_qv = matrix_header_get_index("dat_n", "qv");

	for (my $i=0; $i < @$mat_ref; $i++) {
		my $value = 10**(-($$mat_ref[$i][$idx_qv] / 10));
		$$mat_ref[$i][$idx_p_edit] = $value;
	}
}


sub matrix_calc_likelihood {
	my ($mat_ref, $header_name) = @_;

	my $idx_likelihood = matrix_header_get_index("dat_n", "likelihood");
	my $idx_priorP = matrix_header_get_index("dat_n", "priorP");
	my $idx_p_edit = matrix_header_get_index("dat_n", "p_edit");

	for (my $i=0; $i < @$mat_ref; $i++) {
		my $value = $$mat_ref[$i][$idx_priorP] * $$mat_ref[$i][$idx_p_edit];
		$$mat_ref[$i][$idx_likelihood] = $value;
	}
}


sub matrix_calc_postP {
	my ($mat_ref, $header_name) = @_;

	my $idx_postP = matrix_header_get_index("dat_n", "postP");
	my $idx_likelihood = matrix_header_get_index("dat_n", "likelihood");
	my $sum_likelihood = matrix_sum_cols($mat_ref, $idx_likelihood);

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($sum_likelihood) {
			my $value = $$mat_ref[$i][$idx_likelihood] / $sum_likelihood;
			$$mat_ref[$i][$idx_postP] = $value;
		} else {
			$$mat_ref[$i][$idx_postP] = 0;
		}
	}
}


sub matrix_calc_mtx_extract {
	my ($mat_ref) = @_;

	my $idx_replaced_BC = matrix_header_get_index("MTX_changed", "replaced_BC");
	my @out;
	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_replaced_BC] ne "NO_in_SAVED" &&
			$$mat_ref[$i][$idx_replaced_BC] ne "Not_Sig" &&
			$$mat_ref[$i][$idx_replaced_BC] ne "No_Hamm_1") {

			push @out, $$mat_ref[$i];
		}
	}
	return @out;
}


sub matrix_calc_MTX_changed_extracted_put_bc {
	my ($mat_ref) = @_;

	my $idx_original_BC = matrix_header_get_index("MTX_changed_extracted", "original_BC");
	my $idx_replaced_BC = matrix_header_get_index("MTX_changed_extracted", "replaced_BC");
	for (my $i=0; $i < @$mat_ref; $i++) {
		$$mat_ref[$i][$idx_original_BC] = $$mat_ref[$i][$idx_replaced_BC];
	}
}


sub matrix_save {
	my ($mtx_ref, $file, $mat_name) = @_;
	my @mtx = @$mtx_ref;

	open F, ">$file" or die $!;

	my $cnt = 0;
	my $col_cnt = @{$mtx[0]};
	my %cols_stat;

	if (defined $mat_name) {
		my $ref = $MATRIX_HEADER{$mat_name};
		my @idx = sort { $a <=> $b } keys %$ref;
		my @vals = @$ref{@idx};

		my $col_cnt = @idx;
		$cols_stat{$col_cnt}++;

		my $row_cnt = @mtx;
		print F join "\t", @vals;
		print F "\n";
	}

	no warnings;
	for my $e (@mtx) {
		my $tmp = join "\t", @$e;
		print F $tmp."\n";
		$cnt++;
	}
	use warnings;

	logging("Matrix file saved: $file ($cnt rows x $col_cnt cols)");
	if (keys %cols_stat >= 2) {
		for my $k (sort { $a <=> $b} keys %cols_stat) {
			logging("\t$k columns: $cols_stat{$k}");
		}
	}
}


sub stringdist_hamm {
	my ($str, $mm_max_cnt) = @_;

	my @out;
	my $mm_cnt = 0;

	for my $e (@WhitelistBarcodes) {
		$mm_cnt = ($str ^ $e) =~ tr/\001-\255//; # XOR
		push @out, $e if $mm_cnt  <= $mm_max_cnt;
	}

	@out;
}


sub str_diff_pos_1 {
	my @s1 = split //, $_[0];
	my @s2 = split //, $_[1];

	for (my $i=0; $i < @s1; $i++) {
		if ($s1[$i] ne $s2[$i]) { return $i; }
	}
}


sub arr_max {
	my @tmp = sort { $b <=> $a } @_;
	return $tmp[0];
}


sub min_value {
	return (sort { $a <=> $b } @_)[0];
}


sub sum_value {
	my (@data) = @_;
   my $sum=0;

   for my $e (@data) {
	   $e =~ s/^\s*//;
	   $e =~ s/\s*$//;
	   if (defined $e) {
		   if ($e ne "") {
				$sum += $e;
		   }
	   }
   }
   return $sum;
}


sub matrix_calc_maxdat_n {
	my ($mat_ref, $mat_name, $col_name) = @_;
	my $idx_postP = matrix_header_get_index($mat_name, $col_name);
	my $max = arr_max( matrix_get_array_by_col($mat_ref, $idx_postP ) );

	my @out;
	for my $e_ref (@$mat_ref) {

		if (${$e_ref}[$idx_postP] == $max) {
			push @out, $e_ref;
		}
	}
	@out;
}


sub matrix_remove_col {
	my ($mat_ref, $mat_name, $col_idx) = @_;

	matrix_header_remove($mat_name, $col_idx);

	my @out;
	for my $e_ref (@$mat_ref) {
		my @tmp;
		for (my $i=0; $i < @$e_ref; $i++) {
			if ($i != $col_idx) {
				push @tmp, ${$e_ref}[$i];
			}
		}
		push @out, \@tmp;
	}
	@out;
}


sub matrix_calc_MTX_RE_BARCODE_UMI {
	my ($mat_ref, $mat_name) = @_;

	my $idx_v1 = matrix_header_get_index($mat_name, "V1");
	my $idx_v3 = matrix_header_get_index($mat_name, "V3");
	my $idx_barcode_umi = matrix_header_get_index($mat_name, "BARCODE_UMI");

	for (my $i=0; $i < @$mat_ref; $i++) {
		my $value = $$mat_ref[$i][$idx_v1]."_".$$mat_ref[$i][$idx_v3];
		$$mat_ref[$i][$idx_barcode_umi] = $value;
	}
}


sub matrix_calc_MTX_RE_d_mismatch_WT_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_wt = matrix_header_get_index($header_name, "mismatch.WT");
	my $idx_wt_call = matrix_header_get_index($header_name, "mismatch.WT.call");

	my @out;
	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_wt] <= $WT_ALW) {
			push @out, $$mat_ref[$i];
		}
	}

	if (@out > 0) {
		for (my $i=0; $i < @$mat_ref; $i++) {
			if ($$mat_ref[$i][$idx_wt] <= $WT_ALW) {
				$$mat_ref[$i][$idx_wt_call] = "WT";
			}
		}
	}
}


sub matrix_calc_MTX_MERGE_u_WT_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_wt = matrix_header_get_index($header_name, "mismatch.WT");
	my $idx_wt_call = matrix_header_get_index($header_name, "mismatch.WT.call");

	my @out;
	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_wt] <= $WT_ALW) {
			push @out, $$mat_ref[$i];
		}
	}

	if (@out > 0) {
		for (my $i=0; $i < @$mat_ref; $i++) {
			if ($$mat_ref[$i][$idx_wt] <= $WT_ALW) {
				$$mat_ref[$i][$idx_wt_call] = "WT";
			}
		}
	}
}


sub matrix_calc_MTX_MERGE_u_MUT_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_mut = matrix_header_get_index($header_name, "mismatch.MUT");
	my $idx_mut_call = matrix_header_get_index($header_name, "mismatch.MUT.call");

	my @out;
	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_mut] <= $WT_ALW) {
			push @out, $$mat_ref[$i];
		}
	}

	if (@out > 0) {
		for (my $i=0; $i < @$mat_ref; $i++) {
			if ($$mat_ref[$i][$idx_mut] <= $WT_ALW) {
				$$mat_ref[$i][$idx_mut_call] = "MUT";
			}
		}
	}
}


sub matrix_calc_MTX_MERGE_u_WT_call_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_1 = matrix_header_get_index($header_name, "mismatch.WT.call");
	my $idx_2 = matrix_header_get_index($header_name, "mismatch.call");

	my @out;
	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_1] eq "WT") {
			push @out, $$mat_ref[$i];
		}
	}

	if (@out > 0) {
		for (my $i=0; $i < @$mat_ref; $i++) {
			if ($$mat_ref[$i][$idx_1] eq "WT") {
				$$mat_ref[$i][$idx_2] = "WT";
			}
		}
	}
}


sub matrix_calc_MTX_MERGE_u_MUT_call_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_1 = matrix_header_get_index($header_name, "mismatch.MUT.call");
	my $idx_2 = matrix_header_get_index($header_name, "mismatch.call");

	my @out;
	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_1] eq "MUT") {
			push @out, $$mat_ref[$i];
		}
	}

	if (@out > 0) {
		for (my $i=0; $i < @$mat_ref; $i++) {
			if ($$mat_ref[$i][$idx_1] eq "MUT") {
				$$mat_ref[$i][$idx_2] = "MUT";
			}
		}
	}
}


sub matrix_calc_MTX_MERGE_u_amb_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_1 = matrix_header_get_index($header_name, "mismatch.WT.call");
	my $idx_2 = matrix_header_get_index($header_name, "mismatch.MUT.call");
	my $idx_3 = matrix_header_get_index($header_name, "mismatch.call");

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_1] eq "" && $$mat_ref[$i][$idx_2] eq "") {
			$$mat_ref[$i][$idx_3] = "amb";
		}
	}
}


sub matrix_calc_MTX_RE_d_mismatch_MUT_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_mut = matrix_header_get_index($header_name, "mismatch.MUT");
	my $idx_mut_call = matrix_header_get_index($header_name, "mismatch.MUT.call");

	my @out;
	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_mut] <= $WT_ALW) {
			push @out, $$mat_ref[$i];
		}
	}

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_mut] <= $WT_ALW) {
			$$mat_ref[$i][$idx_mut_call] = "MUT";
		}
	}
}


sub matrix_calc_MTX_RE_d_mismatch_WT_call_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_wt_call = matrix_header_get_index($header_name, "mismatch.WT.call");
	my $idx_mm_call = matrix_header_get_index($header_name, "mismatch.call");

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_wt_call] eq "WT") {
			$$mat_ref[$i][$idx_mm_call] = "WT";
		}
	}
}


sub matrix_calc_MTX_RE_d_mismatch_MUT_call_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_mut_call = matrix_header_get_index($header_name, "mismatch.MUT.call");
	my $idx_mm_call = matrix_header_get_index($header_name, "mismatch.call");

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_mut_call] eq "MUT") {
			$$mat_ref[$i][$idx_mm_call] = "MUT";
		}
	}
}


sub matrix_calc_MTX_RE_d_mismatch_amb_assign {
	my ($mat_ref, $header_name) = @_;

	my $idx_wt_call = matrix_header_get_index($header_name, "mismatch.WT.call");
	my $idx_mm_mut_call = matrix_header_get_index($header_name, "mismatch.MUT.call");
	my $idx_mm_call = matrix_header_get_index($header_name, "mismatch.call");

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ($$mat_ref[$i][$idx_mm_mut_call] eq "MUT" && $$mat_ref[$i][$idx_wt_call] eq "WT" ) {
			$$mat_ref[$i][$idx_mm_call] = "amb";
		}
	}
}


sub matrix_calc_get_mismatch_call {
	my ($mat_ref, $header_name) = @_;

	my $idx_1 = matrix_header_get_index($header_name, "mismatch.call");
	my $cnt_wt=0;
	my $cnt_mut=0;
	my $cnt_amb=0;

	for (my $i=0; $i < @$mat_ref; $i++) {
		if ( $$mat_ref[$i][$idx_1] eq "WT") { $cnt_wt++; }
		if ( $$mat_ref[$i][$idx_1] eq "MUT") { $cnt_mut++; }
		if ( $$mat_ref[$i][$idx_1] eq "amb") { $cnt_amb++; }
	}

	return ($cnt_wt, $cnt_mut, $cnt_amb);

}


sub sed_repeat {
	my ($file_in, $file_out, $line_no_select, $lines_repeat,  $cut_start, $cut_end) = @_;

	$cut_start-- if defined $cut_end;; # Convert into zero-based index
	$cut_end-- if defined $cut_end;
	my @out;

	open FILE, $file_in or die "[$file_in] $!";

	my $line_cnt = 1;
	for my $line (<FILE>) {
		chomp $line;

		if ($line_cnt == $line_no_select) {
			if (defined $cut_start) {
				my $selected = join "", (split //, $line)[$cut_start .. $cut_end];
				push @out, $selected;
			} else {
				push @out, $line;
			}
		}

		my $mod = $line_cnt % $lines_repeat;
		if ( $mod == 0) {
			$line_cnt = 1;
		} else {
			$line_cnt++;
		}
	}
	close FILE;

	if (defined $file_out) {
		open FILE, ">$file_out" or die $!;
		print FILE join "\n", @out;
		close FILE;
	} else {
		return @out;
	}
}


sub paste_cols {
	my @files = @_[ 0 .. $#_ -1];
	my $out_file = $_[ $#_ ];

	my @paste_cols_out = ();
	my $idx = 0;

	for my $file (@files) {
		matrix_header_insert("paste_cols_out", $idx, $idx);
		my $ref = file_read($file);
		my @data = @$ref;
		for (@data) { chomp; }

		matrix_put_values(\@paste_cols_out, "paste_cols_out", "$idx", @data);
		$idx++;
	}

	matrix_save(\@paste_cols_out, $out_file);
}


sub load_phred_scores {

	my $table = <<'__end_of_string__';
Q-Score	P_error	ASCII_Code	Symbol
0	1.00000	33	!
1	0.79433	34	"
2	0.63096	35	#
3	0.50119	36	$
4	0.39811	37	%
5	0.31623	38	&
6	0.25119	39	'
7	0.19953	40	(
8	0.15849	41	)
9	0.12589	42	*
10	0.10000	43	+
11	0.07943	44	,
12	0.06310	45	-
13	0.05012	46	.
14	0.03981	47	/
15	0.03162	48	0
16	0.02512	49	1
17	0.01995	50	2
18	0.01585	51	3
19	0.01259	52	4
20	0.01000	53	5
21	0.00794	54	6
22	0.00631	55	7
23	0.00501	56	8
24	0.00398	57	9
25	0.00316	58	:
26	0.00251	59	;
27	0.00200	60	<
28	0.00158	61	=
29	0.00126	62	>
30	0.00100	63	?
31	0.00079	64	@
32	0.00063	65	A
33	0.00050	66	B
34	0.00040	67	C
35	0.00032	68	D
36	0.00025	69	E
37	0.00020	70	F
38	0.00016	71	G
39	0.00013	72	H
40	0.00010	73	I
41	0.00008	74	J
__end_of_string__

	$table =~ s/^\s*//;
	$table =~ s/\s*$//;
	my @arr = split /\r?\n/, $table;

	my %out;
	for my $line (@arr) {
		chomp $line;
		my ($q_score, $p_err, $ascii, $symbol) = split "\t", $line;
		$out{$symbol}{"q_score"} = $q_score;
		$out{$symbol}{"p_err"} = $p_err;
#		print "$symbol: $q_score / $p_err\n";
	}
	return %out;
}

sub un_gzip {
	my ($f1, $f2) = @_;
	my $f_out;

	if ($f1 =~ /\.gz$/i) {
		$f_out = $f1;
		$f_out =~ s/\.gz$//i;

		unless(gunzip($f1, $f_out)) {
			logging("Error decompressing '$f1': $GunzipError");
			exit;
		}
		$f1 = $f_out;
	}

	if ($f2 =~ /\.gz$/i) {
		$f_out = $f2;
		$f_out =~ s/\.gz$//i;

		unless(gunzip($f2, $f_out)) {
			logging("Error decompressing '$f2': $GunzipError");
			exit;
		}
		$f2 = $f_out;
	}

	return ($f1, $f2);

}


sub stringdist_hamm_pre_calc {
	my $MTX_sub_ref = shift;

	my $batch_amount;

	$batch_amount = int(@$MTX_sub_ref / $thread_max);
	my $left = @$MTX_sub_ref % ($batch_amount * $thread_max);
	$batch_amount += $left;

	my $array_index = 0;

	for my $thread_cnt (1 .. $thread_max) {

		my $array_index_to = $array_index + $batch_amount - 1;
		if ($array_index_to >= @$MTX_sub_ref - 1) { $array_index_to = (@$MTX_sub_ref - 1); }

		my @data_batch = @$MTX_sub_ref[ $array_index .. $array_index_to ];
		my $len = @data_batch;
		last unless @data_batch;

		logging("Thread $thread_cnt: data len($len)  from($array_index)  to($array_index_to)");
		$array_index = $array_index_to + 1;

		my ($t) = threads->new( \&stringdist_hamm_pre_calc_run, @data_batch );
	}

	#####################################
	# Wait for all thread runs to complete
	#####################################
	for my $t (threads->list()) {
		my %hash = $t->join;
		for my $k (keys %hash) {
			$BARCODE_stringdist_hamm{$k} = $hash{$k};
		}
	}
}

sub stringdist_hamm_pre_calc_run {
	my @data = @_;
	my %hash;

	for (my $r=0; $r < @data; $r++) {
		my $ori_BC = $data[$r][0];
		my @Hamm_1 = stringdist_hamm($ori_BC, 1);
		$hash{$ori_BC} = \@Hamm_1;
	}

	return %hash;
}

__END__
