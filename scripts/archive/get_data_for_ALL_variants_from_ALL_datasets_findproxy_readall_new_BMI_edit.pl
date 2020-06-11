#!/usr/bin/perl -w
use DBI;

## Define variables
$PVCUTOFF =    5e-8;		# P-Value cutoff
$LD_SCORE =     0.0;		# LD score
$DISTANCE =  000000;		# Distance

$NUM_TRAIT = 44; 				# number of traits

$SOURCE_PATH = "/humgen/diabetes2/users/mvg/VariantClustering";

## Argument variables
if(@ARGV > 1) { 	# if there are more than 1 argument
	$PH  = "(PH = '$ARGV[1]')";			# take the second argument as phenotype info
	$VER = "";
}
else { 					# if @ARGV <= 1; if there is 1 argument or less
								# use different versions for each phenotype
#	$PH  = "(PH = 'BIP' || PH = 'MDD' || PH = 'SCZ')";
#	$VER = "and (VER = 'mdv37')";

#	$PH  = "(PH = 'allstroke' || PH = 'IS' || PH = 'Stroke_IS')";
#	$VER = "and (VER = 'mdv73')";

	$PH  = "(PH = 'BMI')";
	$VER = "and (VER = 'mdv37')";

#	$PH  = "(PH = 'T2D')";
#	$VER = "and (VER = 'mdv36')";
#	$VER = "and (VER = 'mdv39')";			# T2D related phenotypes - newer version

#	$PH  = "(PH = 'CAD' || PH = 'MI')";
#	$VER = "and (VER = 'mdv96')";

}


#QUERY 1		# join datasets
$query1 = "select distinct META_DATASET_PH.TBL from META_MDV join META_DATASET on META_MDV.Dataset = META_DATASET.Dataset join META_DATASET_PH on META_DATASET.Dataset = META_DATASET_PH.Dataset where PARENT = 'root' and $PH order by SUBJECTS desc;";

#die $query1;

#get Ref DB names			 ?????????????????
$dbh = DBI->connect("DBI:mysql:database=digkb;host=db-dev-seventeen-aws.cxrzznxifeib.us-east-1.rds.amazonaws.com", "diguser", "type2diabetes", {'RaiseError' => 1});
$sth = $dbh->prepare($query1);
$sth->execute();
$arr=$sth->fetchall_arrayref;
foreach(@{$arr}) { foreach(@{$_}) { push @refDBs, $_ } }	# get each of the Ref DB names
die "No refDBs provided\n" if @refDBs == 0;								# die if there is no Ref DB

#read p-value, odds-ratio, and beta column names and sort order --> # get row names and sort number					?????????????
$propFileName = "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/meta_data_file-PROP.txt";
foreach(split /\n/, `awk '\$9=="P_VALUE"'    $propFileName | cut -f1,8`) { /^([^\t]+)\t(\d+)$/ || die; $pval2sort{$1} = $2 }
foreach(split /\n/, `awk '\$9=="ODDS_RATIO"' $propFileName | cut -f1,8`) { /^([^\t]+)\t(\d+)$/ || die; $odds2sort{$1} = $2 }
foreach(split /\n/, `awk '\$9=="BETA"'       $propFileName | cut -f1,8`) { /^([^\t]+)\t(\d+)$/ || die; $beta2sort{$1} = $2 }

# to read variant list from an input file
if(@ARGV>0) {		# if there's an input argument
	$takeVarinatsFromFile = 1;
	open(IN,$ARGV[0]) || die "Cannot open $ARGV[0]\n";
	$HEADER = <IN>;
	$HEADER =~ s/[\r\n]//g;		# split the header by line
	@d = split /\t/, $HEADER;	# split the header by tab
	if(@d > 3 and $d[3] eq "p" and $d[4] eq "dir") {		# if 3rd column is pvalue & 4th column is direction
		$takeDirectionsFromFile = 1;		# then we'll take directions from the input file
	}
	else {
		$takeDirectionsFromFile = 0;
	}
	if($d[1] eq "locus"){
		$takeGenenameFromFile = 1;
	}
	while(<IN>) {		# while there's input from the file
		chomp;	 # remove (usually) any newline character from the end of a string.
		@d = split /\t/, $_;	# split by tab
		$idsToTake{$d[2]} = 1;
		if($takeDirectionsFromFile) {	# when we have direction info in the file
			$id2pv{    $d[2]} = $d[3];	# key = IDs, value = p-value
			$dir{      $d[2]} = $d[4];	# key = IDs, value = direction
			$takeDirectionsFromFile = 1;
		}
		if($takeGenenameFromFile){
			$id2gene{$d[2]} = $d[1];
		}

	}
}
else {	# we don't have an input file to read variant list or directions
	$takeVarinatsFromFile   = 0;
	$takeDirectionsFromFile = 0;
}

#get list of SNPs and directions from PHENTOYPE files
foreach $refDB (@refDBs) {
	last if $takeDirectionsFromFile == 1;	# break if we take directions from file
	($ds,$ph) = split /__/, $refDB;		# dataset/phenotype
	for($i = 1; ; $i++) {
		#take the next file and check if exists
		$fileName = "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/$ds/$ph/DATA/$ds.$ph.$i.txt";
		last if !-e $fileName;

		#read header and find effect-size and p-value column names
		open(IN,$fileName) || die "Cannot open $fileName\n";
		chomp($header = <IN>);
		# initialize the variables to sort pvalue, odds ratio, and beta
		$pvalIndx = -1; $pvalSort = 99999;
		$effcIndx = -1; $effcSort = 99999; $effcNull = -1;

		$i = 0;
		foreach(split /\t/, $header) {		# sort the pvalues, odds ratios, and betas			 # ????????????????????????
			if(exists $pval2sort{$_} && $pval2sort{$_} < $pvalSort) { $pvalIndx = $i ; $pvalSort = $pval2sort{$_} }
			if(exists $odds2sort{$_} && $odds2sort{$_} < $effcSort) { $effcIndx = $i ; $effcSort = $odds2sort{$_} ; $effcNull = 1 }		# OR = 1 for null
			if(exists $beta2sort{$_} && $beta2sort{$_} < $effcSort) { $effcIndx = $i ; $effcSort = $beta2sort{$_} ; $effcNull = 0 } 	# beta = 0 for null
			$i++;
		}
		next if $pvalIndx < 0 || $effcIndx < 0 || $effcNull < 0;

		print STDERR "Processing source file: $fileName\n";
		print "$ds.$ph\n";
		#add reference to te output
		#($tmp=$ds)=~s/^(GWAS_|ExSeq_|ExChip_|WGS_)(.+)/$2/;
		#$takeMe{"$tmp\-$ph\-Z"} = 1;

		if($takeVarinatsFromFile) {
			open(IN,$fileName) || die "Cannot open $fileName\n";
			while(<IN>) {
				/^([^\t]+)/;
                     		next if !exists $idsToTake{$1};			# the first one from the last regular expression -- if the value exists, we consider the var id.
				s/[\r\n]//g;
				@d = split /\t/, $_;		# split by tab
				next if $d[0] !~ /^(\d+)\_(\d+)\_([ACTG])\_([ACTG])$/;  # decimal, decimal, base, base (alleles)
				next if exists $id2pv{$d[0]} && ($d[$pvalIndx] eq "NA" || $id2pv{$d[0]} < $d[$pvalIndx]);			# ??????????????
				$id2pv{$d[0]} = $d[$pvalIndx];
				$dir{  $d[0]} = $d[$effcIndx] > $effcNull ? 1 : -1;			# if the effect is greater than null effect, then direction = 1, otherwise -1.
			}

		}
		else {		# when there is no argument variable for input variant file
			open(IN,"awk '\$".($pvalIndx+1)." <= $PVCUTOFF' $fileName |") || die "Cannot awk $fileName\n";		# ???????????
			while(<IN>) {
				s/[\r\n]//g;
				@d = split /\t/, $_;
				next if $d[0] !~ /^(\d+)\_(\d+)\_([ACTG])\_([ACTG])$/;
				next if exists $id2pv{$d[0]} && $id2pv{$d[0]} < $d[$pvalIndx];
				$id2pv{$d[0]} = $d[$pvalIndx];
				$dir{  $d[0]} = $d[$effcIndx] > $effcNull ? 1 : -1;
			}
		}
	}
}



######################################################### new code starts #################################################

print STDERR "Converting var ID to rsID...\n";

open(IN,"/humgen/diabetes2/users/clairekim/list_VARID_rsID.txt") || die "Cannot read from var ID data\n";
while(<IN>){
	s/[\r\n]//g;
	@d = split /\t/, $_;
	$rs2id{$d[1]} = $d[0];
	$id2rs{$d[0]} = $d[1];
}

#exit;
print STDERR "Reading counts...\n";

open(IN,"/humgen/diabetes2/users/clairekim/VARID_to_counts_44traits.txt") || die "Cannot read from counts file\n";
while(<IN>) {
	s/[\r\n]//g;
	($id,$count)=split /\t/, $_;		# ???????????? what are the two numbers ?
#	next if !exists $id2pv{$id};
	$id2cnt{$id} = $count;
}


print STDERR "Reading counts for variants that are not in current VARID to counts data...\n";

$FILTER = 1;
if($FILTER) {
	#open(IN,"cut -f2 $SOURCE_PATH/miriam-marcin-pheno-match-after-adding-FA.txt |") || die;		# dataset name
#	open(IN,"cut -f2 /humgen/diabetes2/users/clairekim/miriam-marcin-pheno-match-after-adding-FA.txt |") || die;
	open(IN,"cut -f2 /humgen/diabetes2/users/clairekim/miriam-marcin-pheno-match-after-adding-FA-edit-pval-0.01.txt |") || die;
	while(<IN>) { chomp; $takeMe{$_} = 1 }	# $takeMe ??????????
}

#QUERY 2						???????????????????????
$query2 = "select TBL from META_MDV join META_DATASET_PH on META_MDV.Dataset = META_DATASET_PH.Dataset join META_PROP_DATASET_PH on (META_DATASET_PH.Dataset = META_PROP_DATASET_PH.Dataset and META_DATASET_PH.PH = META_PROP_DATASET_PH.PH) where Parent = 'Root' and Ver = 'mdv33' and PROP = 'BETA' and (TBL like '%CKDGenConsortium%' or TBL like '%GIANT_dv%' or TBL like '%GLGC%' or TBL like '%MAGIC%');";

$query2 = "select distinct META_DATASET_PH.TBL from META_MDV join META_DATASET on META_MDV.DATASET=META_DATASET.DATASET join META_DATASET_PH on META_DATASET.DATASET = META_DATASET_PH.DATASET where TECH='GWAS' and SG not like '%metabochip%' and (VER like 'mdv2%' or VER like 'mdv3%') order by PH, META_DATASET_PH.TBL;";


foreach $id (keys %id2pv){
	$id =~ /^(\d+)\_(\d+)\_([ACTG])\_([ACTG])$/;
	$id2chr{$id} = $1;    	# first chunk --> chromosome
	$id2pos{$id} = $2;			# second chunk --> position
	if(($3 eq "A" and $4 eq "T") ||			# check if it is ambiguous pair -- AT or GC
	   ($3 eq "T" and $4 eq "A") ||
	   ($3 eq "C" and $4 eq "G") ||
	   ($3 eq "G" and $4 eq "C") ){
	   	$id2isATorCG{$id} = 1 }
	else {  $id2isATorCG{$id} = 0 }
}

$sth = $dbh->prepare($query2);
$sth->execute();
$arr=$sth->fetchall_arrayref;
foreach(@{$arr}) {
	foreach(@{$_}) {
		($ds,$ph) = split /__/, $_;
		$fileName = "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/$ds/$ph/DATA/$ds.$ph.1.txt";
		$ds =~ s/^GWAS_//;
		$columnName = "$ds-$ph-Z";

		next if $FILTER && !exists $takeMe{$columnName};

		open(IN,$fileName) || die "Cannot open $fileName\n";
		$header = <IN>;
		$header =~ s/[\r\n]//g;
		$i  =  0;
		$bi = -1;
		$ei = -1;
		$ni = -1;
		$pi = -1;
		foreach(split /\t/, $header) {	# index
			$bi = $i if $_ eq "BETA"; 	# beta
			$ni = $i if $_ eq "N_PH";		# sample size
			$ei = $i if $_ eq "SE";			# standard error
			$pi = $i if $_ eq "P_VALUE";
			$i++;
		}
		next if $bi < 0 || $ei < 0;
		%ps = (); #pvalue
		while(<IN>) {
			/^([^\t]+)/;
			next if !exists $id2pv{$1};
			s/[\r\n]//g;
			@d = split /\t/, $_;
			$id2cnt{$d[0]}++;
		}
	}
}



#################

# read counts from file - convert ID - finding variants that need proxy - get LD block - recount counts for proxies - sort proxies

#################

$allrs = "";
print STDERR "Getting variants that need proxies...\n";
foreach $id (keys %id2pv) {
	if( !exists $id2cnt{$id} || $id2cnt{$id} < $NUM_TRAIT*0.8 ){
		print "$id2rs{$id}\t$id\t$id2cnt{$id}\n";
		$id2proxy{$id} = 1;				# to flag var id that needs a proxy
		if (exists $id2rs{$id}) {
			$allrs .= "$id2rs{$id} ";
		}
	}
	else { $idsToUse{$id} = 1; print "$id2rs{$id}\t$id\t$id2cnt{$id}\n"; }	# to flag var id that do not need a proxy
}



print STDERR "Reading LD_EUR...\n";
if ($allrs){
	open(IN,"$SOURCE_PATH/tabix-0.2.6/tabix $SOURCE_PATH/LD_EUR.tsv.bgz $allrs |") || die "Cannot read from LD_EUR\n";
	while(<IN>) {
		s/[\r\n]//g;
		($rs1,$info)=split /\t/, $_;		# ???????????? what are the two numbers ?
		foreach(split  /;/, $info) {	# each rsid that the current rsid is in LD with
			($rs2,$ld)=split /,/,$_;		# completely same as the code above.
			next if $ld < $LD_SCORE || $ld < 0.8;					# if the LD is less than $LD_SCORE threshold
			$rs2rs{$rs1}{$rs2} = $ld;			#  double hashtable / since variants are in LD in a 1-to-many relationship.
			$rs2rs{$rs2}{$rs1} = $ld;
			next if !exists $rs2id{$rs2};
			$ldscore{$rs2id{$rs1}}{$rs2id{$rs2}} = $ld;
			$id = $rs2id{$rs2};
			$proxyids{$id} = 1;
			$id =~ /^(\d+)\_(\d+)\_([ACTG])\_([ACTG])$/;
			$id2chr{$id} = $1;    	# first chunk --> chromosome
			$id2pos{$id} = $2;			# second chunk --> position
			$id2cnt{$id} = 0;
			if($3){
				if(($3 eq "A" and $4 eq "T") ||			# check if it is ambiguous pair -- AT or GC
					 ($3 eq "T" and $4 eq "A") ||
					 ($3 eq "C" and $4 eq "G") ||
					 ($3 eq "G" and $4 eq "C") ){
						$id2isATorCG{$id} = 1 }
				else {  $id2isATorCG{$id} = 0 }
			}
		#	print "$rs1\t$rs2\t$ld\n";
		}
	}

}
#exit;



########################### added code for filtering traits ########################################
print STDERR "Getting counts for proxies...\n";


if ($proxyids > 0){
	$sth = $dbh->prepare($query2);
	$sth->execute();
	$arr=$sth->fetchall_arrayref;
	foreach(@{$arr}) {
		foreach(@{$_}) {
			($ds,$ph) = split /__/, $_;
			$fileName = "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/$ds/$ph/DATA/$ds.$ph.1.txt";
			$ds =~ s/^GWAS_//;
			$columnName = "$ds-$ph-Z";

			next if $FILTER && !exists $takeMe{$columnName};

			open(IN,$fileName) || die "Cannot open $fileName\n";
			$header = <IN>;
			$header =~ s/[\r\n]//g;
			$i  =  0;
			$bi = -1;
			$ei = -1;
			$ni = -1;
			$pi = -1;
			foreach(split /\t/, $header) {	# index
				$bi = $i if $_ eq "BETA"; 	# beta
				$ni = $i if $_ eq "N_PH";		# sample size
				$ei = $i if $_ eq "SE";			# standard error
				$pi = $i if $_ eq "P_VALUE";
				$i++;
			}
			next if $bi < 0 || $ei < 0;
			%ps = (); #pvalue
			while(<IN>) {
				/^([^\t]+)/;
				next if !exists $proxyids{$1};
				s/[\r\n]//g;
				@d = split /\t/, $_;
				$id2cnt{$d[0]}++;
			}
		}
	}
}



####################################### added code ######################################

print STDERR "Sorting potential proxies...\n";																									# takes forever...
foreach $rs1id (keys %id2proxy) {
	$rs1 = $id2rs{$rs1id};
	# sort by ATCG / count / r2 ($score)
	$ldscore{$rs1id}{$rs1id} = 1;
	if (!exists $id2cnt{$rs1id}) { $id2cnt{$rs1id} = 0; $id2isATorCG{$rs1id} = 1; }
	@proxies = sort { $id2isATorCG{$a} <=> $id2isATorCG{$b} || $id2cnt{$b} <=> $id2cnt{$a} || $ldscore{$rs1id}{$b} <=> $ldscore{$rs1id}{$a} || $a cmp $b } keys %{$ldscore{$rs1id}};
	# select one variant and set $idsToTake $dir $id2pv
	foreach $id (@proxies) {				# we want to use one variant per each LD block
		if ($ldscore{$rs1id}{$id} >= 0.8 && $id2cnt{$id} >= $NUM_TRAIT*0.8 && $id2chr{$rs1id} == $id2chr{$id}){
			$proxy{$rs1id} = $id;
			$dir{$id} = $dir{$rs1id}; 	# set direction
			$id2pv{$id} = $id2pv{$rs1id};
			$idsToUse{$id} = 1;
			last;
		}
		else { next; }
	}

}



# exit;


#################################### added code end ############################################ below edited a bit


foreach $id (keys %idsToUse) {		# see if rsid is in a right format
	$id2rsnum{$id} = $id2rs{$id}=~/^rs(\d+)$/ ? $1 : 999999999;		# ???????????
}
$allrs = "";			# quick sort for AT/CG, count number, pvalue, rsnumber .. count / rsnum / a and b -- to make sure we get the same variant every time we run the code
foreach $id (keys %idsToUse){
	if (!exists $id2cnt{$id}){
		$id2cnt{$id} = 0;
		$id2isATorCG{$id} = 1;
	}
	if (!exists $id2pv{$id}) { $id2pv{$id} = 1; }
}
foreach $id (sort { $id2isATorCG{$a} <=> $id2isATorCG{$b} || $id2cnt{$b} <=> $id2cnt{$a} ||
                     $id2pv{$a} <=> $id2pv{$b} || $id2rsnum{$a} <=> $id2rsnum{$b} || $a cmp $b } keys %idsToUse) { # get keys from %id2rs hashtable first and then sort
	if(0 < $DISTANCE) {
		for($i=0; $i < @maybeIDs; $i++) {			# check if they're in LD block		# ??????????????
			last if $id2chr{$id} == $id2chr{$maybeIDs[$i]} and
		    	    abs($id2pos{$id}  - $id2pos{$maybeIDs[$i]}) < $DISTANCE
		}
		next if $i < @maybeIDs;
	}
	push @maybeIDs, $id;
	$rs2yesLD{$id2rs{$id}} = 0;
	next if !exists $id2rs{$id};
	$allrs .= "$id2rs{$id} ";
}

print STDERR "Reading LD_EUR...\n";
if($allrs){
	open(IN,"$SOURCE_PATH/tabix-0.2.6/tabix $SOURCE_PATH/LD_EUR.tsv.bgz $allrs |") || die "Cannot read from LD_EUR\n";
	while(<IN>) {
		s/[\r\n]//g;
		($rs1,$info)=split /\t/, $_;		# ???????????? what are the two numbers ?
		$rs2yesLD{$rs1}++;
		#foreach $myVariable (split  /;/, $info) {
		#	($rs2,$ld)=split /,/, $myVariable;
		foreach(split  /;/, $info) {	# each rsid that the current rsid is in LD with
			($rs2,$ld)=split /,/,$_;		# completely same as the code above.
			next if $ld < $LD_SCORE;					# if the LD is less than $LD_SCORE threshold
			next if !exists $rs2yesLD{$rs2};	# continue; if $rs2yesLD has no LD info on rs2
			$rs2rs{$rs1}{$rs2} = 1;			#  double hashtable / since variants are in LD in a 1-to-many relationship.
			$rs2rs{$rs2}{$rs1} = 1;
		}
	}
}

foreach $id (@maybeIDs) {				# we want to use one variant per each LD block
	$rs = $id2rs{$id};
	next if $rs2yesLD{$rs} == 0;			# we want to use one variant per each LD block -- continue if current rsid has no LD info
	for($i=0; $i < @takeIDs; $i++) {		# see if IDs in takeIDs in LD with maybeIDs				??????????????????
		last if exists $rs2rs{$rs}{$id2rs{$takeIDs[$i]}};			# ?????????????????
	}
	next if $i < @takeIDs;
	push @takeIDs, $id;
	if (exists $id2gene{$id}){
		$gene2count{$id2gene{$id}}++;
	}
	else {
		$sth = $dbh->prepare("select CLOSEST_GENE from common_dv18 where VAR_ID='$id'");
		$sth->execute();
		$arr=$sth->fetchall_arrayref;
		$gene = $$arr[0][0];
		$gene = "XXXXX" if !defined $gene;
		$id2gene{$id} = $gene;
		$gene2count{$gene}++;		# why do we update the count ? -- we need to name the variants in the same gene differently by indexing them.
	}
}
#foreach $id (@takeIDs){
#	if (!exists $id2chr{$id}) {
#		$id =~ /^(\d+)\_(\d+)\_([ACTG])\_([ACTG])$/;
#		$id2chr{$id} = $1;
#		$id2pos{$id} = $2;
#	}
#}
foreach $id (sort { $id2chr{$a} <=> $id2chr{$b} || $id2pos{$a} <=> $id2pos{$b} } @takeIDs) {
	$ids{$id} = 1;			# to record if the id is present
	$gene = $id2gene{$id};
	if($gene2count{$gene} > 1) {
		$gene2number{$gene}++;
		$gene .= "_$gene2number{$gene}";
	}
	push @ids, $id;
	push @OUT, "$id2rs{$id}\t$gene\t$id\t$id2pv{$id}";
}


#########################################################################################
$HEADER = "SNP\tlocus\tVAR_ID_hg19\tPVALUE";


#read sample size for phenotypes
open(IN,"/humgen/diabetes2/users/clairekim/phenotype-sample-size_edit.txt") || die;
while(<IN>) {
	chomp;
	($pheno,$samplesize) = split /\t/, $_;
	$pheno2samplesize{$pheno} = $samplesize;
}

#QUERY 2						???????????????????????
$query2 = "select TBL from META_MDV join META_DATASET_PH on META_MDV.Dataset = META_DATASET_PH.Dataset join META_PROP_DATASET_PH on (META_DATASET_PH.Dataset = META_PROP_DATASET_PH.Dataset and META_DATASET_PH.PH = META_PROP_DATASET_PH.PH) where Parent = 'Root' and Ver = 'mdv33' and PROP = 'BETA' and (TBL like '%CKDGenConsortium%' or TBL like '%GIANT_dv%' or TBL like '%GLGC%' or TBL like '%MAGIC%');";

$query2 = "select distinct META_DATASET_PH.TBL from META_MDV join META_DATASET on META_MDV.DATASET=META_DATASET.DATASET join META_DATASET_PH on META_DATASET.DATASET = META_DATASET_PH.DATASET where TECH='GWAS' and SG not like '%metabochip%' and (VER like 'mdv2%' or VER like 'mdv3%') order by PH, META_DATASET_PH.TBL;";

#get seffect size for phenotypes
$sth = $dbh->prepare($query2);
$sth->execute();
$arr=$sth->fetchall_arrayref;
foreach(@{$arr}) {
	foreach(@{$_}) {
		($ds,$ph) = split /__/, $_;
		$fileName = "/humgen/diabetes2/users/mvg/portal/scripts/VARIANTS/PHENOTYPES/$ds/$ph/DATA/$ds.$ph.1.txt";
	#	next if !exists $trait{$fileName};
		$ds =~ s/^GWAS_//;
		$columnName = "$ds-$ph-Z";
	#	print "$ds-$ph-Z\n";
		next if $FILTER && !exists $takeMe{$columnName};

		print STDERR "Processing data file: $fileName\n";

		open(IN,$fileName) || die "Cannot open $fileName\n";
		$header = <IN>;
		$header =~ s/[\r\n]//g;
		$i  =  0;
		$bi = -1;
		$ei = -1;
		$ni = -1;
	#	$pi = -1;
		foreach(split /\t/, $header) {	# index
			$bi = $i if $_ eq "BETA"; 	# beta
			$ni = $i if $_ eq "N_PH";		# sample size
			$ei = $i if $_ eq "SE";			# standard error
	#		$pi = $i if $_ eq "P_VALUE";
			$i++;
		}
		next if $bi < 0 || $ei < 0;
#		print $fileName,"\n" if $ni < 0;
#		next;
		%es=(); #standar error
		%bs=(); #betas
		%ns=(); #sample size
	#	%ps =(); #pvalue
	#	$minpval = 1;
	#	$snpcnt = 0;
		while(<IN>) {
			/^([^\t]+)/;
			next if !exists $ids{$1};
			s/[\r\n]//g;
			@d = split /\t/, $_;
			$es{$d[0]} = $d[$ei];
			$bs{$d[0]} = $d[$bi] * $dir{$d[0]};		# beta * direction
			$ns{$d[0]} = $d[$ni] if $ni >= 0;
	#		$ps{$d[0]} = $d[$pi];
	#		$snpcnt++;
	#		if ($d[$pi] < $minpval) { $minpval = $d[$pi]; }
		}
	#	next if $minpval > 0.01;
	#	$HEADER .= "\t$ds-$ph-ZN";
		#$HEADER .= "\t$ds-$ph-Z";
		if($ni >= 0) {	# sample size >= 0
			$HEADER .= "\t$ds-$ph-ZN";
			$I = 0;
			foreach(@ids) {					# why 258 ? -- mean of the sqrt(sample size)
			#	$OUT[$I++] .= "\t" . (exists $ns{$_} ? 258*$bs{$_}/$es{$_}/($ns{$_}**(1/3)) : ".");
				$OUT[$I++] .= "\t" . (exists $ns{$_} ? 258*$bs{$_}/$es{$_}/sqrt($ns{$_}) : "."); #print Z-score/sqrt(N)			# if there is no sample size info, then print "."
			#	$OUT[$I++] .= "\t" . (exists $ps{$_} ? $ps{$_} : ".");													 #print p-value
			#	$OUT[$I++] .= "\t" . (exists $bs{$_} ? $bs{$_} : ".");                          #print BETA
			#	$OUT[$I++] .= "\t" . (exists $ns{$_} ? $ns{$_} : ".");
			# $OUT[$I++] .= "\t" . (exists $ns{$_} ? $bs{$_}/$es{$_} : ".");
			}
		}
		else {	# when sample size info is missing
			$pheno = "$ds-$ph";
			next if !exists $pheno2samplesize{$pheno};
			$HEADER .= "\t$ds-$ph-ZN";
		#	die "No N for: $pheno\n" if !exists $pheno2samplesize{$pheno};
			$phenoNS = $pheno2samplesize{$pheno};
			$I = 0;
			foreach(@ids) {
			#	$OUT[$I++] .= "\t" . (exists $bs{$_} ? 258*$bs{$_}/$es{$_}/($phenoNS**(1/3)) : ".");
				$OUT[$I++] .= "\t" . (exists $bs{$_} ? 258*$bs{$_}/$es{$_}/sqrt($phenoNS) : "."); #print Z-score/sqrt(N)
			#	$OUT[$I++] .= "\t" . (exists $ps{$_} ? $ps{$_} : ".");														#print p-value
			#	$OUT[$I++] .= "\t" . (exists $bs{$_} ? $bs{$_} : ".");                           #print BETA
			#	$OUT[$I++] .= "\t" . (exists $pheno2samplesize{$pheno} ? $phenoNS : ".");
			# $OUT[$I++] .= "\t" . (exists $bs{$_} ? $bs{$_}/$es{$_} : ".");
			}
		}
	}
}



#exit;
print $HEADER,"\n";
foreach(@OUT) { print "$_\n" };
