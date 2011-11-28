#!/usr/bin/perl -w

####
#
# Script che produce il grafo delle isoforme (STDOUT) dal file GTF dei trascritti (STDIN)
# Produce anche in STDOUT i percorsi del grafo che corrispondono ai full-length transcripts contenuti nel GTF
# Attenzione: nel gtf deve essere presente un solo gene! La genomica in input deve sempre essere sullo strand +
#
####

use Time::HiRes;
use strict;
use Getopt::Long;

our $debug = 0;
my $contract_path;
my $wdir;
my $gtfFile;
my $gen_file;
my $outIsoformFile;
my $outIsoformFileRegions;
my $outIsoformFileGDL;
my $outIsoformPaths;
my $min_length;
my $force_strand;

GetOptions (
			'contract-path=i' => \$contract_path,
			'min-length=i' => \$min_length,
			'working-dir=s' => \$wdir,
            'gtf=s' => \$gtfFile,
            'genomic=s' => \$gen_file, 
            'isoforms=s' => \$outIsoformFile,
            'isoforms-regions=s' => \$outIsoformFileRegions,
            'isoforms-gdl=s' => \$outIsoformFileGDL,
            'isoforms-paths=s' => \$outIsoformPaths,
            'force-on-plus-strand=i' => \$force_strand,
            'debug' => \$debug,
           );
my $usage="Usage: perl BuildIsoformGraph.pl [options]
 --contract-path= <integer> (mandatory) type of path contraction
 							(0=codingregions; 1=blocks; 2=RNA-seq-graph-like)
 --min-length= <integer> (optional) minimum length for retaining a node (only for type of path contraction = 2)
 							    (default: 0)
 --working-dir=<file>		 (optional) path of the working directory
 								(default: ./)
 --gtf= <file>          	(optional) name of the input GTF file
 								(default: exons.gtf)
 --genomic= <file>           (optional) name of the genomic sequence file (5'3')
 								(default: genomic.txt)
 --isoforms= <file>          (optional) name of the output Isoform Graph file
 								(default: isoform-graph.txt)
 --isoforms-regions= <file>  (optional) name of the output Isoform Graph file with regions
 								(default: isoform-graph-regions.txt)
 --isoforms-gdl= <file>  	(optional) name of the output Isoform Graph file in GDL format
 								(default: isoform-graph.gdl)
 --isoforms-paths= <file>  	(optional) name of the output Isoform Paths file
 								(default: isoform-paths.txt)
 --force-on-plus-strand= <integer>  (optional) flag for forcing the isoform graph on plus strand
 								(1: the graph is given on plus strand; 0: the graph is given on the gene strand)
 								(default: 0)
 --debug                     emits debugging information
";

die $usage unless (defined $contract_path and $contract_path >= 0 and $contract_path <= 2);

$wdir="./" if (not defined $wdir or $wdir eq '' or $wdir eq ".");
print "Working directory: ", $wdir, "\n";
chop $wdir if(substr($wdir, length($wdir)-1, 1) eq "/");

$gtfFile="exons.gtf" if (not defined $gtfFile or $gtfFile eq '');
print "GTF file: ", $gtfFile, "\n";

$gen_file="genomic.txt" if (not defined $gen_file or $gen_file eq '');
print "Genomic sequence file: ", $gen_file, "\n";

$outIsoformFile="isoform-graph.txt" if (not defined $outIsoformFile or $outIsoformFile eq '');
print "Output isoform graph file: ", $outIsoformFile, "\n";

$outIsoformFileRegions="isoform-graph-regions.txt" if (not defined $outIsoformFileRegions or $outIsoformFileRegions eq '');
print "Output isoform graph file with regions: ", $outIsoformFileRegions, "\n";

$outIsoformFileGDL="isoform-graph.gdl" if (not defined $outIsoformFileGDL or $outIsoformFileGDL eq '');
print "Output GDL isoform graph file: ", $outIsoformFileGDL, "\n";

$outIsoformPaths="isoform-paths.txt" if (not defined $outIsoformPaths or $outIsoformPaths eq '');
print "Output isoform paths file: ", $outIsoformPaths, "\n";

print "Path contraction: ", $contract_path, "\n";
$min_length=0 if (not defined $min_length or $min_length < 0);
if($contract_path == 2){
	print "Nodes of length < $min_length will be discarded!\n";
}

$force_strand=0 if (not defined $force_strand or ($force_strand != 0 and $force_strand != 1));

$gtfFile=$wdir."/".$gtfFile;
$gen_file=$wdir."/".$gen_file;
$outIsoformFile=$wdir."/".$outIsoformFile;
$outIsoformFileRegions=$wdir."/".$outIsoformFileRegions;
$outIsoformFileGDL=$wdir."/".$outIsoformFileGDL;
$outIsoformPaths=$wdir."/".$outIsoformPaths;

for my $f (($gen_file, $gtfFile)){
    die "File $f does not exists\n" unless (-f $f);
}

open OUT, ">$outIsoformFile" or die "Could not open $outIsoformFile: $!\n";
open OUTREG, ">$outIsoformFileRegions" or die "Could not open $outIsoformFileRegions: $!\n";
open OUTGDL, ">$outIsoformFileGDL" or die "Could not open $outIsoformFileGDL: $!\n";
open OUTIP, ">$outIsoformPaths" or die "Could not open $outIsoformPaths: $!\n";

open GEN, "<", $gen_file or die "Could not read $gen_file: $!\n";
my $header=<GEN>;
die "Genomic file $gen_file not in FASTA format!" unless($header =~ m/^>/);
my $gen_seq="";
while(<GEN>){
	chomp;
	$gen_seq.=$_;
}
close GEN;
$gen_seq=uc($gen_seq);
my $gen_length=length($gen_seq);

open GTF, $gtfFile or die "Could not open $gtfFile: $!\n";

my $strand=1;
my %exon_hash=();	#key=exon_start, value=hash_ref1(key=exon_end, value=hash_ref2(key=transcript_ID, value=1))

my %full_lengths=(); #key=transcript_ID value list_ref(exon1_start, exon1_end, exon2_start, exon2_end, ...)
my %exon_map=();	#key=exon_start, value=hash_ref1(key=exon_end, value=hash_ref(key=node value=dimensione))
					#node è il nodo del grafo che risulta mappato all'esone
while(<GTF>){
	chomp;
	my @record_list=split(/\s/, $_);
	
	#print $record_list[2], "\n";
	#exit;
	
	if($record_list[2] eq "exon"){
		my $transcript_id_in_record=-1;
		my $k=0;
		my $found=0;
		while($k <= $#record_list && $found == 0){
			if($record_list[$k] eq "transcript_id"){
				$found=1;
			}
			$k++;
		}
		if($found == 0 || $k == $#record_list+1){
			die "No transcript_id field found in the GTF file!\n";
		}		
		
		my $exon_start=$record_list[3];
		my $exon_end=$record_list[4];
		
		if($record_list[6] eq "-"){
			$strand=-1;
		}
                
		#Inversione per gene -
		if($strand == -1 && $force_strand != 1){
			$exon_start=$gen_length-$record_list[4]+1;
			$exon_end=$gen_length-$record_list[3]+1;
		}
		else{
			$exon_start=$record_list[3];
			$exon_end=$record_list[4];
                }

                if(exists $full_lengths{$record_list[$k]}){
			if($strand == 1){
				push @{$full_lengths{$record_list[$k]}}, ($exon_start, $exon_end);
			}else{
				unshift @{$full_lengths{$record_list[$k]}}, ($exon_start, $exon_end);
			}
		}else{
			$full_lengths{$record_list[$k]}=[$exon_start, $exon_end];
		}
		my %exon_map1;
		if(exists $exon_map{$exon_start}){
			%exon_map1=%{$exon_map{$exon_start}};
			if(!exists $exon_map1{$exon_end}){
				my %exon_map2=();
				$exon_map1{$exon_end}=\%exon_map2;
			}
		}
		else{
			%exon_map1=();
			my %exon_map2=();
			$exon_map1{$exon_end}=\%exon_map2;
		}
		$exon_map{$exon_start}=\%exon_map1;
                
		my %exon_hash1;
		if(exists $exon_hash{$exon_start}){
			#print "\tEsiste start\n";
			%exon_hash1=%{$exon_hash{$exon_start}};
			my %exon_hash2;
			if(exists $exon_hash1{$exon_end}){
				#print "\tEsiste end\n";
				%exon_hash2=%{$exon_hash1{$exon_end}};
			}
			else{
				#print "\tNON Esiste end\n";
				%exon_hash2=();
			}
			#$exon_hash2{$record_list[9]}=1;
			$exon_hash2{$record_list[$k]}=1;
			
			$exon_hash1{$exon_end}=\%exon_hash2;
			
		}
		else{
			#print "\tNON Esiste start\n";
			%exon_hash1=();
			my %exon_hash2=();
			#$exon_hash2{$record_list[9]}=1;
			$exon_hash2{$record_list[$k]}=1;
			
			$exon_hash1{$exon_end}=\%exon_hash2;
		}
		$exon_hash{$exon_start}=\%exon_hash1;
	}
}

close GTF;

#foreach my $transcr_ID(keys %full_lengths){
#	print $transcr_ID;
#	print "\t", join("-", @{$full_lengths{$transcr_ID}}), "\n";
#}


print "The gene strand is ", ($strand == 1)?("PLUS"):("MINUS");
my $cong=($strand == -1 && $force_strand == 1)?(" but"):(" and");
print $cong;
print " the isoform graph is ", ($force_strand == 1)?("forced on PLUS"):("on the gene"), " strand\n";


if($strand == -1){
	$gen_seq=reverse_complement($gen_seq);
}

my %B_hash=();	#key=position, value=0|1|2 (0:leftend, 1:rightend, 2:both)
foreach my $key(keys %exon_hash){
	if(exists $B_hash{$key}){
		if($B_hash{$key} == 1){
			$B_hash{$key}=2;
		}
	}
	else{
		$B_hash{$key}=0;
	}
	my %exon_hash1=%{$exon_hash{$key}};
	foreach my $key1(keys %exon_hash1){
		if(exists $B_hash{$key1}){
			if($B_hash{$key1} == 0){
				$B_hash{$key1}=2;
			}
		}
		else{
			$B_hash{$key1}=1;
		}
		
		#my %exon_hash2=%{$exon_hash1{$key1}};
		#print $key, "-", $key1, "\n";
		#my @keys_id=keys %exon_hash2;
		#print "\t".join(",", @keys_id), "\n";
	}
}

my @B_list=();	#[genomic_position, 0|1]	0:leftend, 1:rightend
foreach my $key(sort {$a <=> $b} (keys %B_hash)){
	if($B_hash{$key}==2){
		push @B_list, [$key, 0];
		push @B_list, [$key, 1];
	}
	else{
		push @B_list, [$key, $B_hash{$key}];
	}
	#print "POS=".$key." FLAG=".$B_hash{$key}."\n";
}

#@G_list=(...,ref,...)
#ref: reference to list=(reference to l1, f1, f2, trs_hash_ref ,f3:contracted, list_of_out_ref, list_of_in_ref) ==> list of no-variation regions
#l1=(l1,r1,l2,r2,...) where li-ri is a region
#f1: flag of the block left type ==> 0 for start and 1 for end
#f2: flag of the block right type ==> 0 for start and 1 for end
#f1 and f2 after contraction may be not consistent
#trs_hash_ref: reference to the hash of the supporting transcripts
#f3: contraction flag ==> 1 if the nodes was contracted, 2 if was deleted
#list_of_out_ref: reference to the list of outcoming nodes (each node is represented by the reference ref in G_list)
#list_of_in_ref: reference to the list of incoming nodes (each node is represented by the reference ref in G_list)
my @G_list=();

for(my $p=0; $p < $#B_list; $p++){
	my @pos1=@{$B_list[$p]};
	my @pos2=@{$B_list[$p+1]};
	
	my @region=();
	my $add_region=0;
	#Region of type start-* ==> no-variation region
	if($pos1[1] == 0){
		$add_region=1;
		my @chunks=();
		push @chunks, $pos1[0];
		if($pos2[1] == 0){
			push @chunks, $pos2[0]-1;
		}
		else{
			push @chunks, $pos2[0];
		}
		push @region, \@chunks;
		push @region, $pos1[1], $pos2[1];
	}
	else{
		#Region of type end-end ==> no-variation region
		if($pos1[1] == $pos2[1]){
			$add_region=1;
			my @chunks=();
			push @chunks, $pos1[0]+1;
			push @chunks, $pos2[0];
			push @region, \@chunks;
			push @region, $pos1[1], $pos2[1];
		}
		else{
			#The region is not null
			if($pos2[0]-$pos1[0]-2 >= 0){
				my @exon_starts=sort {$a <=> $b} keys %exon_hash;
				my $i=0;
				my $found=0;
				while($i <= $#exon_starts && ($exon_starts[$i] <= $pos1[0] && $found == 0)){
					my %exon_hash1=%{$exon_hash{$exon_starts[$i]}};
					my @exon_ends=sort {$b <=> $a} keys %exon_hash1;
					if($exon_ends[0] >= $pos2[0]){
						$found=1;	
					}
					$i++;
				}
				if($found == 1){
					$add_region=1;
					my @chunks=();
					push @chunks, $pos1[0]+1;
					push @chunks, $pos2[0]-1;
					push @region, \@chunks;
					push @region, $pos1[1], $pos2[1];
				}
			}
		}
	}
	if($add_region == 1){
		my %region_trs=();
		my @exon_starts=sort {$a <=> $b} keys %exon_hash;
		my $i=0;
		while($i <= $#exon_starts && $exon_starts[$i] <= ${$region[0]}[0]){
			my %exon_hash1=%{$exon_hash{$exon_starts[$i]}};
			my @exon_ends=sort {$b <=> $a} keys %exon_hash1;
			my $j=0;
			while($j <= $#exon_ends && ($exon_ends[$j] >= ${$region[0]}[1])){
				my %exon_hash2=%{$exon_hash1{$exon_ends[$j]}};
				my @ids=keys %exon_hash2;
				
				#push @region_trs, @ids;
				foreach my $id(@ids){
					$region_trs{$id}=1;
				}
				
				$j++;
			}
			$i++;
		}
		
		push @region, \%region_trs;	
		push @region, 0;
		push @region, [];
		push @region, [];
		push @G_list, \@region;
	}
}

my $adj_matx;
for(my $i=0; $i <= $#G_list; $i++){
	for(my $j=0; $j <= $#G_list; $j++){
		$adj_matx->[$i][$j]=0;
	}
}

for(my $i=0; $i <= $#G_list; $i++){
	my @region=@{$G_list[$i]};
	
	#my @chunks=@{$region[0]};
	#print $i, "->", $chunks[0]."-".$chunks[1]."-".$region[1]."-".$region[2]."\n";
	
	if($region[1] == $region[2]){
		if($region[1] == 0){
			$adj_matx->[$i][$i+1]=1;
		}
		else{
			$adj_matx->[$i-1][$i]=1;
		}
	}
	else{
		if($region[1] == 1){
			$adj_matx->[$i][$i+1]=1;
			$adj_matx->[$i-1][$i]=1;
		}
	}
}

for(my $i=0; $i < $#G_list; $i++){
	my @region=@{$G_list[$i]};

	my @chunks=@{$region[0]};
	
	#print "Node ".$i.": ".$chunks[0]."-".$chunks[1]."-".$region[1]."-".$region[2]."\n";

	my %hash_region_trs=%{$region[3]};

	#print "\t".join(",", keys %hash_region_trs)."\n";

	if($region[1] == 0 && $region[2] == 0){
		if($adj_matx->[$i][$i+1] != 1){
			die "Region ".$chunks[0]."-".$chunks[1]." must be linked to the next one!\n";
		}
	}
	else{
		if($region[1] == 1 && $region[2] == 0){
			if($adj_matx->[$i][$i+1] != 1 || $adj_matx->[$i-1][$i] != 1){
				die "Region ".$chunks[0]."-".$chunks[1]." must be linked to the next one!\n";
			}
		}
		else{
			if($region[1] == 1 && $region[2] == 1){
				if($adj_matx->[$i-1][$i] != 1){
					die "Region ".$chunks[0]."-".$chunks[1]." must be linked to the previous one!\n";
				}
			}

			for(my $j=$i+1; $j <= $#G_list; $j++){
				
				my @next_region=@{$G_list[$j]};
				my %next_region_trs=%{$next_region[3]};
				my $count_inters=0;
				foreach my $tr(keys %next_region_trs){
					
					if(exists $hash_region_trs{$tr}){
						if($hash_region_trs{$tr} == 1){
							$hash_region_trs{$tr}=0;
							$count_inters++;
						}
					}
				}				
			
				if($adj_matx->[$i][$j] == 1){
					if($j > $i+1 || $count_inters == 0){
						die "Problem 2!\n";
					}
				}
				else{
					if($count_inters > 0){
						$adj_matx->[$i][$j]=1;
					}
				}
			}
		}		
	}
}

#Add the adjacency list
for(my $i=0; $i <= $#G_list; $i++){
	my @adj_node=();
	push @adj_node, $G_list[$i];
	
	for(my $j=0; $j <= $#G_list; $j++){
		if($adj_matx->[$i][$j] == 1){
			push @{${$G_list[$i]}[5]}, $G_list[$j];
			push @{${$G_list[$j]}[6]}, $G_list[$i];
		}	
	}
}

#Path contraction
my $i=0;
while($contract_path >= 1 && $i <= $#G_list){
	my $region_ref=$G_list[$i];
	#Non-contracted node
	if(${$region_ref}[4] == 0){
		#Region of type *-end	
		if($contract_path == 2 || ${$region_ref}[2] == 1){
			#print "Region ".${${$region_ref}[0]}[0]."-".${${$region_ref}[0]}[1]."\n";
			my $stop=0;
			while($stop == 0){
				my @outs=@{${$region_ref}[5]};
				#Just one adj
				if($#outs == 0){
					my $out_ref=$outs[0];
					if(${$out_ref}[4] == 1){
						die "Problem in contraction flag!\n";
					}
					#Indegree = 1
					if(scalar(@{${$out_ref}[6]}) == 1){
						#Type start-*
						if($contract_path == 2 || ${$out_ref}[1] == 0){
							
							#print "\tto region ".${${$out_ref}[0]}[0]."-".${${$out_ref}[0]}[1]."\n";
														
							my $link_nodes=1;
							
							#Check if the transcript sets are the same
							if($contract_path == 1){
								my %tr_first=%{${$region_ref}[3]};
								my %tr_second=%{${$out_ref}[3]};
								
								#RAFFA
								my @chunks_from=@{${$region_ref}[0]};
								my $right_from=$chunks_from[1];
								my @chunks_to=@{${$out_ref}[0]};
								my $left_to=$chunks_to[0];
								my $p=0;
								while($p <= $#G_list && $link_nodes == 1){
									my $region_ref_cfr=$G_list[$p];
									my @chunks_cfr=@{${$region_ref_cfr}[0]};
									my $end_cfr=$chunks_cfr[0];
									if($end_cfr > $right_from && $end_cfr < $left_to){
										$link_nodes=0;
									}
									$p++;
								}
								
								#RAFFA
								if($link_nodes == 1){
							
									my @keys_first=keys %tr_first;
									foreach my $key(@keys_first){
										if(!exists $tr_second{$key}){
											$link_nodes=0;
										}
									}
								#RAFFA
								}
								
								if($link_nodes == 1){
									my @keys_second=keys %tr_second;
									foreach my $key(@keys_second){
										if(!exists $tr_first{$key}){
											$link_nodes=0;
										}
									}
								}
							}
							
							if($link_nodes == 1){													
								${$out_ref}[4]=1; #Contract node
								${$G_list[$i]}[5]=${$out_ref}[5];
								push @{${$G_list[$i]}[0]}, @{${$out_ref}[0]};
								
								foreach my $out_n(@{${$out_ref}[5]}){
									for(my $w=0; $w < scalar(@{${$out_n}[6]}); $w++){
										if(${${$out_n}[6]}[$w] == $out_ref){
											${${$out_n}[6]}[$w]=$G_list[$i];
										}
									}
								}
							
								if($contract_path == 1 && ${$out_ref}[2] == 0){
									$stop=1;	
								}
								else{
									$region_ref=$out_ref;	
								}
							}
							else{
								$stop=1;	
							}
						}
						else{
							$stop=1;	
						}
					}
					else{
						$stop=1;	
					}
				}
				else{
					$stop=1;	
				}
			}
			
		}
	}
	$i++;
}

#Eliminazione dei nodi iniziali/terminali (cioé con solo un arco uscente/entrante) di lunghezza < $min_length
#[list_ref:list of regions, 0|1, 0|1, trs_hash_ref ,0|1|2:contracted|deleted, list_of_out_ref, list_of_in_ref]	list of no-variation regions
my $k=0;
while($k <= $#G_list){
	if(${$G_list[$k]}[4] == 0){
		my @outs=@{${$G_list[$k]}[5]};
		my @ins=@{${$G_list[$k]}[6]};
		if((scalar(@outs) == 1 && scalar(@ins) == 0) || (scalar(@outs) == 0 && scalar(@ins) == 1)){
			my @regions=@{${$G_list[$k]}[0]};
			die "Problem 1!\n" if(scalar(@regions) % 2 != 0);
			my $region_length=0;
			my $p=0;
			while($p <= $#regions){
				$region_length+=($regions[$p+1]-$regions[$p]+1);
				$p+=2;
			}
			if($region_length < $min_length){
				${$G_list[$k]}[4]=2; #Metto a deleted
				if(scalar(@outs) == 1){	#nodo iniziale
					my $out_ref=$outs[0];
					my @out_region=@{$out_ref};
					my @n_ins=@{$out_region[6]};
					
					my @new_n_ins=();
					foreach my $n(@n_ins){
						if($n != $G_list[$k]){
							push @new_n_ins, $n;
						}
					}
					#Non aggiorno i nodi entranti del nodo uscente perché risolvo in fase di stampa
					${${${$G_list[$k]}[5]}[0]}[6]=\@new_n_ins;
					${$G_list[$k]}[5]=[];
				}
				else{	#nodo terminale
					my $in_ref=$ins[0];
					my @in_region=@{$in_ref};
					my @n_outs=@{$in_region[5]};
				
					my @new_n_outs=();
					foreach my $n(@n_outs){
						if($n != $G_list[$k]){
							push @new_n_outs, $n;
						}
					}
					#Aggiorno i nodi uscenti del nodo entrante perché risolvo in fase di stampa
					${${${$G_list[$k]}[6]}[0]}[5]=\@new_n_outs;			
					${$G_list[$k]}[6]=[];
				}
			}
		}
	}
	$k++;
}

#Eliminazione degli altri nodi di lunghezza < $min_length
#[list_ref:list of regions, 0|1, 0|1, trs_hash_ref ,0|1:contracted, list_of_out_ref, list_of_in_ref]	list of no-variation regions
#$k=0;
#while($k <= $#G_list){
#	if(${$G_list[$k]}[4] == 0){
#		my @regions=@{${$G_list[$k]}[0]};
#		die "Problem 1!\n" if(scalar(@regions) % 2 != 0);
#		my $region_length=0;
#		my $p=0;
#		while($p <= $#regions){
#			$region_length+=($regions[$p+1]-$regions[$p]+1);
#			$p+=2;
#		}
#		if($region_length < $min_length){
#			my @outs=@{${$G_list[$k]}[5]};
#			my @ins=@{${$G_list[$k]}[6]};
#			
#			#${$G_list[$k]}[4]=2;	#Metto a deleted
#			my $out_ref=0;
#			my $in_ref=0;
#		}
#	}
#	$k++;
#}

#Node labelling
$k=0;
my $label=1;
my $label_contr=1;
while($k <= $#G_list){
	if(${$G_list[$k]}[4] == 0){
		push @{$G_list[$k]}, $label;
		$label++;
	}
	else{
		push @{$G_list[$k]}, "Contracted ".$label_contr;
		$label_contr++;
	}
	$k++;
}

print OUTGDL "graph: \{\n";
print OUTGDL "node.shape\t: circle\n";
print OUTGDL "node.color\t: blue\n";
print OUTGDL "node.height\t: 80\n";
print OUTGDL "node.width\t: 80\n";

my $max_index_node=1;
foreach my $node_ref(@G_list){
	my @region=@{$node_ref};

	#Non-contracted node
	if($region[4] == 0){
		my @chunks=@{$region[0]};
		my $title="";
		my $i=0;
		my $dimension=0;
		my $chunk_seq="";
		while($i <= $#chunks){
			$title.=$chunks[$i]."-".$chunks[$i+1].";";
			$dimension=$dimension+$chunks[$i+1]-$chunks[$i]+1;
			
			$chunk_seq=$chunk_seq.substr($gen_seq, $chunks[$i]-1, $chunks[$i+1]-$chunks[$i]+1);
                        
                        #print "Chunk ", $chunks[$i]."-".$chunks[$i+1], "(of node ", $region[7], ")\n";
			
			#Mapping agli esoni
			foreach my $exon_start(keys %exon_map){
				if($chunks[$i] >= $exon_start){
					my %hash_temp=%{$exon_map{$exon_start}};
					foreach my $exon_end(keys %hash_temp){
						if($chunks[$i+1] <= $exon_end){
							my %hash_temp2=%{$hash_temp{$exon_end}};
							if(exists $hash_temp2{$region[7]}){
								$hash_temp2{$region[7]}=$hash_temp2{$region[7]}+$dimension;
							}else{
								$hash_temp2{$region[7]}=$dimension;
							}
							
							$hash_temp{$exon_end}=\%hash_temp2;
							$exon_map{$exon_start}=\%hash_temp;	
							#print "\tmap to exon $exon_start-$exon_end\n";
						}					
					}									
				}
			}
			
			$i+=2;	
                }
                
		$max_index_node=($max_index_node < $region[7])?($region[7]):($max_index_node);
		print OUTIP $region[7], "(", length($chunk_seq), ")\t";
                
		print OUT "node\#".$region[7]." ".$chunk_seq."\n";
		print OUTREG "node\#".$region[7]." ".$title."\n";
		
		#print OUTGDL "\tnode: \{\n\t\ttitle: "."\"".$title."\""."\n\t\}\n";
		print OUTGDL "\tnode: \{\n\t\ttitle: "."\"".$region[7]."\"\n";
		print OUTGDL "\t\tlabel: "."\"".$region[7]." - ".$dimension."\""."\n\t\}\n";
	}
}

print OUTIP "\n";

my $edge_index=1;
foreach my $node_ref(@G_list){
	my @region=@{$node_ref};

	#Non-contracted node
	if($region[4] == 0){
		my @chunks=@{$region[0]};
		my $title="";
		my $i=0;
		while($i <= $#chunks){
			$title.=$chunks[$i]."-".$chunks[$i+1].",";
			$i+=2;	
		}
		my @outs=@{$region[5]};
		
		foreach my $adj_ref(@outs){
			my @dest_region=@{$adj_ref};
			
			my @dest_chunks=@{$dest_region[0]};
			my $dest_title="";
			my $j=0;
			while($j <= $#dest_chunks){
				$dest_title.=$dest_chunks[$j]."-".$dest_chunks[$j+1].",";
				$j+=2;	
			}
			
			print OUT "edge\#".$edge_index." ".$region[7].";".$dest_region[7]."\n";
			print OUTREG "edge\#".$edge_index." ".$region[7].";".$dest_region[7]."\n";

			#print OUTGDL "\tedge: \{\n\t\tsource: "."\"".$title."\""."\n";
			#print OUTGDL "\t\ttarget: "."\"".$dest_title."\""."\n\t\}\n";
			print OUTGDL "\tedge: \{\n\t\tsource: "."\"".$region[7]."\""."\n";
			print OUTGDL "\t\ttarget: "."\"".$dest_region[7]."\""."\n\t\}\n";
			
			$edge_index++;
		}
	}
}
print OUTGDL "\}\n";

close OUT;
close OUTREG;
close OUTGDL;

#Stampa percorsi compatibili con i full-lengths
foreach my $transcript(keys %full_lengths){
	
	my @exons=@{$full_lengths{$transcript}};
	my $i=0;
	
	my %path=();
	
	while($i <= $#exons){
		my $exon_start=$exons[$i];
		my $exon_end=$exons[$i+1];
		
		my %hash_map1=%{$exon_map{$exon_start}};
		my %hash_map2=%{$hash_map1{$exon_end}};
		my @nodes=keys %hash_map2;
		
		foreach my $node(@nodes){
			$path{$node}=$hash_map2{$node};
		}
				
		$i+=2;
	}
	
	my @nodes=sort {$a <=> $b} keys %path;
	
	my $prev=1;
	my $tabs=-1;
	my $tab="\t";
	foreach my $i(0..$#nodes){
		$tabs=$nodes[$i]-$prev;
		print OUTIP $tab x $tabs;
		print OUTIP "x";
		$prev=$nodes[$i];
	}
	$tabs=$max_index_node-$prev+1;
	
	print OUTIP $tab x $tabs;
	print OUTIP $transcript, "\n";
}

close OUTIP;

exit;

sub reverse_complement{
	my $gen_to_be_rev=shift;
	
	my $rev_seq=reverse $gen_to_be_rev;
	my $new_seq="";

	my $i=0;
	
	$rev_seq =~ tr/acgtACGT/tgcaTGCA/;
	return $rev_seq;
}

sub complement{
	my $car=shift;
	
	my $new_car;
	if($car eq "a"){
		$new_car="t";
	}elsif($car eq "A"){
		$new_car="T";
	}elsif($car eq "t"){
		$new_car="a";
	}elsif($car eq "T"){
		$new_car="A";
	}elsif($car eq "c"){
		$new_car="g";
	}elsif($car eq "C"){
		$new_car="G";
	}elsif($car eq "g"){
		$new_car="c";
	}elsif($car eq "G"){
		$new_car="C";
	}else{
		$new_car=$car;
	}
	
	return $new_car;
	
}
