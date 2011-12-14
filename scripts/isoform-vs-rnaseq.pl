#!/usr/bin/perl -w

use strict;
use POSIX;
use Getopt::Long;

our $debug = 0;
my $wdir;
my $fileIsoformG;
my $RNASeqG;
my $out_file;
my $log_file;
my $f_length;
my $f_offset;
my $mapping_threshold;
my $gnuplot_file = "gnuplot-script.gnu";

GetOptions (
			'working-dir=s' => \$wdir,
            'isoforms=s' => \$fileIsoformG,
            'rna-seq=s' => \$RNASeqG,
            'output=s' => \$out_file,
            'log-file=s' => \$log_file,
            'fingerp-length=i' => \$f_length,
            'fingerp-offset=i' => \$f_offset,
            'mapping-threshold=i' => \$mapping_threshold,
            'gnuplot-file=s' => \$gnuplot_file,
            'debug' => \$debug,
           );
my $usage="Usage: perl isoform-vs-rnaseq2.pl [options]
 --fingerp-length= <integer> (mandatory) length for fingerprint the GF nodes (>= 1)
 --working-dir=<file>		 (optional) path of the working directory
 								(default: ./)
 --isoforms= <file>          (optional) name of the Isoform Graph file
 								(default: isoform-graph.txt)
 --rna-seq= <file>           (optional) name of the RNA-seq Graph file
 								(default: RNA-seq-graph.txt)
 --output=<file>			 (optional) name of the output file
 								(default: GRvsGFoutput.txt)
 --log-file=<file>			 (optional) name of the log file
 								(default: GRvsGFlog.txt)
 --mapping-threshold	     (optional) maximum threshold for border mapping (>= 0)
 								(default: 2)
 --gnuplot-file=<file>       name of the Gnuplot script file
                                         (default: gnuplot-script.gnu)
 --debug                     emits debugging information
";

die $usage unless (defined $f_length and $f_length >= 1);

#$f_offset=-1 unless(defined $f_offset and $f_offset >= 0);
$f_offset=-1;

$mapping_threshold=2 unless(defined $mapping_threshold and $mapping_threshold >= 0);
  
$wdir="./" if (not defined $wdir or $wdir eq '' or $wdir eq ".");
print "Working directory: ", $wdir, "\n";
chop $wdir if(substr($wdir, length($wdir)-1, 1) eq "/");
            
$fileIsoformG="isoform-graph.txt" if (not defined $fileIsoformG or $fileIsoformG eq '');
print "Isoform graph file: ", $fileIsoformG, "\n";

$RNASeqG="RNA-seq-graph.txt" if (not defined $RNASeqG or $RNASeqG eq '');
print "RNA-seq graph file: ", $RNASeqG, "\n";

$out_file="GRvsGFoutput.txt" if (not defined $out_file or $out_file eq '');
print "Output file: ", $out_file, "\n";

$log_file="GRvsGFlog.txt" if (not defined $log_file or $log_file eq '');
print "Log file: ", $log_file, "\n";

print "Fingerprinting length: ", $f_length, "\n"; 
print "Fingerprinting offset: ", ($f_offset == -1)?("all the sequence"):($f_offset), "\n"; 
print "Max mapping threshold: ", $mapping_threshold, "\n"; 

print "Gene\t#GRnodes\t#GFnodes\t#GRGFnodes\t#GRnodesINGF\t#GFnodesINGR\t";
print "#GRedges\t#GFedges\t#GRGFedges\t#GRedgesINGF\tGFedgesINGR\t";
print "\n";

my $dirh;
opendir ($dirh, "$wdir") or die "Cannot open $wdir!\n";
my @lists=readdir($dirh);
foreach my $dir(@lists){
	
	if($dir =~ m/out_(.+)/){
	#if($dir =~ m/out_(STAG2)/){
	#if($dir =~ m/(.+)_no_error/){
	#if($dir =~ m/(.+)_real/){
		
		my $current_dir=$wdir."/".$dir;
		#my $current_dir=$wdir."/".$dir."/corretta";
		
		#print "Current ", $current_dir, "\n";
		my $gene=$1;
		
		print STDERR $gene, "\t";
		
		my $current_fileIsoformG=$current_dir."/".$fileIsoformG;
		my $current_RNASeqG=$current_dir."/".$RNASeqG;
		my $current_out_file=$current_dir."/".$out_file;
		my $current_log_file=$current_dir."/".$log_file;
                my $current_gnuplot_file=$current_dir."/".$gnuplot_file;
		
		my $next_iter=0;
		
		for my $f (($current_fileIsoformG, $current_RNASeqG)){
		    #die "File $f does not exists\n" unless (-f $f);
		   if(!(-f $f)){
		   		print STDERR "File $f does not exists\n";
		    	$next_iter=1;
		    }	
		    
		}
		
		if($next_iter){
			print STDOUT $gene, "\tnot found\n";
			next
		};
		
		open OUT, ">$current_out_file" or die "Could not open $current_out_file: $!\n";
		open LOG, ">$current_log_file" or die "Could not open $current_log_file: $!\n";
                open GNUPLOT, ">$current_gnuplot_file" or die "Could not open $current_gnuplot_file: $!\n";
		
		open GF, "<", $current_fileIsoformG or die "Could not read $current_fileIsoformG: $!\n";
		
		#Build GF node list and GF adj matrix
		#GF node list ==> (list_ref1); list1 ==> (list_ref2); list2 ==> list2[0]: sequence
		
		my $GF_adj_matx;
		my @GF_node_list=();
		
		my $count_GF_node=1;
		my $count_GF_edge=1;
		
		my @GF_edges=();
		
		my $GF_small_blocks_64=0;
		my $GF_small_blocks_32=0;
		my $GF_small_blocks_16=0;
		my $GF_min_length=-1;
		
		while(<GF>){
			chomp;
			my @record=split(/ /, $_);
			die "Wrong node\edge format in $current_fileIsoformG!\n" unless(scalar(@record) == 2);
			if($record[0] =~ /node\#(\d+)/i){
				my $node_index=$1;
				die "Wrong node indexing in $current_fileIsoformG!\n" unless ($node_index == $count_GF_node);
				my @GF_node=();
				
				$record[1]=uc($record[1]);
				push @GF_node, $record[1];
				
				if(length($record[1]) < 64){
					$GF_small_blocks_64++;
					if(length($record[1]) < 32){
						$GF_small_blocks_32++;
						if(length($record[1]) < 16){
							$GF_small_blocks_16++;							
						}
					}
					if($GF_min_length == -1 || $GF_min_length > length($record[1])){
						$GF_min_length=length($record[1]);
					}
				}	
				push @GF_node_list, \@GF_node;
				$count_GF_node++;
			}
			elsif($record[0] =~ /edge\#(\d+)/i){
				my $edge_index=$1;
				die "Wrong edge indexing in $current_fileIsoformG!\n" unless ($edge_index == $count_GF_edge);
				push @GF_edges, $record[1];
				$count_GF_edge++;
					
			}
		}
		
		close GF;
		
		print STDERR $GF_small_blocks_64, "\t", $GF_small_blocks_32, "\t", $GF_small_blocks_16, "\t", $GF_min_length;
		print STDERR "\n";
		
		$count_GF_node-=1;
		$count_GF_edge-=1;
		
		print LOG "GF_GRAPH \#nodes=$count_GF_node, \#edges=$count_GF_edge\n";
		
		for(my $i=0; $i < $count_GF_node; $i++){
			for(my $j=0; $j < $count_GF_node; $j++){
				$GF_adj_matx->[$i][$j]=0;
			}
		}
		
		foreach my $edge_str(@GF_edges){
			my @edge=split(/;/, $edge_str);
			die "\tWrong edge format in $fileIsoformG!\n" unless (scalar(@edge) == 2);
			my $start=$edge[0]-1;
			my $end=$edge[1]-1;
			$GF_adj_matx->[$start][$end] = 1;
		}
		
		#PRINT @GF_node_list
		print_node_list(\@GF_node_list, 1, \*LOG);
		
		#PRINT $GF_adj_matx
		print_adj_matrix($GF_adj_matx, $count_GF_node, 1, \*LOG);
		
		#Build left and right fingerprint for GF
		#key=nt sequence, value=list_ref; list of ref to (node, offset)
		my $GF_finger_ref=get_fingerprinting(\@GF_node_list, $f_length, $f_offset, $debug);
		my %GF_left_fingerp=%{${$GF_finger_ref}[0]};
		my %GF_right_fingerp=%{${$GF_finger_ref}[1]};
		
		#PRINT %GF_left_fingerp
		print_fingerprinting(\%GF_left_fingerp, 1, \*LOG);
		
		#PRINT %GF_right_fingerp
		print_fingerprinting(\%GF_right_fingerp, 1, \*LOG);
		
		open GR, "<", $current_RNASeqG or die "Could not read $current_RNASeqG: $!\n";
		
		#Build GR node list and GR adj matrix
		#GR node list ==> (list_ref1); list1 ==> (list_ref2); list2 ==> list2[0]: nt sequence
		
		my $GR_adj_matx;
		my @GR_node_list=();
		
		my $count_GR_node=1;
		my $count_GR_edge=1;
		
		my @GR_edges=();
		
		while(<GR>){
			chomp;
			my @record=split(/ /, $_);
			die "Wrong node\edge format in $RNASeqG!\n" unless(scalar(@record) == 2);
			if($record[0] =~ /node\#(\d+)/i){
				
				$record[1]=uc($record[1]);
				
				my $node_index=$1;
				die "Wrong node indexing in $RNASeqG!\n" unless ($node_index == $count_GR_node);
				die "Fingerprint length too small wrt node $count_GR_node!\n" unless($f_length <= length($record[1]));
				my @GR_node=();
				chop $record[1] if(substr($record[1], length($record[1])-1, 1) eq ";");
				push @GR_node, $record[1];
				push @GR_node_list, \@GR_node;
				$count_GR_node++;
			}
			elsif($record[0] =~ /edge\#(\d+)/i){
				my $edge_index=$1;
				die "Wrong edge indexing in $RNASeqG!\n" unless ($edge_index == $count_GR_edge);
				push @GR_edges, $record[1];
				$count_GR_edge++;
					
			}
		}
		
		close GR;
		
		$count_GR_node-=1;
		$count_GR_edge-=1;
		
		print LOG "GR_GRAPH \#nodes=$count_GR_node, \#edges=$count_GR_edge\n";
		
		for(my $i=0; $i < $count_GR_node; $i++){
			for(my $j=0; $j < $count_GR_node; $j++){
				$GR_adj_matx->[$i][$j]=0;
			}
		}
		
		foreach my $edge_str(@GR_edges){
			my @edge=split(/;/, $edge_str);
			die "\tWrong edge format in $RNASeqG!\n" unless (scalar(@edge) == 2);
			my $start=$edge[0]-1;
			my $end=$edge[1]-1;
			$GR_adj_matx->[$start][$end] = 1;
		}
		
		#PRINT @GR_node_list
		print_node_list(\@GR_node_list, 0, \*LOG);
		
		#PRINT $GR_adj_matx
		print_adj_matrix($GR_adj_matx, $count_GR_node, 0, \*LOG);
		
		#Build left and right fingerprint for GR
		#key=nt sequence, value=list_ref; list of ref to (node, offset)
		my $GR_finger_ref=get_fingerprinting(\@GR_node_list, $f_length, $f_offset, $debug);
		my %GR_left_fingerp=%{${$GR_finger_ref}[0]};
		my %GR_right_fingerp=%{${$GR_finger_ref}[1]};
		
		#PRINT %GR_left_fingerp
		print_fingerprinting(\%GR_left_fingerp, 0, \*LOG);
		
		#PRINT %GR_right_fingerp
		print_fingerprinting(\%GR_right_fingerp, 0, \*LOG);
		
		#Map GR nodes to GF nodes (in @GR_node_list, two more list references are added to each node. These
		#lists contain (1) all the GF nodes which a GR node maps to, and for each mapping GF node,
		#the position, on the sequence of the GF node, of all the occurrences of the sequence of
		#the GR node (see %hash_temp). (2) All the cut mappings.
		map_nodes(\@GR_node_list, \@GF_node_list, \%GF_left_fingerp, \%GF_right_fingerp, $GF_adj_matx,
			$f_offset, $debug, \*LOG);
		    
		#Map GF nodes to GR nodes (in @GF_node_list, two more list references are added to each node. These
		#lists contain (1) all the GR nodes which a GF node maps to, and for each mapping GR node,
		#the position, on the sequence of the GR node, of all the occurrences of the sequence of
		#the GF node (see %hash_temp). (2) All the cut mappings.
		map_nodes(\@GF_node_list, \@GR_node_list, \%GR_left_fingerp, \%GR_right_fingerp, $GR_adj_matx,
			$f_offset, $debug, \*LOG);
		
		#Update GR nodes with list of contained GF nodes
		update_with_contained_nodes(\@GR_node_list, \@GF_node_list, $debug);
		
		#Update GF nodes with list of contained GR nodes    	
		update_with_contained_nodes(\@GF_node_list, \@GR_node_list, $debug);
			
		print "$gene\t";
		
		#Rimuove il mapping per i nodi piccoli (modifica introdotta il 13 dic 2011)
		my $GF_node_list_ref1=mapremove_micro_nodes(\@GF_node_list, 5, 1, \*LOG);
		my $GR_node_list_ref1=mapremove_micro_nodes(\@GR_node_list, 5, 0, \*LOG);
		@GF_node_list=@{$GF_node_list_ref1};
		@GR_node_list=@{$GR_node_list_ref1};
		
		#Sistema i nodi piccoli (modifica introdotta il 12 dic 2011)
		my $GF_node_list_ref2=handle_small_nodes(\@GF_node_list, 10, 1, \*LOG);
		my $GR_node_list_ref2=handle_small_nodes(\@GR_node_list, 10, 0, \*LOG);
		@GF_node_list=@{$GF_node_list_ref2};
		@GR_node_list=@{$GR_node_list_ref2};		
		
		#PRINT GR to GF node mapping
		print OUT "GR to GF node mapping\n";
		my $GRtoGF_match_ref=print_node_mapping(\@GR_node_list, 0, \*OUT, $gene);
		#Due nodi sono in relazione di match se sono $mapping_threshold-trim equivalenti

		#PRINT GF to GR node mapping
		print OUT "GF to GR node mapping\n";
		my $GFtoGR_match_ref=print_node_mapping(\@GF_node_list, 1, \*OUT, $gene);
		#Due nodi sono in relazione di match se sono $mapping_threshold-trim equivalenti

		my %GRtoGF_match_hash=%{$GRtoGF_match_ref};
		my %GFtoGR_match_hash=%{$GFtoGR_match_ref};

		#Solo per stampa a video
		my $GRtoGF_match=0;	#Biunivoci
		my $GRtoGF_map=0;
		foreach my $node(keys %GRtoGF_match_hash){			
			my %temp_hash=%{$GRtoGF_match_hash{$node}};
			my @keys=keys %temp_hash;
			if(scalar(@keys) > 0){
				$GRtoGF_map++;
				if(scalar(@keys) == 1){
					my %temp_hash2=%{$GFtoGR_match_hash{$keys[0]}};
					my @keys2=keys %temp_hash2;
					if(scalar(@keys2) == 1 && $node == $keys2[0]){
						$GRtoGF_match++;
					}
				}
			}
		}
		
		#Solo per stampa a video
		my $GFtoGR_map=0;
		foreach my $node(keys %GFtoGR_match_hash){			
			my %temp_hash=%{$GFtoGR_match_hash{$node}};
			my @keys=keys %temp_hash;
			if(scalar(@keys) > 0){
				$GFtoGR_map++;
			}
		}
		
		#Stampa a video
		print scalar(@GR_node_list), "\t", scalar(@GF_node_list), "\t", $GRtoGF_match, "\t", $GRtoGF_map, "\t", $GFtoGR_map, "\t";
		
		my %completely_new_edges=();
		
		#GR edges validation
		print OUT "GR edge validation\n";
		my $GR_match_edges=0;	#Biunivoci
		my $GR_map_edges=0;
		my $GR_edge_num=0;
		for(my $i=0; $i < $count_GR_node; $i++){
			for(my $j=0; $j < $count_GR_node; $j++){
				if($GR_adj_matx->[$i][$j] == 1){
					print OUT "\tGR edge (", $i+1, ",", $j+1, ") =>";
					my $GR_edge="(".($i+1).",".($j+1).")";
					$completely_new_edges{$GR_edge}="";
							
					$GR_edge_num++;
					
					my %left_match_nodes=%{$GRtoGF_match_hash{$i}};
					my @left_keys=keys %left_match_nodes;
					my %right_match_nodes=%{$GRtoGF_match_hash{$j}};
					my @right_keys=keys %right_match_nodes;
					
					my $map=0;
					for(my $k=0; $k <= $#left_keys; $k++){
						for(my $p=0; $p <= $#right_keys; $p++){
							if($GF_adj_matx->[$left_keys[$k]][$right_keys[$p]] == 1){
								$map=1;
							}
						}
					}
					$GR_map_edges+=$map;
					
					#Conteggio dei matching biunivoci tra gli edge in GR e quelli GF solo per stampa a video
					if(scalar(@left_keys) == 1 && scalar(@right_keys) == 1){
						my %left_temp_hash=%{$GFtoGR_match_hash{$left_keys[0]}};
						my %right_temp_hash=%{$GFtoGR_match_hash{$right_keys[0]}};
						my @left_keys2=keys %left_temp_hash;
						my @right_keys2=keys %right_temp_hash;
						if((scalar(@left_keys2) == 1 && $i == $left_keys2[0]) && (scalar(@right_keys2) == 1 && $j == $right_keys2[0])){
							if($GF_adj_matx->[$left_keys[0]][$right_keys[0]] == 1){
								$GR_match_edges++;
							}
						}
					}
												
					my %start_mapping_GF_nodes_h=();
					my %end_mapping_GF_nodes_h=();
		
					foreach my $key(keys %{${$GR_node_list[$i]}[1]}){
						$start_mapping_GF_nodes_h{$key}=1;
					}
					
					foreach my $key(keys %{${$GR_node_list[$j]}[1]}){
						$end_mapping_GF_nodes_h{$key}=1;
					}

					foreach my $key(keys %{${$GR_node_list[$i]}[3]}){
						$start_mapping_GF_nodes_h{$key}=1;
					}
					foreach my $key(keys %{${$GR_node_list[$j]}[3]}){
						$end_mapping_GF_nodes_h{$key}=1;
					}
					
					
					my @start_mapping_GF_nodes=keys %start_mapping_GF_nodes_h;
					my @end_mapping_GF_nodes=keys %end_mapping_GF_nodes_h;
					
					for(my $k=0; $k<=$#start_mapping_GF_nodes; $k++){
						for(my $p=0; $p<=$#end_mapping_GF_nodes; $p++){
							my $sGF_node=$start_mapping_GF_nodes[$k];
							my $eGF_node=$end_mapping_GF_nodes[$p];
							if($sGF_node == $eGF_node){
								print OUT " contr(", $sGF_node+1, ")";
								
							}
							else{
								if($GF_adj_matx->[$sGF_node][$eGF_node] == 1){
									print OUT " val(", $sGF_node+1, ",", $eGF_node+1, ") ";
										
								}
								else{
									my $new_GF_edge="(".($sGF_node+1).",".($eGF_node+1).")";
									#$new_GF_edge does not exists in the GF graph and is labelled as new
									$completely_new_edges{$GR_edge}=$new_GF_edge;
								}
							}
						}
					}
					print OUT "\n";
				}
			}
		}
		
		#GF edges validation
		print OUT "GF edge validation\n";
		my $GF_match_edges=0;
		my $GF_map_edges=0;
		my $GF_edge_num=0;
		for(my $i=0; $i < $count_GF_node; $i++){
			for(my $j=0; $j < $count_GF_node; $j++){
				if($GF_adj_matx->[$i][$j] == 1){
					print OUT "\tGF edge (", $i+1, ",", $j+1, ") =>";
					
					$GF_edge_num++;

					my %left_match_nodes=%{$GFtoGR_match_hash{$i}};
					my @left_keys=keys %left_match_nodes;
					my %right_match_nodes=%{$GFtoGR_match_hash{$j}};
					my @right_keys=keys %right_match_nodes;
					
					my $map=0;
					for(my $k=0; $k <= $#left_keys; $k++){
						for(my $p=0; $p <= $#right_keys; $p++){
							if($GR_adj_matx->[$left_keys[$k]][$right_keys[$p]] == 1){
								$map=1;
							}
						}
					}
					$GF_map_edges+=$map;
					
					my %start_mapping_GR_nodes_h=();
					my %end_mapping_GR_nodes_h=();
		
					foreach my $key(keys %{${$GF_node_list[$i]}[1]}){
						$start_mapping_GR_nodes_h{$key}=1;
					}
					foreach my $key(keys %{${$GF_node_list[$j]}[1]}){
						$end_mapping_GR_nodes_h{$key}=1;
					}
					
					foreach my $key(keys %{${$GF_node_list[$i]}[3]}){
						$start_mapping_GR_nodes_h{$key}=1;
					}
					foreach my $key(keys %{${$GF_node_list[$j]}[3]}){
						$end_mapping_GR_nodes_h{$key}=1;
					}
		
					my @start_mapping_GR_nodes=keys %start_mapping_GR_nodes_h;
					my @end_mapping_GR_nodes=keys %end_mapping_GR_nodes_h;
					
					for(my $k=0; $k<=$#start_mapping_GR_nodes; $k++){
						for(my $p=0; $p<=$#end_mapping_GR_nodes; $p++){
							my $sGR_node=$start_mapping_GR_nodes[$k];
							my $eGR_node=$end_mapping_GR_nodes[$p];
							if($sGR_node == $eGR_node){
								print OUT " contr(", $sGR_node+1, ")";
							}
							else{
								if($GR_adj_matx->[$sGR_node][$eGR_node] == 1){
									print OUT " val(", $sGR_node+1, ",", $eGR_node+1, ") ";
								}
							}
						}
					}
					print OUT "\n";
				}
			}
		}
                #Print Gnuplot script file
                print GNUPLOT "set term post dash size 9,6 \"Verdana\"\n";
                print GNUPLOT "set title \"Gene - $gene\"\n";
		print GNUPLOT "set key off\n";
		print GNUPLOT "unset ytics\n";
		my $ps_file = $gnuplot_file;
		$ps_file =~ s/\.[^\.]*$//;
                print GNUPLOT "set output '$ps_file.ps'\n\n";
                print GNUPLOT "set style line 1 lt rgb \"black\" lw 3\n";
                print GNUPLOT "set style line 2 lt rgb \"red\" lw 3\n\n";
                my $tot_len = 0;
                foreach my $q(0..$#GF_node_list){
                  my @list=@{$GF_node_list[$q]};
                  $tot_len = $tot_len + length($list[0]);
                }
                my $gnuplot_objs = 0;
                
		print OUT "#GF_BLOCKS\n";
                my $last_inserted = 0;
                my $space_len= POSIX::ceil($tot_len/100*2);
                my @c_begin = ();
                my @c_end = ();
		foreach my $q(0..$#GF_node_list){
                  my @list=@{$GF_node_list[$q]};
                  print OUT $q+1, "(", length($list[0]), ") ";
                  my $this_end= $last_inserted + length($list[0]);
                  $gnuplot_objs = $gnuplot_objs + 1;
                  print GNUPLOT "set object $gnuplot_objs rect from $last_inserted,0 to $this_end,1 fc lt 3 lw 1 front\n";
                  push @c_begin, $last_inserted;
                  push @c_end, $this_end;
                  $last_inserted = $this_end + $space_len;
                }
                my $max_intron = 0;
                for(my $i=0; $i < $count_GF_node; $i++){
                  for(my $j=0; $j < $count_GF_node; $j++){
                    if($GF_adj_matx->[$i][$j] == 1){
                      my $diff = ($c_begin[$j] - $c_end[$i]);
                      if($diff > $max_intron){
                        $max_intron = $diff;
                      }
                    }
                  }
                }
                my $max_h= 0;
                for(my $i=0; $i < $count_GF_node; $i++){
                  for(my $j=0; $j < $count_GF_node; $j++){
                    if($GF_adj_matx->[$i][$j] == 1){
                      my $diff = ($c_begin[$j] - $c_end[$i])/2;
                      my $height = ($diff / $max_intron)*10 + 1;
                      print GNUPLOT "set arrow from $c_end[$i],1 to ", $c_end[$i]+$diff,",$height nohead ls 1\n";
                      print GNUPLOT "set arrow from ",$c_end[$i]+$diff,",$height to $c_begin[$j],1 nohead ls 1\n";
                      if($height > $max_h){
                        $max_h = $height;
                      }
                    }
                  }
                }
		print OUT "\n#GR_BLOCKS\n";
		foreach my $q(0..$#GR_node_list){
			my @list=@{$GR_node_list[$q]};
			print OUT $q+1, "(", length($list[0]), ") ";
		}
		print OUT "\n";
		my $max_new_h = 0;
		print_node_variations(\@GF_node_list, 1, $GR_adj_matx, \*OUT, $gene, $count_GR_node, \@GR_node_list, \@c_begin, \@c_end, $gnuplot_objs, \$last_inserted, $max_intron, \$max_new_h);
		
		print OUT "#NEW_LINKS\n";	#In GR wrt GF
                
		my @compl_new_edges=keys %completely_new_edges;
		foreach my $e(@compl_new_edges){
                  if($completely_new_edges{$e} ne ""){
                    print OUT $completely_new_edges{$e}, " ";
                    $completely_new_edges{$e} =~ m/\(([0-9]+),([0-9]+)\)/;
                    my $start= ($1)-1;
                    my $end= ($2)-1;
                    if($start <= $#GF_node_list && $end <= $#GF_node_list){
                      if($start > $end){
                        my $tmp = $start;
                        $start = $end;
                        $end = $tmp;
                      }
                      my $diff = ($c_begin[$end] - $c_end[$start])/2;
                      my $height = ($diff / $max_intron)*10 + 1;
                      print GNUPLOT "set arrow from $c_end[$start],1 to ", $c_end[$start]+$diff,",$height nohead ls 2\n";
                      print GNUPLOT "set arrow from ",$c_end[$start]+$diff,",$height to $c_begin[$end],1 nohead ls 2\n";
                      if($height > $max_h){
                        $max_h = $height;
                      }
                    }
                  }
		}
		print OUT "\n";
		#print OUT join(",", @compl_new_edges), "\n";
		
		print $GR_edge_num, "\t", $GF_edge_num, "\t", $GR_match_edges, "\t", $GR_map_edges, "\t", $GF_map_edges, "\t";
		print "\n";

                print GNUPLOT "plot [-1:$last_inserted][-",$max_new_h+1,":",$max_h+1,"] 0 lt rgb \"blue\" lw 3\n";
					
		close OUT;
		close LOG;
	}
}	
exit;


#SUBROUTINES

sub get_node_sequence {
    my $left_right_ref=shift;
    my $gen_seq=shift;
    
    my @region_left_right=@{$left_right_ref};
    die "Wrong node regions!\n" unless (scalar(@region_left_right) % 2 == 0);
    
    my $node_seq="";
    my $i=0;
    while($i < $#region_left_right){
    	$node_seq.=substr($gen_seq, $region_left_right[$i]-1, $region_left_right[$i+1]-$region_left_right[$i]+1);
		$i+=2;
    }
        
    return $node_seq;
}

sub print_node_mapping {
    my $node_list_ref=shift;
    my @node_list=@{$node_list_ref};
    my $is_GF_node=shift;
    local *FH=shift;
    my $gene=shift;
    
    my %match_nodes=();
            
 	for(my $i=0; $i <= $#node_list; $i++){
		my @node=@{$node_list[$i]};
		my $node_seq="";
		$node_seq=$node[0];
		
		my %mapping_hash=%{$node[1]};
		my %mapping_cut_hash=%{$node[2]};
		my %mapping_contained_hash=%{$node[3]};
		
		if($is_GF_node){
			print FH "\tGF node ";
		}
		else{
			print FH "\tGR node ";
		}
		print FH $i+1, " =>";
		
		my %match_hash=();
			
		if(scalar(keys %mapping_hash) > 0){

			foreach my $key(keys %mapping_hash){
				my %hash_temp=%{$mapping_hash{$key}};
				foreach my $offs(keys %hash_temp){
					my @off_list=@{$hash_temp{$offs}};
					if(scalar(@off_list) == 2 || $off_list[2] == 1){
						print FH " map(";
						print FH $key+1;
						print FH "[";
						print FH $offs, ",", $off_list[0], ",", $off_list[1];
						print FH "]";
						print FH ")";
					
						if($offs <= $mapping_threshold && $off_list[1] <= $mapping_threshold){
							$match_hash{$key}=1 if(scalar(@off_list) == 2 || $off_list[2] == 1);
						}
					}	
				}
			}
		}
		
		if(scalar(keys %mapping_cut_hash) > 0){
			foreach my $key(keys %mapping_cut_hash){
				my %hash_temp=%{$mapping_cut_hash{$key}};
				foreach my $next_node(keys %hash_temp){
					my @cuts=@{$hash_temp{$next_node}};
					my $first_cut=$key+1;
					if($first_cut == 0){
						$first_cut="?";
					}
					my $second_cut=$next_node+1;
					if($second_cut == 0){
						$second_cut="?";
					}
					print FH " cut(";
					print FH "\{", $first_cut, "[", $cuts[0], ",", $cuts[1], "],", $second_cut, "[", $cuts[2], ",", $cuts[3], "]", "\}";
					print FH ")";
				}
			}
		}
		
		if(scalar(keys %mapping_contained_hash) > 0){
			foreach my $key(keys %mapping_contained_hash){
				my %hash_temp=%{$mapping_contained_hash{$key}};
				foreach my $offs(keys %hash_temp){
					my @off_list=@{$hash_temp{$offs}};
					if(scalar(@off_list) == 2 || $off_list[2] == 1){
						print FH " contain(";
						print FH $key+1;
						print FH "[";
						print FH $offs, ",", $off_list[0], ",", $off_list[1];
						print FH "]";
						print FH ")";
					
						if($offs <= $mapping_threshold && $off_list[1] <= $mapping_threshold){
							$match_hash{$key}=1;
						}
					}
				}
			}
		}
		
		$match_nodes{$i}=\%match_hash;
		
		print FH "\n";
	}
	
	#Un nodo mappa a un altro se e solo se le stringhe sono $mapping_threshold-trim equivalenti
	return \%match_nodes;
}

sub handle_small_nodes {
    my $node_list_ref=shift;
    my @node_list=@{$node_list_ref};
    my $small_length=shift;
    my $is_GF=shift;
    local *FH=shift;
          
 	for(my $i=0; $i <= $#node_list; $i++){
		my @node=@{$node_list[$i]};
		my $node_seq=$node[0];

		#DA METTERE COME PARAMETRO??
		if(length($node_seq) < $small_length){
			if($is_GF){
				print FH "GF#";
			}else{
				print FH "GR#";
			}
			print FH ($i+1), " is a small node\n";

			my %mapping_hash=%{$node[1]};
			my %mapping_cut_hash=%{$node[2]};
			my %mapping_contained_hash=%{$node[3]};

			my $key_min_length=-1;
			my $min_length=-1;
			foreach my $key(keys %mapping_hash){
				
				my %hash_temp=%{$mapping_hash{$key}};
				my @all_keys=keys %hash_temp;
				my $first_key=$all_keys[0];
				my @off_list=@{$hash_temp{$first_key}};
				my $length=$first_key+$off_list[0]+$off_list[1];
				if($key_min_length == -1 || $length < $min_length){
					$key_min_length=$key;
					$min_length=$length;
				}
			}
			my $all=1;
			
			foreach my $key(keys %mapping_hash){
				
				if($key != $key_min_length){
					delete ($mapping_hash{$key});
				}else{
					
					$all=0;
				}
			}
			if($all == 1){
				%mapping_hash=();
			}

			$node[1]=\%mapping_hash;

			my $key_max_length=-1;
			my $max_length=-1;
			foreach my $key(keys %mapping_contained_hash){
				my %hash_temp=%{$mapping_contained_hash{$key}};
				my @all_keys=keys %hash_temp;
				my $first_key=$all_keys[0];
				my @off_list=@{$hash_temp{$first_key}};
				my $length=$off_list[0];
				if($key_max_length == -1 || $length > $max_length){
					$key_max_length=$key;
					$max_length=$length;
				}
			}
			$all=1;
			
			foreach my $key(keys %mapping_contained_hash){
				if($key != $key_max_length){
					delete ($mapping_contained_hash{$key});
				}else{
					$all=0;
				}
			}
			if($all == 1){
				%mapping_contained_hash=();
			}

			$node[3]=\%mapping_contained_hash;
			$node[2]=\%mapping_cut_hash;
		}
		$node_list[$i]=\@node;
	}
	
	return \@node_list;
}

sub mapremove_micro_nodes {
    my $node_list_ref=shift;
    my @node_list=@{$node_list_ref};
    my $micro_length=shift;
    my $is_GF=shift;
    local *FH=shift;
            
 	for(my $i=0; $i <= $#node_list; $i++){
		my @node=@{$node_list[$i]};
		my $node_seq=$node[0];

		my %mapping_hash=%{$node[1]};
		my %mapping_cut_hash=%{$node[2]};
		my %mapping_contained_hash=%{$node[3]};

		if(length($node_seq) < $micro_length){
			if($is_GF){
				print FH "GF#";
			}else{
				print FH "GR#";
			}
			print FH ($i+1), " is a micro node and has been removed from the mapping\n";
			%mapping_hash=();
			%mapping_cut_hash=();
			%mapping_contained_hash=();
		}else{
			my $all=1;
			foreach my $key(keys %mapping_contained_hash){
				my %hash_temp=%{$mapping_contained_hash{$key}};
				my @all_keys=keys %hash_temp;
				my $first_key=$all_keys[0];
				my @off_list=@{$hash_temp{$first_key}};
				my $length=$off_list[0];
				if($length < $micro_length){
					if($is_GF){
						print FH "GR#";
					}else{
						print FH "GF#";
					}
					print FH $key, " is a micro node and has been removed from the containment of ";
					if($is_GF){
						print FH "GF#", ($i+1), "\n";
					}else{
						print FH "GR#", ($i+1), "\n";
					}
					delete ($mapping_contained_hash{$key});
				}
				else{
					$all=0;
				}
			}
			if($all == 1){
				%mapping_contained_hash=();
			}
		}
		$node[1]=\%mapping_hash;
		$node[2]=\%mapping_cut_hash;
		$node[3]=\%mapping_contained_hash;

		$node_list[$i]=\@node;
	}
	return \@node_list;
}


sub print_node_variations {
    my $node_list_ref=shift;
    my @node_list=@{$node_list_ref};
    my $is_GF_node=shift;
    my $adj_matx_to=shift; #if $is_GF_node is true, this is the GR matrix (and viceversa)
    local *FH=shift;
    my $gene=shift;
    my $count_node=shift;	#if $is_GF_node is true, this is the number of nodes in GR (and viceversa)
    my $node_list_ref2=shift; #Node list of GR if $is_GF_node is true (and viceversa)
    my $start_GF=shift;
    my $end_GF=shift;
    my $n_obj=shift;
    my $last_inserted=shift;
    my $max_introns=shift;
    my $max_new_h=shift;
    my @node_list2=@{$node_list_ref2};
    
    #for $is_GF_node=1 the new nodes are in GR (and viceversa)
    my %completely_new_nodes=();
    foreach my $n(1..$count_node){
    	$completely_new_nodes{$n}=1;
    }
    
    if($is_GF_node){
    	print FH "#NEW_GR_BLOCKS\n";
    }else{
    	print FH "#NEW_GF_BLOCKS\n\n";
    }
    
 	for(my $i=0; $i <= $#node_list; $i++){
		my @node=@{$node_list[$i]};
		
		my $length=length($node[0]);
		
		my %mapping_hash=%{$node[1]};
		my %mapping_contained_hash=%{$node[3]};
		
		my @variant_node_list=();
		my %variant_node_offs=();
			
		if(scalar(keys %mapping_hash) > 0){
			foreach my $key(keys %mapping_hash){
				my %hash_temp=%{$mapping_hash{$key}};
				foreach my $offs(keys %hash_temp){
					my @off_list=@{$hash_temp{$offs}};
					$offs=-$offs;
					if(exists $variant_node_offs{$offs}){
						push @{$variant_node_offs{$offs}}, $key+1;
					}else{
						$variant_node_offs{$offs}=[$key+1];
					}
				}
								
				push @variant_node_list, ($key+1);
				$completely_new_nodes{$key+1}=0;
			}
		}
		
		if(scalar(keys %mapping_contained_hash) > 0){
			foreach my $key(keys %mapping_contained_hash){
				#MODIFICARE
				my $found=0;
				foreach my $n(@variant_node_list){
					if($n == $key+1){
						$found=1;
					}
				}

				if($found == 0){
					my %hash_temp=%{$mapping_contained_hash{$key}};
					foreach my $offs(keys %hash_temp){
						my @off_list=@{$hash_temp{$offs}};
						if(exists $variant_node_offs{$offs}){
							push @{$variant_node_offs{$offs}}, $key+1;
						}else{
							$variant_node_offs{$offs}=[$key+1];
						}
					}
				}				
				
				push @variant_node_list, ($key+1);
				$completely_new_nodes{$key+1}=0;
			}
		}

		my %new_blocks=();
		my %new_edges=();
		my $p=0;
		while($p < $#variant_node_list){
			my $j=$p+1;
			
			while($j <= $#variant_node_list){
				my $node1=$variant_node_list[$p];
				my $node2=$variant_node_list[$j];
				
				if($adj_matx_to->[$node1-1][$node2-1] == 1){
					my $new_edge="(".$node1.",".$node2.")";
					$new_blocks{$node1}=1;
					$new_blocks{$node2}=1;
					$new_edges{$new_edge}=1;
				}else{
					if($adj_matx_to->[$node2-1][$node1-1] == 1){
						my $new_edge="(".$node2.",".$node1.")";
						$new_blocks{$node1}=1;
						$new_blocks{$node2}=1;
						$new_edges{$new_edge}=1;
					}
				}
				$j++;
			}
			$p++;
		}
                
		my @offsets=sort {$a <=> $b} keys %variant_node_offs;
		#my @new_blocks=();
		#foreach my $off(@offsets){
		#	push @new_blocks, @{$variant_node_offs{$off}};
		#}
		
		#if(scalar(@new_blocks) > 0){
			if($is_GF_node){
				print FH "GF#";
			}
			else{
				print FH "GR#";
			}
			print FH $i+1, "(", $length, ")\t=>\t";
			
			my @new_edges=keys %new_edges;
			my @start_GR;
			my @end_GR;
			foreach my $off(@offsets){
				my @node_index=@{$variant_node_offs{$off}};
				foreach my $n(@node_index){
					my $index_n=$n-1;
					my @node_r=@{$node_list2[$index_n]};		
					my $l=length($node_r[0]);
					if($is_GF_node){
						print FH "GR";
					}
					else{
						print FH "GF";
					}
					print FH "#", $n, "[", ($off+1), ",", ($off+$l), "] ";
                                        $n_obj = $n_obj + 1;
                                        print GNUPLOT "set object $n_obj rect from ", ($start_GF->[$i])+($off+1),",-1 to ",($start_GF->[$i])+($off+$l),",0 fc lt 2 lw 1 front\n";
					$start_GR[$n] = ($start_GF->[$i])+($off+1);
					$end_GR[$n] = ($start_GF->[$i])+($off+$l);
                                        if(($start_GF->[$i])+($off+1) > $$last_inserted){
                                          $$last_inserted = ($start_GF->[$i])+($off+1);
                                        }
				}
			}
			
			if(scalar(@new_edges) > 0){
				if($is_GF_node){
					print FH "GR";
				}
				else{
					print FH "GF";
				}
				print FH "links(", join(";", @new_edges), ")";
				foreach my $e(@new_edges){
					$e =~ m/\(([0-9]+),([0-9]+)\)/;
                    			my $start= ($1);
                    			my $end= ($2);
					#print "Nodi GR Start $1 end $2\n";
                      			if($start > $end){
                        			my $tmp = $start;
                        			$start = $end;
                        			$end = $tmp;
                      			}
                      			my $diff = ($start_GR[$end] - $end_GR[$start])/2;
                      			my $height = ($diff / $max_introns)*10 + 1;
					$height=2 if $height<2;
					print GNUPLOT "set arrow from $end_GR[$start],-1 to ", $end_GR[$start]+$diff,",-$height nohead ls 2\n";
	                        	print GNUPLOT "set arrow from ",$end_GR[$start]+$diff,",-$height to $start_GR[$end],-1 nohead ls 2\n";
					if($height > $$max_new_h){
						$$max_new_h = $height;
					}
				}
			}

			print FH "\n";
			
				
		#}
	}	
	
	if($is_GF_node){
    	print FH "#UNMAPPED_GR_BLOCKS\n";
    }else{
    	print FH "#UNMAPPED_GF_BLOCKS\n";
    }
    
    my @news=sort {$a <=> $b} keys %completely_new_nodes;
    foreach my $new(@news){
    	if($completely_new_nodes{$new} == 1){
    		my $index_n=$new-1;
			my @node_r=@{$node_list2[$index_n]};		
			my $l=length($node_r[0]);
    		print FH $new, "(", $l, ")", ",";
    	}    	
    }
    print FH "\n";
  }

sub get_fingerprinting {
    my $node_list_ref=shift;
    my @node_list=@{$node_list_ref};
    my $f_length=shift;
    my $f_offset=shift;
    my $debug=shift;
	
	my %left_fingerp=();
	my %right_fingerp=();
	my $counter=0;

	foreach my $node_ref(@node_list){
		my @node=@{$node_ref};
		my $node_seq=$node[0];
	
		#Left fingerprint
		my $start_offset=0;
		while(($f_offset == -1 || $start_offset <= $f_offset) && $start_offset <= length($node_seq)-$f_length){
			my $l_fingerp=substr($node_seq, $start_offset, $f_length);
		
			my @finger_node=($counter, $start_offset);
			my @finger_list;
			if(exists $left_fingerp{$l_fingerp}){
				@finger_list=@{$left_fingerp{$l_fingerp}};
			}
			else{
				@finger_list=();
			}
			push @finger_list, \@finger_node;
			$left_fingerp{$l_fingerp}=\@finger_list;
			$start_offset++;
		}
	
		#Right fingerprint
		#if($f_offset != -1){
			my $end_offset=length($node_seq)-1;
			while(($f_offset == -1 || $end_offset >= length($node_seq)-$f_offset-1) && $end_offset >= $f_length-1){
				my $r_fingerp=substr($node_seq, $end_offset-$f_length+1, $f_length);

				my @finger_node=($counter, $end_offset);
				my @finger_list;
				if(exists $right_fingerp{$r_fingerp}){
					@finger_list=@{$right_fingerp{$r_fingerp}};
				}
				else{
					@finger_list=();
				}
				push @finger_list, \@finger_node;
				$right_fingerp{$r_fingerp}=\@finger_list;
				$end_offset--;
			}
		#}
	
		$counter++;
	}

	return [\%left_fingerp, \%right_fingerp];
	
}

sub print_node_list {
	my $node_list_ref=shift;
    my @node_list=@{$node_list_ref};
    my $is_GF_node=shift;
    local *FH=shift;
    
	my $count_n=1;
	foreach my $node_ref(@node_list){
		my @node=@{$node_ref};
		my $node_seq=$node[0];
		if($is_GF_node){
			print FH "GF_GRAPH";			
		}
		else{
			print FH "GR_GRAPH";
			
		}
		print FH " node\#", $count_n;
		print FH " sequence(", $node_seq, ")\n";
		$count_n++;
	}
}

sub print_adj_matrix {
	my $adj_matx=shift;
	my $count_node=shift;
    my $is_GF_node=shift;
    local *FH=shift;
    
    my $tag;        
	if($is_GF_node){
		$tag="GF_ADJ";
	}
	else{
		$tag="GR_ADJ";
	}
	
	for(my $i=0; $i < $count_node; $i++){
		print FH "$tag ";
		for(my $j=0; $j < $count_node; $j++){
			print FH $adj_matx->[$i][$j], " ";
		}
		print FH "\n";
	}
}

sub print_fingerprinting {
	my $fingerp_ref=shift;
	my %fingerp=%{$fingerp_ref};
	my $is_GF_node=shift;
	local *FH=shift;
	
	foreach my $key(keys %fingerp){
		if($is_GF_node){
			print FH "GF_FNPR ", $key." in(";
		}
		else{
			print FH "GR_FNPR ", $key." in(";
		}
		my @finger_list=@{$fingerp{$key}};
		my $comma=0;
		foreach my $finger_node_ref(@finger_list){
			if($comma){
				print FH ";";
			}
			$comma=1;
			my @finger_node=@{$finger_node_ref};
			print FH "node(".($finger_node[0]+1). "),offset(".$finger_node[1].")";
		}
		print FH ")\n";
	}
}

sub map_nodes {
	my $node_tbm_list_ref=shift;
	my $node_to_list_ref=shift;
    my @node_tbm_list=@{$node_tbm_list_ref};
    my @node_to_list=@{$node_to_list_ref};
    my $left_fingerp_to_ref=shift;
    my $right_fingerp_to_ref=shift;
    my %left_fingerp_to=%{$left_fingerp_to_ref};
    my %right_fingerp_to=%{$right_fingerp_to_ref};
    my $adj_matx_to=shift;
    my $f_offset=shift;
    my $debug=shift;
	local *FH=shift;
	
	for(my $i=0; $i <= $#node_tbm_list; $i++){
		my @node_tbm=@{$node_tbm_list[$i]};
		my $node_tbm_seq=$node_tbm[0];
		
		my $node_tbm_prefix=substr($node_tbm_seq, 0, $f_length);
		my $node_tbm_suffix=substr($node_tbm_seq, length($node_tbm_seq)-$f_length, $f_length);
	
		my %mapping_hash=();
		push @{${$node_tbm_list_ref}[$i]}, \%mapping_hash;
	
		my %mapping_cut_hash=();
		push @{${$node_tbm_list_ref}[$i]}, \%mapping_cut_hash;
	
		if(exists $left_fingerp_to{$node_tbm_prefix}){
			my @finger_list=@{$left_fingerp_to{$node_tbm_prefix}};
			foreach my $finger_node_ref(@finger_list){
			
				my @finger_node=@{$finger_node_ref};
				my @node_to=@{$node_to_list[$finger_node[0]]};
				my $node_to_seq=$node_to[0];
							
				if(length($node_to_seq)-$finger_node[1] >= length($node_tbm_seq)){
					my $sub_to=substr($node_to_seq, $finger_node[1], length($node_tbm_seq));
				
					if($node_tbm_seq eq $sub_to){
								
						my %hash_temp=%{${$node_tbm_list[$i]}[1]};
						my %hash_temp2;

						if(exists $hash_temp{$finger_node[0]}){
							%hash_temp2=%{$hash_temp{$finger_node[0]}};
						}
						else{
							%hash_temp2=();
						}
					
						my $cut_suffix=length($node_to_seq)-$finger_node[1]-length($node_tbm_seq);
						$hash_temp2{$finger_node[1]}=[length($node_tbm_seq), $cut_suffix];
						$hash_temp{$finger_node[0]}=\%hash_temp2;
						
						${${$node_tbm_list_ref}[$i]}[1]=\%hash_temp;
					}
				}
				else{
					#Controllo inutile
					if(length($node_to_seq)-$finger_node[1] > 0){
										
						my $suff_to=substr($node_to_seq, $finger_node[1], length($node_to_seq)-$finger_node[1]);
						my $pref_tbm=substr($node_tbm_seq, 0, length($suff_to));
						my $suff_length=length($node_tbm_seq)-length($pref_tbm);
						if($suff_to eq $pref_tbm){
							my $suff_tbm=substr($node_tbm_seq, length($node_tbm_seq)-$suff_length, $suff_length);
							my $cut_ok=0;
							for(my $q=0; $q <= $#node_to_list; $q++){
								if($q != $finger_node[0] && $adj_matx_to->[$finger_node[0]][$q]){
									my @next_node_to=@{$node_to_list[$q]};
									my $next_node_to_seq=$next_node_to[0];
									
									if($suff_length <= length($next_node_to_seq)){
								
										my $next_pref_to=substr($next_node_to_seq, 0, $suff_length);

										if($suff_tbm eq $next_pref_to && length($suff_tbm) > $mapping_threshold){
										#if($suff_tbm eq $next_pref_to){
											$cut_ok=1;
											my %hash_temp=%{${$node_tbm_list[$i]}[2]};
											my %hash_temp2;
											if(exists $hash_temp{$finger_node[0]}){
												%hash_temp2=%{$hash_temp{$finger_node[0]}};
											}
											else{
												%hash_temp2=();
											}
											$hash_temp2{$q}=[length($node_to_seq)-length($pref_tbm), length($pref_tbm), length($suff_tbm), length($next_node_to_seq)-length($suff_tbm)];
											$hash_temp{$finger_node[0]}=\%hash_temp2;
											${${$node_tbm_list_ref}[$i]}[2]=\%hash_temp;
										}
									}
								}
							}
							
							if($cut_ok == 0  && length($suff_tbm) <= $mapping_threshold){
#								my %hash_temp=%{${$node_tbm_list[$i]}[2]};
#								my %hash_temp2;
#								if(exists $hash_temp{$finger_node[0]}){
#									%hash_temp2=%{$hash_temp{$finger_node[0]}};
#								}
#								else{
#									%hash_temp2=();
#								}
#								my $false_node=-1;
#								$hash_temp2{$false_node}=[length($node_to_seq)-length($pref_tbm), length($pref_tbm), length($suff_tbm), -1];
#								$hash_temp{$finger_node[0]}=\%hash_temp2;
#								${${$node_tbm_list_ref}[$i]}[2]=\%hash_temp;

								my %hash_temp=%{${$node_tbm_list[$i]}[1]};
								my %hash_temp2;
								if(exists $hash_temp{$finger_node[0]}){
									%hash_temp2=%{$hash_temp{$finger_node[0]}};
								}
								else{
									%hash_temp2=();
								}
												
								$hash_temp2{length($node_to_seq)-length($pref_tbm)}=[length($pref_tbm), -length($suff_tbm)];
								$hash_temp{$finger_node[0]}=\%hash_temp2;
								${${$node_tbm_list_ref}[$i]}[1]=\%hash_temp;									
							}
						}
					}
				}
			}
		}
	
		#if($f_offset != -1){
			if(exists $right_fingerp_to{$node_tbm_suffix}){
				my @finger_list=@{$right_fingerp_to{$node_tbm_suffix}};
				foreach my $finger_node_ref(@finger_list){
					my @finger_node=@{$finger_node_ref};
					my @node_to=@{$node_to_list[$finger_node[0]]};
					my $node_to_seq=$node_to[0];
					
					#Inutile per $f_offset = -1
					if($finger_node[1]-length($node_tbm_seq)+1 >= 0){
						
						my $sub_to=substr($node_to_seq, $finger_node[1]-length($node_tbm_seq)+1, length($node_tbm_seq));
							
						if($node_tbm_seq eq $sub_to){
							
							my %hash_temp=%{${$node_tbm_list[$i]}[1]};
							my %hash_temp2;
							if(exists $hash_temp{$finger_node[0]}){
								%hash_temp2=%{$hash_temp{$finger_node[0]}};
							}
							else{
								%hash_temp2=();
							}
							my $cut_suffix=length($node_to_seq)-$finger_node[1]-1;
							$hash_temp2{$finger_node[1]-length($node_tbm_seq)+1}=[length($node_tbm_seq), $cut_suffix];
							$hash_temp{$finger_node[0]}=\%hash_temp2;
							${${$node_tbm_list_ref}[$i]}[1]=\%hash_temp;
						}
					}
					else{
						#Controllo inutile
						if($finger_node[1] >= 0){
							my $pref_to=substr($node_to_seq, 0, $finger_node[1]+1);
							my $suff_tbm=substr($node_tbm_seq, length($node_tbm_seq)-length($pref_to), length($pref_to));
							
							my $pref_length=length($node_tbm_seq)-length($suff_tbm);
							
							if($pref_to eq $suff_tbm){
								my $pref_tbm=substr($node_tbm_seq, 0, $pref_length);
								my $cut_ok=0;
								for(my $q=0; $q <= $#node_to_list; $q++){
									if($q != $finger_node[0] && $adj_matx_to->[$q][$finger_node[0]]){
										my @prev_node_to=@{$node_to_list[$q]};
										my $prev_node_to_seq=$prev_node_to[0];
	
										if($pref_length <= length($prev_node_to_seq)){
											my $prev_suff_to=substr($prev_node_to_seq, length($prev_node_to_seq)-$pref_length, $pref_length);
											
											if($pref_tbm eq $prev_suff_to && length($pref_tbm) > $mapping_threshold){
											#if($pref_tbm eq $prev_suff_to){
												$cut_ok=1;	
												my %hash_temp=%{${$node_tbm_list[$i]}[2]};
												my %hash_temp2;
												if(exists $hash_temp{$q}){
													%hash_temp2=%{$hash_temp{$q}};
												}
												else{
													%hash_temp2=();
												}
												$hash_temp2{$finger_node[0]}=[length($prev_node_to_seq)-length($pref_tbm), length($pref_tbm), length($suff_tbm), length($node_to_seq)-length($suff_tbm)];
												$hash_temp{$q}=\%hash_temp2;
												${${$node_tbm_list_ref}[$i]}[2]=\%hash_temp;
											}
										}
									}
								}
								
								if($cut_ok == 0 && length($pref_tbm) <= $mapping_threshold){
#									my %hash_temp=%{${$node_tbm_list[$i]}[2]};
#									my %hash_temp2;
#									my $false_node=-1;
#									if(exists $hash_temp{$false_node}){
#										%hash_temp2=%{$hash_temp{$false_node}};
#									}
#									else{
#										%hash_temp2=();
#									}
#													
#									$hash_temp2{$finger_node[0]}=[-1, length($pref_tbm), length($suff_tbm), length($node_to_seq)-length($suff_tbm)];
#									$hash_temp{$false_node}=\%hash_temp2;
#									${${$node_tbm_list_ref}[$i]}[2]=\%hash_temp;
									
									my %hash_temp=%{${$node_tbm_list[$i]}[1]};
									my %hash_temp2;
									if(exists $hash_temp{$finger_node[0]}){
										%hash_temp2=%{$hash_temp{$finger_node[0]}};
									}
									else{
										%hash_temp2=();
									}
													
									$hash_temp2{-length($pref_tbm)}=[length($suff_tbm), length($node_to_seq)-length($suff_tbm)];
									$hash_temp{$finger_node[0]}=\%hash_temp2;
									${${$node_tbm_list_ref}[$i]}[1]=\%hash_temp;									
								}
							}
						}
					}
				}
			}
		#}
	}
}

sub update_with_contained_nodes {
	my $node_tbu_list_ref=shift;
	my $node_u_list_ref=shift;
	my @node_tbu_list=@{$node_tbu_list_ref};
    my @node_u_list=@{$node_u_list_ref};
    my $debug=shift;
    
    for(my $i=0; $i <= $#node_tbu_list; $i++){
		my %contained_hash=();
		push @{${$node_tbu_list_ref}[$i]}, \%contained_hash;
    }
    
    for(my $i=0; $i<=$#node_u_list; $i++){
    	my %hash_temp=%{${$node_u_list[$i]}[1]};
     	foreach my $map_node(keys %hash_temp){
    		my %contained_hash_temp=%{${$node_tbu_list[$map_node]}[scalar(@{${$node_tbu_list_ref}[$map_node]})-1]};
    		my %hash_temp2=%{$hash_temp{$map_node}};
    		foreach my $pref_length(keys %hash_temp2){
    			my @off_list=@{$hash_temp2{$pref_length}};
    			my $map_length=$off_list[0];
    			my $suff_length=$off_list[1];
    			my %contained_hash_temp2;
				if(exists $contained_hash_temp{$i}){
					%contained_hash_temp2=%{$contained_hash_temp{$i}};
				}
				else{
					%contained_hash_temp2=();
				}
					
				$contained_hash_temp2{$pref_length}=[$map_length, $suff_length];
				$contained_hash_temp{$i}=\%contained_hash_temp2;
				${${$node_tbu_list_ref}[$map_node]}[scalar(@{${$node_tbu_list_ref}[$map_node]})-1]=\%contained_hash_temp;
    		}
    	}
    }
}
